#!/usr/bin/env python3
"""Phase 21.2 -- Source-audit cross-walk.

Reads docs/parity/feature_inventory.tsv produced by Pass 21.1 and fills in
the `our_status` column for every row by indexing the entire src/ tree
once into a token -> [(subsystem, file, line)] map and matching each
feature's salient keywords against it.

Status assignment rules
-----------------------
  n/a       row is bookkeeping (Copyright/Trademark/Table of Contents/
            Prerequisites/Problem Description/Typographical Conventions/
            Using This Manual/section-number-only headings/etc.)
  present   >=3 distinct DIAGNOSTIC tokens match the src/ index
            (diagnostic = token has <=80 src/ hits; super-common tokens
             like "point", "data", "surface" carry no parity signal)
  partial   exactly 1 or 2 diagnostic tokens match
  absent    no diagnostic token matches

A 6th column `citations` is added with up to 3 "subsystem/file:line"
refs from the strongest-matching token (highest src/ frequency).  The
spec text said "populate our_status with file/line citations" but the
enum {present,partial,absent,n/a} is mutually exclusive with free-form
citations, so we keep the enum in `our_status` and add `citations` as
a sibling column.  Pass 21.3 (gap classification) consumes both.

Stopword strategy: aggressive.  CFD docs are full of generic words
("flow", "model", "value", "setup", ...) that match everything in
src/ and produce false positives.  STOPWORDS below kills the worst
offenders; the remaining tokens are domain-meaningful nouns/verbs.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "src"
TSV_IN = ROOT / "docs" / "parity" / "feature_inventory.tsv"
TSV_OUT = TSV_IN  # in-place rewrite

SRC_SUFFIXES = {".hpp", ".h", ".cpp", ".cc", ".cxx", ".cmake"}
SRC_NAMES = {"CMakeLists.txt"}

TOKEN_RE = re.compile(r"[a-z][a-z0-9_]{3,}")   # >=4 chars, must start alpha

# Stopwords: too generic to be diagnostic of any specific capability.
# Keep this tight -- killing real terms here causes false "absent".
STOPWORDS = {
    # generic English nouns/verbs
    "this", "that", "with", "from", "into", "your", "their", "have", "been",
    "will", "shall", "must", "more", "less", "than", "then", "also", "such",
    "some", "many", "most", "only", "each", "both", "other", "either",
    "between", "before", "after", "above", "below", "while", "where", "when",
    "what", "which", "these", "those", "here", "there", "about", "over",
    "under", "through", "during", "within", "without", "should", "would",
    "could", "given", "using", "used", "uses", "make", "makes", "made",
    "set", "sets", "get", "gets", "got", "see", "sees", "seen", "show",
    "shown", "case", "cases", "step", "steps", "task", "tasks",
    # doc structural noise
    "manual", "guide", "tutorial", "tutorials", "chapter", "section",
    "page", "pages", "appendix", "preface", "foreword", "contents",
    "introduction", "overview", "summary", "conclusion", "abstract",
    "table", "tables", "figure", "figures", "list", "lists", "index",
    "reference", "references", "bibliography", "glossary", "license",
    "copyright", "trademark", "trademarks", "disclaimer", "notice",
    "publisher", "author", "version", "release", "edition", "isbn",
    "convention", "conventions", "typographical", "prerequisites",
    "problem", "description", "setup", "solution", "result", "results",
    # vendor / product names
    "ansys", "fluent", "comsol", "openfoam", "workbench", "multiphysics",
    "simall", "beta",
    # ultra-generic CFD nouns that match too broadly in src/
    "flow", "fluid", "model", "models", "modeling", "modelling", "value",
    "values", "data", "files", "file", "name", "names", "type", "types",
    "method", "methods", "option", "options", "parameter", "parameters",
    "system", "systems", "domain", "domains", "field", "fields",
    "general", "default", "settings", "specify", "specified",
    "specifying", "create", "creating", "creation", "define", "defining",
    "definition", "select", "selecting", "selection", "configure",
    "configuration", "compute", "computed", "computation", "computing",
    "calculate", "calculation", "input", "output", "inputs", "outputs",
    "function", "functions", "level", "levels", "size", "sized",
    "simulation", "simulations", "user", "users",
}


def tokenise(text: str) -> set[str]:
    """Lowercase token set, stopwords dropped."""
    return {t for t in TOKEN_RE.findall(text.lower()) if t not in STOPWORDS}


# ---------------------------------------------------------------------------
# n/a heuristics
# ---------------------------------------------------------------------------
_SECTION_NUM_RE = re.compile(r"^\s*\d+(\.\d+)*\.?\s+")  # "3.2.1. Foo" -> "Foo"

_NA_FEATURE_PATTERNS = [
    re.compile(r"^copyright", re.I),
    re.compile(r"^trademark", re.I),
    re.compile(r"^license", re.I),
    re.compile(r"^contents$", re.I),
    re.compile(r"^table of contents", re.I),
    re.compile(r"^foreword$", re.I),
    re.compile(r"^preface$", re.I),
    re.compile(r"^introduction$", re.I),
    re.compile(r"^overview$", re.I),
    re.compile(r"^summary$", re.I),
    re.compile(r"^prerequisites?$", re.I),
    re.compile(r"^problem description$", re.I),
    re.compile(r"^typographical conventions", re.I),
    re.compile(r"^using this manual", re.I),
    re.compile(r"^how to use", re.I),
    re.compile(r"^where to find", re.I),
    re.compile(r"^what'?s\s+in\s+this", re.I),
    re.compile(r"^the contents of", re.I),
    re.compile(r"^for the beginner$", re.I),
    re.compile(r"^for the experienced user$", re.I),
    re.compile(r"^acknowledge?ments?$", re.I),
    re.compile(r"^references?$", re.I),
    re.compile(r"^bibliography$", re.I),
    re.compile(r"^index$", re.I),
    re.compile(r"^appendix\b", re.I),
    re.compile(r"^chapter\s+\d+:?\s*$", re.I),     # bare "Chapter 7:" with no title
    re.compile(r"^\d+(\.\d+)*\.?\s*$", re.I),       # pure section numbering "3.2.1."
    re.compile(r"^step\s+\d+\b", re.I),
    # tutorial-step boilerplate (after numbering strip)
    re.compile(r"^preparation$", re.I),
    re.compile(r"^obtaining the solution$", re.I),
    re.compile(r"^calculating a new solution", re.I),
    re.compile(r"^launching ansys fluent$", re.I),
    re.compile(r"^reading the (mesh|case) file$", re.I),
    re.compile(r"^setting up( the)?$", re.I),
    re.compile(r"^saving (the )?(case|mesh|data)", re.I),
    re.compile(r"^displaying( the)? preliminary solution$", re.I),
    re.compile(r"^summary of the (problem|solution)$", re.I),
    re.compile(r"^network of computers$", re.I),
    re.compile(r"^further improvements?$", re.I),
    re.compile(r"^settings?$", re.I),    # too generic to score
    re.compile(r"^solid$", re.I),
    re.compile(r"^liquid$", re.I),
    re.compile(r"^gas$", re.I),
    re.compile(r"^results?$", re.I),
    re.compile(r"^discussion$", re.I),
    re.compile(r"^postprocessing$", re.I),
    re.compile(r"^materials?$", re.I),
    re.compile(r"^mesh$", re.I),
    re.compile(r"^model$", re.I),
    re.compile(r"^geometry$", re.I),
    re.compile(r"^physics$", re.I),
    re.compile(r"^study$", re.I),
    re.compile(r"^definitions?$", re.I),
    re.compile(r"^domain$", re.I),
    re.compile(r"^global definitions?$", re.I),
    re.compile(r"^representations?,\s*warranties", re.I),
]


def strip_section_num(s: str) -> str:
    return _SECTION_NUM_RE.sub("", s).strip()


def is_na(feature: str, chapter: str) -> bool:
    f = strip_section_num(feature.strip())
    for rx in _NA_FEATURE_PATTERNS:
        if rx.search(f):
            return True
    # very short headings with no real keywords also n/a
    toks = tokenise(f"{chapter} {feature}")
    if not toks:
        return True
    return False


# ---------------------------------------------------------------------------
# Src indexing
# ---------------------------------------------------------------------------
def iter_src_files(root: Path) -> Iterable[Path]:
    for p in root.rglob("*"):
        if p.is_file() and (p.suffix.lower() in SRC_SUFFIXES or
                            p.name in SRC_NAMES):
            yield p


def subsystem_of(path: Path) -> str:
    rel = path.relative_to(SRC_DIR)
    return rel.parts[0] if rel.parts else "?"


def build_index() -> tuple[dict[str, list[tuple[str, str, int]]], int]:
    """token -> list of (subsystem, relpath, lineno).  Capped per token.

    Returns (index, nfiles_scanned).
    """
    PER_TOKEN_CAP = 50    # avoid blowing memory on common tokens
    idx: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    nfiles = 0
    for p in iter_src_files(SRC_DIR):
        nfiles += 1
        subsys = subsystem_of(p)
        rel = p.relative_to(ROOT).as_posix()
        try:
            for lineno, raw in enumerate(p.open("r", encoding="utf-8",
                                                errors="replace"), start=1):
                seen_in_line: set[str] = set()
                for tok in TOKEN_RE.findall(raw.lower()):
                    if tok in STOPWORDS or tok in seen_in_line:
                        continue
                    seen_in_line.add(tok)
                    bucket = idx[tok]
                    if len(bucket) < PER_TOKEN_CAP:
                        bucket.append((subsys, rel, lineno))
        except OSError:
            continue
    return idx, nfiles


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------
# Tokens with more index entries than this are dropped as diagnostic
# signal -- they match everywhere in src/ and prove nothing.  Per-token
# entries are capped at PER_TOKEN_CAP=50 during indexing so the threshold
# below acts on the *capped* count; using 40 means "appears in 40+ files"
# which is the noise floor we observed ("point", "data", "vector", ...).
DIAGNOSTIC_HIT_CAP = 40


def classify(feature: str,
             chapter: str,
             idx: dict[str, list[tuple[str, str, int]]]
             ) -> tuple[str, str]:
    if is_na(feature, chapter):
        return "n/a", ""

    tokens = tokenise(f"{chapter} {feature}")
    if not tokens:
        return "n/a", ""

    diagnostic: list[tuple[str, list[tuple[str, str, int]]]] = []
    for t in tokens:
        hits = idx.get(t)
        if not hits:
            continue
        if len(hits) > DIAGNOSTIC_HIT_CAP:
            continue   # too common to mean anything specific
        diagnostic.append((t, hits))

    if not diagnostic:
        return "absent", ""

    n = len(diagnostic)
    if n >= 3:
        status = "present"
    else:
        status = "partial"

    # Citations: smallest hit list first (most specific token).
    diagnostic.sort(key=lambda kv: len(kv[1]))
    cites: list[str] = []
    for tok, hits in diagnostic:
        _, rel, ln = hits[0]
        cites.append(f"{rel}:{ln}")
        if len(cites) >= 3:
            break
    return status, "; ".join(cites)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print(f"Phase 21.2 -- source-audit cross-walk")
    print(f"  Indexing {SRC_DIR.relative_to(ROOT)}/ ...")
    idx, nfiles = build_index()
    print(f"  Indexed {nfiles} files, {len(idx)} unique tokens")

    rows_in = TSV_IN.read_text(encoding="utf-8").splitlines()
    if not rows_in:
        print(f"ERROR: {TSV_IN} is empty")
        return 1
    header = rows_in[0].split("\t")
    # Old header: source chapter feature area our_status
    # New header: source chapter feature area our_status citations
    if header[:5] != ["source", "chapter", "feature", "area", "our_status"]:
        print(f"ERROR: unexpected header: {header}")
        return 1

    counts = {"present": 0, "partial": 0, "absent": 0, "n/a": 0}
    out_lines = ["\t".join(["source", "chapter", "feature", "area",
                            "our_status", "citations"])]
    for raw in rows_in[1:]:
        cols = raw.split("\t")
        # Pad to at least 5 columns
        while len(cols) < 5:
            cols.append("")
        source, chapter, feature, area, _old_status = cols[:5]
        status, cites = classify(feature, chapter, idx)
        counts[status] += 1
        out_lines.append("\t".join([source, chapter, feature, area,
                                    status, cites]))

    TSV_OUT.write_text("\n".join(out_lines) + "\n", encoding="utf-8",
                       newline="\n")

    total = sum(counts.values())
    print(f"\n  Total rows classified: {total}")
    for k in ("present", "partial", "absent", "n/a"):
        pct = 100.0 * counts[k] / total if total else 0.0
        print(f"    {k:8s}  {counts[k]:5d}   ({pct:5.1f}%)")
    print(f"\n  Wrote -> {TSV_OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
