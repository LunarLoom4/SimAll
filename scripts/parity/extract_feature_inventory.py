#!/usr/bin/env python3
"""Phase 21.1 -- Table-of-contents extractor for the competitor parity matrix.

Streams the six advertised source documents at the workspace root, walks
every Markdown heading, and emits one row per heading into
docs/parity/feature_inventory.tsv with columns

    source  chapter  feature  area  our_status

`source` is the short-code (Fluent-TG, OF-UG, ...).
`chapter` is the nearest H1/H2 ancestor heading text.
`feature` is the heading line itself (any level).
`area` is a coarse bucket inferred by keyword scan
       (geometry, mesh, bc, turbulence, combustion, multiphase,
        dynamic, postprocessing, parallel, materials, solver,
        scripting, io, ui, workflow, validation, theory, other).
`our_status` is left blank on purpose -- Pass 21.2 (source-audit
       cross-walk) fills it in with file/line citations and one of
       {present, partial, absent, n/a}.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "docs" / "parity"
OUT_TSV = OUT_DIR / "feature_inventory.tsv"

SOURCES = {
    "ANSYS_Fluent_18_Tutorial_Guide.md":    "Fluent-TG",
    "COMSOL_MultiPhysics_User_Guide.md":    "CMP-UG",
    "COMSOL_CFD_Module_User_Guide.md":      "CFD-UG",
    "COMSOL_CFD_Module_Software_Guide.md":  "CFD-SG",
    "OpenFOAM_User_Guide.md":               "OF-UG",
    "OpenFOAM_Programmers_Guide.md":        "OF-PG",
}

# Ordered so earlier rules win on ambiguous text.
AREA_KEYWORDS = [
    ("geometry",       ["geometry", "sketch", "feature tree", "cad", "step file",
                        "iges", "brep", "solid model", "surface model", "primitive",
                        "extrud", "revolv", "loft", "sweep", "fillet", "chamfer",
                        "boolean", "imprint", "share topology", "defeatur", "heal"]),
    ("mesh",           ["mesh", "grid", "tessell", "prism layer", "hex-core",
                        "hexcore", "polyhedral", "tetrahedral", "tetrahed",
                        "octree", "wrap", "snap", "remesh", "size field",
                        "smoothing", "edge swap"]),
    ("bc",             ["boundary condition", "boundary cond", "inlet", "outlet",
                        "wall function", "wall treatment", "symmetry", "periodic",
                        "non-conformal", "far-field", "porous jump", "radiator",
                        "interior", "axis"]),
    ("turbulence",     ["turbul", "rans", " les", "des ", " sst", "k-epsilon",
                        "k-omega", "spalart", "reynolds stress", "y+", "y-plus",
                        "hybrid model"]),
    ("combustion",     ["combust", "reaction", "ignit", "spark", "flame",
                        "species transport", "chemistry", "flamelet", " pdf ",
                        "mixture fraction", "equilibrium", "pollutant", "soot",
                        "nox", "fgm", "slfm"]),
    ("multiphase",     ["multiphase", " vof", "volume of fluid", "eulerian",
                        "mixture model", " dpm", "lagrangian", "cavit",
                        "wall film", "population balance", "drift"]),
    ("dynamic",        ["dynamic mesh", "moving mesh", "sliding mesh",
                        "mesh motion", "rotating", "mrf", "6dof", "sdof",
                        "rigid body", "overset", "morphing", "rbf"]),
    ("postprocessing", ["post-process", "postprocess", "contour", "iso-surface",
                        "isosurface", "streamline", "pathline", "vortex core",
                        "monitor", "report definit", "animation", "section plane",
                        "probe", "scalar bar"]),
    ("parallel",       ["parallel", "mpi", "cluster", "domain decomposit",
                        "gpu", "cuda", "openmp", "scaling", "partition"]),
    ("materials",      ["material", "fluid propert", "solid propert",
                        "equation of state", "viscos", "thermal conductivity",
                        "specific heat", "density model"]),
    ("solver",         ["solver", "scheme", "pressure-velocity", "pressure based",
                        "density based", "coupled", "simple", "piso", "pimple",
                        "amg ", "multigrid", "preconditioner", "convergence",
                        "residual", "timestep", "time step", "cfl",
                        "numerics", "discretiz", "discretis", "linear solver"]),
    ("scripting",      ["scripting", "journal", " udf", "python", " tui ",
                        "macro", "expression", "parameter", "udm"]),
    ("io",             ["import file", "export file", "case file", "write data",
                        "read data", "file format", " vtk", " cgns", "ensight",
                        "tecplot", "fluent format"]),
    ("ui",             ["panel", "dialog", "menu", "ribbon", "tree view",
                        "schematic", "workbench", "graphical user", "interface",
                        "window", "navigator", "outline view"]),
    ("workflow",       ["workflow", "task page", "wizard", "setup task",
                        "guided"]),
    ("validation",     ["validation", "verification", "benchmark", "tutorial",
                        "example case", "test case"]),
    ("theory",         ["theory", "governing equation", "conservation",
                        "derivation"]),
]

HEADING_RE = re.compile(r"^(#{1,6})\s+(.+?)\s*$")


def classify_area(text: str) -> str:
    s = text.lower()
    for area, kws in AREA_KEYWORDS:
        for kw in kws:
            if kw in s:
                return area
    return "other"


def clean_cell(s: str) -> str:
    return (
        s.replace("\t", " ")
         .replace("\r", " ")
         .replace("\n", " ")
         .strip()
    )


def extract(md_path: Path, source_code: str):
    rows = []
    h1 = ""
    h2 = ""
    with md_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            m = HEADING_RE.match(raw.rstrip("\n"))
            if not m:
                continue
            level = len(m.group(1))
            title = m.group(2).strip().strip("*").strip()
            if not title:
                continue
            if level == 1:
                h1 = title
                h2 = ""
            elif level == 2:
                h2 = title

            if level <= 2:
                chapter = h2 if h2 else h1
            else:
                chapter = h2 or h1 or "(root)"

            area = classify_area(f"{chapter} {title}")
            rows.append((source_code,
                         clean_cell(chapter),
                         clean_cell(title),
                         area,
                         ""))
    return rows


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_rows = []
    totals = []
    for fname, code in SOURCES.items():
        p = ROOT / fname
        if not p.exists():
            print(f"WARN: missing source file {fname}")
            continue
        rows = extract(p, code)
        totals.append((code, fname, len(rows)))
        all_rows.extend(rows)

    with OUT_TSV.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("source\tchapter\tfeature\tarea\tour_status\n")
        for r in all_rows:
            fh.write("\t".join(r) + "\n")

    print("Phase 21.1 -- feature inventory extraction")
    print("-" * 60)
    for code, fname, n in totals:
        print(f"  {code:10s}  {n:5d} rows   {fname}")
    print("-" * 60)
    print(f"  TOTAL       {len(all_rows):5d} rows -> "
          f"{OUT_TSV.relative_to(ROOT).as_posix()}")
    if len(all_rows) < 600:
        print(f"  WARNING: target was >=600 rows, got {len(all_rows)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
