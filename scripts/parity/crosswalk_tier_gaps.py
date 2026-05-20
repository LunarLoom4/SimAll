"""crosswalk_tier_gaps.py — Phase 21 Pass 21.3.

Reads ``docs/parity/feature_inventory.tsv`` (produced by Pass 21.2),
selects every row whose ``our_status`` is *absent* or *partial*, and
assigns each one a gap tier:

* **Tier-1 — table-stakes industrial CFD.**  Things a general-purpose
  solver simply must have: geometry import, meshing primitives, the
  standard set of boundary conditions, RANS turbulence (k-eps, k-omega,
  SST, Spalart-Allmaras), pressure-based segregated/coupled solvers,
  basic postprocessing (contour/vector/streamline/surface integrals),
  materials, parallel decomposition, scripting/journal, and the
  workbench-style UI flow.
* **Tier-2 — vertical-market depth.**  Features that distinguish a
  serious industrial code from a textbook one: density-based /
  compressible solvers, LES/DES/SAS scale-resolving turbulence,
  combustion (flamelet, PDF, EDC, NOx, soot), multiphase (VOF,
  Eulerian, mixture, cavitation, wall film), DPM/Lagrangian particles,
  radiation models (P1, DO, S2S, Monte Carlo), dynamic / sliding /
  overset mesh, MRF rotating frames, porous media, conjugate heat
  transfer, AMR.
* **Tier-3 — multi-physics extensions.**  Adjacent physics that ride on
  top of the CFD core: acoustics / aeroacoustics (FW-H), MHD /
  electromagnetics, fluid-structure interaction, shape optimisation /
  adjoint sensitivity, reduced-order models, immersed boundary,
  lattice-Boltzmann, electrochemistry / battery / fuel-cell, plasma.

The tier is decided by trying these in priority order:

1. Tier-3 keyword hit on chapter / feature text.
2. Tier-2 keyword hit on chapter / feature text.
3. Tier-1 keyword hit on chapter / feature text.
4. Fall-back by ``area`` (geometry/mesh/bc/... → Tier-1, turbulence/
   combustion/multiphase/dynamic → Tier-2, theory/other → Tier-3).

Output: ``docs/parity/gap_classification.tsv`` with schema::

    source  chapter  feature  area  our_status  tier  citations

Only *absent* and *partial* rows are written (the "present" rows have
no gap to classify, and "n/a" rows are doc-structural noise).
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "docs" / "parity" / "feature_inventory.tsv"
OUTPUT    = ROOT / "docs" / "parity" / "gap_classification.tsv"


# ---------------------------------------------------------------------------
# Keyword rules.  Compiled once.  First-match-wins inside each tier; tiers
# are tried in order T3 → T2 → T1 so the more specific physics overrides
# the generic CFD baseline.
# ---------------------------------------------------------------------------
def _compile(patterns: Iterable[str]) -> re.Pattern[str]:
    return re.compile("|".join(patterns), re.IGNORECASE)


_T3_KEYWORDS = _compile([
    r"\bacoustic", r"\baeroacoust", r"\bfw[-\s]?h\b", r"ffowcs",
    r"\bmhd\b", r"magnetohydro", r"electromagnet", r"electrochem",
    r"\bbattery\b", r"fuel[-\s]?cell", r"\bplasma\b",
    r"adjoint", r"shape\s+optim", r"topology\s+optim",
    r"reduced[-\s]?order", r"\brom\b",
    r"immersed\s+boundary", r"\bibm\b",
    r"fluid[-\s]?structure", r"\bfsi\b",
    r"lattice[-\s]?boltzmann", r"\blbm\b",
    r"\boptimization\b", r"\boptimisation\b",
    r"sensitivity\s+analysis", r"design\s+of\s+experiment", r"\bdoe\b",
])


_T2_KEYWORDS = _compile([
    # Scale-resolving turbulence.
    r"large[-\s]?eddy", r"\bles\b", r"detached[-\s]?eddy", r"\bdes\b",
    r"\bsas\b", r"wmles", r"scale[-\s]?adaptive",
    r"reynolds[-\s]?stress", r"\brsm\b",
    r"transition\s+model", r"gamma[-\s]?re[-\s]?theta",
    # Radiation.
    r"radiation", r"\bp1\b", r"\bp-1\b", r"discrete[-\s]?ordinate",
    r"\bs2s\b", r"surface[-\s]?to[-\s]?surface", r"monte[-\s]?carlo",
    r"view\s+factor",
    # Multiphase / particles / film.
    r"\bdpm\b", r"discrete[-\s]?phase", r"lagrangian\s+particle",
    r"\bspray\b", r"\binjection", r"breakup",
    r"eulerian", r"\bvof\b", r"volume[-\s]?of[-\s]?fluid",
    r"level[-\s]?set", r"mixture\s+model", r"cavitat", r"wall[-\s]?film",
    r"phase[-\s]?change", r"\bmelt", r"\bsolidif", r"evaporat", r"condens",
    # Dynamic mesh / rotating / overset.
    r"moving\s+mesh", r"sliding\s+mesh", r"\bmrf\b",
    r"multiple\s+reference\s+frame", r"dynamic\s+mesh",
    r"\b6\s*dof\b", r"six[-\s]?dof", r"overset", r"chimera",
    r"turbomachinery", r"\brotor\b", r"\bstator\b",
    r"fan\s+model", r"blade\s+passage",
    # Reactions / combustion (beyond area=combustion).
    r"chemkin", r"reaction\s+mechanism", r"kinetic\s+scheme",
    r"\bpdf\b", r"flamelet", r"eddy[-\s]?dissipation", r"\bedc\b",
    r"\bspark\b", r"\bsoot\b", r"\bnox\b",
    # Porous / heat exchanger.
    r"porous", r"darcy", r"forchheimer",
    r"heat\s+exchanger", r"conjugate\s+heat", r"\bcht\b",
    # AMR / density-based.
    r"adaptive\s+mesh", r"hanging[-\s]?node",
    r"density[-\s]?based", r"coupled\s+solver",
    r"compressible", r"supersonic", r"hypersonic", r"\bshock\b",
])


_T1_KEYWORDS = _compile([
    # Geometry / mesh basics.
    r"\bhexa\b", r"\btetra\b", r"\bprism\b", r"\bpolyhedra",
    r"inflation\s+layer", r"prism\s+layer", r"boundary\s+layer\s+mesh",
    r"\bquad\b", r"\btri\b", r"surface\s+mesh", r"volume\s+mesh",
    r"meshing", r"refinement", r"sizing",
    # Boundary / initial conditions.
    r"\binlet\b", r"\boutlet\b", r"\bwall\b", r"symmetry", r"periodic",
    r"pressure[-\s]?outlet", r"velocity[-\s]?inlet", r"mass[-\s]?flow",
    r"boundary\s+condition", r"initial\s+condition",
    # Baseline RANS.
    r"k-?\s*epsilon", r"k-?\s*omega", r"\bsst\b",
    r"spalart", r"\bsa\s+model\b",
    # Pressure-based solver loop.
    r"\bsimple\b", r"\bsimplec\b", r"\bpiso\b", r"\bcoupled\b",
    r"pressure[-\s]?based",
    r"residual", r"convergence", r"\biteration", r"courant", r"\bcfl\b",
    r"time[-\s]?step",
    # Postprocessing basics.
    r"\breport\b", r"\bmonitor\b", r"\bplot\b", r"contour\s+plot",
    r"vector\s+plot", r"streamline", r"surface\s+integral",
    r"volume\s+integral", r"\bxy\s+plot\b",
    # Parallel.
    r"\bparallel\b", r"\bpartition", r"\bmetis\b", r"\bmpi\b",
    r"domain\s+decomposition",
    # Scripting / IO.
    r"\bjournal\b", r"\btui\b", r"\bbatch\b",
    r"\bpython\s+script", r"\bmacro\b",
    r"case\s+file", r"data\s+file", r"\bcgns\b", r"\bensight\b",
    r"read\s+mesh", r"write\s+mesh", r"export", r"import",
    # Materials.
    r"material\s+propert", r"newton", r"viscosity", r"density",
    r"specific\s+heat", r"thermal\s+conductivity",
])


# Area fall-back when no keyword fires.
_AREA_TIER: dict[str, str] = {
    "geometry":       "T1",
    "mesh":           "T1",
    "bc":             "T1",
    "materials":      "T1",
    "solver":         "T1",
    "parallel":       "T1",
    "scripting":      "T1",
    "io":             "T1",
    "ui":             "T1",
    "workflow":       "T1",
    "postprocessing": "T1",
    "validation":     "T1",
    "turbulence":     "T2",
    "combustion":     "T2",
    "multiphase":     "T2",
    "dynamic":        "T2",
    "theory":         "T3",
    "other":          "T3",
}


def classify(chapter: str, feature: str, area: str) -> str:
    """Return ``"T1"`` / ``"T2"`` / ``"T3"`` for a single inventory row."""
    text = f"{chapter}  {feature}"
    if _T3_KEYWORDS.search(text):
        return "T3"
    if _T2_KEYWORDS.search(text):
        return "T2"
    if _T1_KEYWORDS.search(text):
        return "T1"
    return _AREA_TIER.get(area, "T3")


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def main() -> int:
    if not INVENTORY.exists():
        sys.stderr.write(
            f"[crosswalk_tier_gaps] inventory not found: {INVENTORY}\n"
            "  run scripts/parity/extract_feature_inventory.py first\n"
        )
        return 2

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    rows_in  = 0
    rows_out = 0
    by_status_tier: Counter[tuple[str, str]] = Counter()
    by_source_tier: Counter[tuple[str, str]] = Counter()

    with INVENTORY.open("r", encoding="utf-8", newline="") as fh_in, \
         OUTPUT.open("w", encoding="utf-8", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")
        writer = csv.writer(fh_out, delimiter="\t",
                            lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["source", "chapter", "feature", "area",
                         "our_status", "tier", "citations"])

        for row in reader:
            rows_in += 1
            status = row["our_status"]
            if status not in ("absent", "partial"):
                continue
            tier = classify(row["chapter"], row["feature"], row["area"])
            writer.writerow([
                row["source"], row["chapter"], row["feature"],
                row["area"], status, tier, row.get("citations", ""),
            ])
            rows_out += 1
            by_status_tier[(status, tier)] += 1
            by_source_tier[(row["source"], tier)] += 1

    # -- diagnostics --------------------------------------------------------
    sys.stdout.write(
        f"[crosswalk_tier_gaps] wrote {OUTPUT.relative_to(ROOT).as_posix()} "
        f"({rows_out} rows from {rows_in} scanned)\n"
    )
    sys.stdout.write("\n  by status x tier:\n")
    for (status, tier), n in sorted(by_status_tier.items()):
        sys.stdout.write(f"    {status:<8} {tier}: {n:>5}\n")
    sys.stdout.write("\n  by source x tier:\n")
    for (source, tier), n in sorted(by_source_tier.items()):
        sys.stdout.write(f"    {source:<10} {tier}: {n:>5}\n")

    # Sanity gate: at least the documented 155 absent rows must survive.
    absent_total = sum(n for (s, _), n in by_status_tier.items() if s == "absent")
    if absent_total < 100:
        sys.stderr.write(
            f"[crosswalk_tier_gaps] WARN: only {absent_total} absent rows "
            "tiered (expected ~155 from Pass 21.2)\n"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
