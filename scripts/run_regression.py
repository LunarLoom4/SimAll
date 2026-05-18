#!/usr/bin/env python3
"""SimAll Beta — regression driver.

Runs the ctest label group `regression`, parses CTest XML, and emits
a JUnit-style summary to build/regression-results.xml so CI can
report fine-grained failures.

Usage:
    python scripts/run_regression.py [--build-dir build] [--preset default]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def run(cmd: list[str], cwd: Path | None = None) -> int:
    print("$", " ".join(cmd), flush=True)
    return subprocess.call(cmd, cwd=cwd)


def collect_ctest_xml(build_dir: Path) -> Path | None:
    tag = build_dir / "Testing" / "TAG"
    if not tag.exists():
        return None
    timestamp = tag.read_text().splitlines()[0].strip()
    out = build_dir / "Testing" / timestamp / "Test.xml"
    return out if out.exists() else None


def to_junit(test_xml: Path, junit_xml: Path) -> tuple[int, int]:
    root = ET.parse(test_xml).getroot()
    cases = []
    failed = 0
    for t in root.findall(".//Test"):
        status = t.get("Status", "passed")
        name   = t.findtext("Name") or "?"
        time   = t.find(".//NamedMeasurement[@name='Execution Time']/Value")
        secs   = float(time.text) if time is not None else 0.0
        case   = ET.Element("testcase", {"name": name, "classname": "simall.regression", "time": f"{secs:.3f}"})
        if status != "passed":
            failed += 1
            msg = t.findtext(".//Measurement/Value") or status
            ET.SubElement(case, "failure", {"message": status}).text = msg
        cases.append(case)
    suite = ET.Element("testsuite", {
        "name":  "simall-regression",
        "tests": str(len(cases)),
        "failures": str(failed),
    })
    suite.extend(cases)
    ET.ElementTree(suite).write(junit_xml, encoding="utf-8", xml_declaration=True)
    return len(cases), failed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--preset",    default=None,
                        help="optional ctest preset name")
    parser.add_argument("--label",     default="regression")
    args = parser.parse_args()

    cmd = ["ctest", "--output-on-failure", "-L", args.label]
    if args.preset:
        cmd = ["ctest", "--preset", args.preset, "--output-on-failure", "-L", args.label]
    rc = run(cmd, cwd=args.build_dir)

    xml = collect_ctest_xml(args.build_dir)
    if xml is None:
        print("No CTest XML produced — was the build configured for testing?", file=sys.stderr)
        return rc or 2

    junit = args.build_dir / "regression-results.xml"
    n, fails = to_junit(xml, junit)
    print(f"Regression summary: {n} tests, {fails} failures -> {junit}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
