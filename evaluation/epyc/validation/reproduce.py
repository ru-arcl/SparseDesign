#!/usr/bin/env python3
"""Rebuild EPYC validation summaries and manuscript evidence without running solvers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true", help="Regenerate the PDF/PNG with matplotlib")
    parser.add_argument("--paper", type=Path, help="Copy summary, claims and figure into this paper directory")
    args = parser.parse_args()
    base = Path(__file__).resolve().parent
    manifest = json.loads((base / "raw-manifest.json").read_text())
    for name, expected in manifest["sha256"].items():
        if digest(base / name) != expected:
            raise RuntimeError("Recorded input/source hash mismatch: " + name)
    build = json.loads((base / "build.json").read_text())
    for name, expected in build["source_sha256"].items():
        if digest(base / "source/solver" / name) != expected:
            raise RuntimeError("Executed solver source mismatch: " + name)
    for name, key in (("publication_validation.py", "script_sha256"),
                      ("publication_validation.cc", "driver_sha256")):
        if digest(base / "source" / name) != build[key]:
            raise RuntimeError("Executed validation source mismatch: " + name)

    analysis = base / "source/analysis/publication_validation.py"
    result = runpy.run_path(str(analysis))["summary"](base)
    if (result["failures"] or result["loops"]["passed"] != result["loops"]["cases"]
            or result["exhaustive"]["passed"] != result["exhaustive"]["comparisons"]
            or not result["constraint_feasible_set_monotonicity_pass"]
            or not result["codon_precision_common_objective_nonnegative_regret"]):
        raise RuntimeError("EPYC validation checks failed")

    plotter = base / "source/analysis/plot_publication_validation.py"
    if args.plot:
        env = dict(os.environ)
        env["SOURCE_DATE_EPOCH"] = "1790366400"  # 2026-09-25 20:00 UTC; stable PDF metadata.
        subprocess.run([sys.executable, str(plotter), "--validation", str(base)],
                       check=True, env=env)
    if args.paper:
        plot = json.loads((base / "plot-provenance.json").read_text())
        for key, path in (("summary_sha256", base / "summary.json"),
                          ("designs_sha256", base / "designs.jsonl"),
                          ("script_sha256", plotter)):
            if plot[key] != digest(path):
                raise RuntimeError("Figure provenance mismatch; regenerate with --plot: " + key)
        for source, destination in (("summary.json", "data/validation-summary.json"),
                                    ("claims.tex", "generated/validation-claims.tex"),
                                    ("validation-robustness.pdf", "figures/validation-robustness.pdf")):
            target = args.paper / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(base / source, target)
    print(f"EPYC validation: {result['passed_design_cases']}/{result['total_design_cases']} designs, "
          f"{result['loops']['passed']}/{result['loops']['cases']} loops, "
          f"{result['exhaustive']['passed']}/{result['exhaustive']['comparisons']} exhaustive checks")


if __name__ == "__main__":
    main()
