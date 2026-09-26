#!/usr/bin/env python3
"""Deterministic end-to-end checks for inverse lambda inference."""
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))

from frontier import infer_lambda  # noqa: E402


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else REPO / "lineardesign-clean").resolve()
    table = REPO / "data" / "codon_usage_freq_table_human.csv"
    options = {"bin_path": str(binary), "table": str(table), "threads": "1"}
    protein = "MDEYWHKR"

    mfe_design = infer_lambda(protein, "AUGGAUGAGUACUGGCACAAACGU", **options)
    check(mfe_design["supported"], "lambda=0 design should be supported")
    check(abs(mfe_design["lambda_min"]) < 1e-12, "MFE design should start at lambda=0")
    check(0.4 < mfe_design["lambda_max"] < 0.6, "unexpected finite support boundary")

    cai_design = infer_lambda(protein, "AUGGACGAGUACUGGCACAAGCGG", **options)
    check(cai_design["supported"], "CAI=1 design should be supported")
    check(0.4 < cai_design["lambda_min"] < 0.6, "unexpected CAI design boundary")
    check(cai_design["lambda_max"] is None and cai_design["unbounded"], "CAI=1 design should remain optimal indefinitely")

    dominated = infer_lambda(protein, "AUGGAUGAAUAUUGGCAUAAACGU", **options)
    check(not dominated["supported"], "dominated mRNA should have no supporting lambda")
    print("inverse lambda tests: 3 passed")


if __name__ == "__main__":
    main()
