#!/usr/bin/env python3
"""Make a vector research figure from completed publication validation records.

Optional plotting dependency: matplotlib==3.10.6. This script never runs solvers.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", type=Path, required=True)
    args = parser.parse_args()
    base = args.validation.resolve()
    summary = json.loads((base / "summary.json").read_text())
    if summary["failures"]:
        raise RuntimeError("Inspect failed validation cases before producing this figure")
    records = [json.loads(x) for x in (base / "designs.jsonl").read_text().splitlines()]
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="sparsedesign-matplotlib-"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
                         "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.75), layout="constrained")
    colors = {"S": "#17658c", "YI": "#b55931"}
    for family in ("S", "YI"):
        for lam in (0, 4):
            values = next(x for x in summary["synthetic_family_statistics"] if x["family"] == family and x["lambda_user"] == lam)
            axes[0].plot([3*n for n in values["amino_acid_lengths"]], values["candidates_over_n_squared"],
                         marker="o" if lam == 0 else "s", linestyle="-" if lam == 0 else "--",
                         color=colors[family], markersize=4, label=rf"{family}, $\lambda={lam}$")
    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].set_xticks([96, 192, 384, 768], labels=["96", "192", "384", "768"])
    axes[0].set_xlabel("RNA length n (nt)")
    axes[0].set_ylabel(r"Retained candidates $Z/n^2$")
    axes[0].set_title("A  Finite synthetic stress families", loc="left")
    axes[0].legend(frameon=False, ncols=2, loc="center left", bbox_to_anchor=(0.01, 0.5))
    axes[0].grid(axis="y", alpha=.18)
    proteins = sorted({r["case"]["protein"] for r in records if r["case"]["group"] == "constraints"}, key=len)
    for protein, color, marker in zip(proteins, ["#17658c", "#76608a"], ["o", "s"]):
        values = sorted((r for r in records if r["case"]["group"] == "constraints" and r["case"]["protein"] == protein and r["case"]["lambda_user"] == 0), key=lambda r:len(r["case"]["motifs"]))
        last = values[-1]
        axes[1].plot(range(len(values)), [r["nodes"]/len(r["rna"]) for r in values], marker=marker,
                     color=color, markersize=4, label=f"{len(protein)} aa (final width {last['width']})")
    axes[1].set_xticks(range(5), labels=["0", "1", "4", "8", "16"])
    axes[1].set_xlabel("Number of forbidden motifs (nested sets)")
    axes[1].set_ylabel("Automaton nodes per nucleotide N/n")
    axes[1].set_title("B  Constraint state growth", loc="left")
    axes[1].legend(frameon=False, loc="upper left")
    axes[1].grid(axis="y", alpha=.18)
    fig.savefig(base / "validation-robustness.pdf")
    fig.savefig(base / "validation-robustness.png", dpi=200)
    plt.close(fig)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    (base / "plot-provenance.json").write_text(json.dumps({"matplotlib": matplotlib.__version__, "script_sha256": digest(__file__), "summary_sha256": digest(base / "summary.json"), "designs_sha256": digest(base / "designs.jsonl")}, indent=2) + "\n")


if __name__ == "__main__":
    main()
