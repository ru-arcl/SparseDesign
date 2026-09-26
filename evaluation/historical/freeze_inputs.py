#!/usr/bin/env python3
"""Package the explicitly selected original records, never solver code/binaries."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-repo", required=True, type=Path)
    args = ap.parse_args()
    root = args.original_repo.resolve()
    out = Path(__file__).resolve().parent
    result = root / "research/results"
    selections = {}
    for group, dirname in [
        ("july13", "dystrophin-2026-07-13"),
        ("july23", "dystrophin-node1-j16-2026-07-23"),
        ("packed", "packed-layout-2026-07-14"),
    ]:
        for p in sorted((result / dirname).iterdir()):
            if p.is_file() and (p.suffix in {".out", ".time", ".err", ".meta", ".md", ".json", ".sh", ".log"}
                               or p.name in {"SHA256SUMS", "provenance.txt", "status.txt", "historical-dystseq.fasta", "codon-usage-clean.csv", "codon-usage-dense.csv"}
                               or p.name.startswith("frozen-")):
                selections[f"{group}/{p.name}"] = p
    for name in ["sparse-ablation-P30281-Q6UXY8-3x-2026-07-14.csv", "sparse-ablation-final-square-build.md", "sparse-ablation-build.md", "sparse-ablation-build-SHA256SUMS", "ablation-gxx-version.txt", "ablation-gxx-verbose.txt", "ablation-harness-ldd.txt", "ablation-python-version.txt", "hardware-arrakis-lscpu.json", "hardware-arrakis-cpu-node.txt", "hardware-arrakis-uname.txt", "hardware-arrakis-memory.txt"]:
        selections[f"ablation/{name}"] = result / name
    for name in ["packed-layout-2026-07-14-run.md", "packed-layout-2026-07-14-validation.json", "packed-layout-2026-07-14-validation.md"]:
        selections[f"packed/{name}"] = result / name
    for name in ["NP_000100.3_dystrophin_Dp427c.fasta", "NP_000100.3_dystrophin_Dp427c.md"]:
        selections[f"inputs/{name}"] = root / "research/data" / name
    files = []
    for logical, p in sorted(selections.items()):
        raw = p.read_bytes()
        stored = f"raw/{logical}.gz"
        compressed = gzip.compress(raw, mtime=0)
        dest = out / stored
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(compressed)
        files.append({"logical_path": logical, "original_path": str(p.relative_to(root)), "stored_path": stored,
                      "original_sha256": sha(raw), "original_bytes": len(raw),
                      "stored_sha256": sha(compressed), "stored_bytes": len(compressed)})
    source = []
    for name in ["analyze_sparse_ablation.py", "analyze_dystrophin_bench.py", "analyze_packed_layout.py"]:
        p = root / "research" / name
        raw = p.read_bytes()
        dest = out / "source" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
        source.append({"stored_path": f"source/{name}", "original_path": f"research/{name}", "sha256": sha(raw)})
    manifest = {"schema": "historical-input-bundle-v1", "files": files, "analysis_sources": source,
                "scope": "Original numeric, output, timing and provenance records only; no solver binaries or solver sources. Historical launch scripts are compressed documentation, never executed by the analyzer."}
    (out / "inputs.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Packaged {len(files)} original files and {len(source)} analysis sources")


if __name__ == "__main__":
    main()
