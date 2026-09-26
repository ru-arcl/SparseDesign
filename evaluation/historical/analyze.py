#!/usr/bin/env python3
"""Check frozen historical evidence and regenerate claims; Python 3.10+, stdlib only."""
from __future__ import annotations
import argparse
from datetime import datetime
from decimal import Decimal
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import statistics
import sys
import tempfile
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "source" / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def unpack(target):
    manifest = json.loads((ROOT / "inputs.json").read_text())
    require(manifest["schema"] == "historical-input-bundle-v1", "unknown input schema")
    seen = set()
    for entry in manifest["files"]:
        logical = Path(entry["logical_path"])
        stored = Path(entry["stored_path"])
        require(not logical.is_absolute() and ".." not in logical.parts and not stored.is_absolute() and ".." not in stored.parts, "unsafe bundle path")
        require(str(logical) not in seen, "duplicate bundle member")
        seen.add(str(logical))
        packed = (ROOT / stored).read_bytes()
        require(sha(packed) == entry["stored_sha256"] and len(packed) == entry["stored_bytes"], f"compressed input hash/size mismatch: {stored}")
        raw = gzip.decompress(packed)
        require(sha(raw) == entry["original_sha256"] and len(raw) == entry["original_bytes"], f"original input hash/size mismatch: {logical}")
        dest = target / logical
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
    for entry in manifest["analysis_sources"]:
        require(sha((ROOT / entry["stored_path"]).read_bytes()) == entry["sha256"], "analysis source hash mismatch")
    return manifest


def jul23(root, dyst):
    provenance = (root / "provenance.txt").read_text()
    require("node1 CPUs16-31 j16 sequential-lineardesign-then-clean" in provenance, "July23 protocol mismatch")
    events = re.findall(r"\[([^]]+)\] (START|END) (lineardesign|clean)([^\n]*)", (root / "progress.log").read_text())
    require([(e[1], e[2]) for e in events] == [("START", "lineardesign"), ("END", "lineardesign"), ("START", "clean"), ("END", "clean")], "July23 chronology incomplete")
    times = [datetime.fromisoformat(e[0]) for e in events]
    require(times == sorted(times), "July23 chronology nonsequential")
    require(all("status=0" in events[i][3] for i in [1, 3]), "July23 exit event failure")
    hashes = {Path(p).name: h for h, p in re.findall(r"^([a-f0-9]{64})  (.+)$", provenance, flags=re.M)}
    require(hashes["LinearDesign_2D"] == "50040b7e7956c8c69358a992762a82909d5b0361acd8a0d0ff4f5c920ba2bcbe", "July23 dense binary identity changed")
    require(hashes["lineardesign-clean"] == "3772e1e1109711b266c8eabf52eb0dab2c9e60689e1d0f9fd99ceeec9d71f885", "July23 sparse binary identity changed")
    rows = []
    for idx, name in enumerate(["lineardesign", "clean"]):
        stem = name + "_j16"
        require((root / (stem + ".err")).read_bytes() == b"", f"{stem}: nonempty stderr")
        raw = (root / (stem + ".out")).read_bytes()
        design = dyst.parse_design(raw, stem)
        timing = dyst.parse_timing((root / (stem + ".time")).read_bytes(), stem)
        args = dyst.command_arguments(timing.command, stem)
        binder_index = next(i for i, item in enumerate(args) if Path(item).name == "numabind_epyc_node")
        require(args[binder_index + 1] == "1", "July23 node mismatch")
        if name == "lineardesign":
            require(args[-4:] == ["0", "0", "./codon_usage_freq_table_human.csv", "16"], "July23 dense CLI mismatch")
        else:
            require(args[-6:-3] == ["-l", "0", "-c"] and args[-2:] == ["-j", "16"], "July23 sparse CLI mismatch")
        require(all(item in args for item in ["OMP_PROC_BIND=close", "OMP_PLACES=cores", "OMP_DYNAMIC=false", "OMP_WAIT_POLICY=active", "OMP_THREAD_LIMIT=16"]), "July23 OpenMP protocol mismatch")
        window = (times[2 * idx + 1] - times[2 * idx]).total_seconds()
        require(abs(window - float(timing.wall_seconds)) <= 2, "July23 timing/chronology mismatch")
        require(len(design.sequence) == 11031 and design.mfe == Decimal("-7161.40"), "July23 output length/MFE mismatch")
        rows.append({"run_id": stem, "wall_seconds": float(timing.wall_seconds), "max_rss_kib": timing.max_rss_kib,
                     "max_rss_gib": timing.max_rss_kib / 1048576, "mfe_kcal": str(design.mfe), "cai": str(design.cai),
                     "rna_sha256": sha(design.sequence.encode()), "structure_sha256": sha(design.structure.encode()),
                     "command": timing.command, "threads": 16, "node": 1, "binary_sha256": hashes["LinearDesign_2D" if name == "lineardesign" else "lineardesign-clean"]})
    return {"protocol": "July23 node1 sequential dense-then-sparse pair; one run per implementation", "runs": rows,
            "dense_over_sparse_wall": rows[0]["wall_seconds"] / rows[1]["wall_seconds"], "recorded_dependency_sha256": hashes}


def analyze():
    with tempfile.TemporaryDirectory(prefix="historical-evidence-") as tmp:
        root = Path(tmp)
        manifest = unpack(root)
        abl = module("analyze_sparse_ablation")
        dyst = module("analyze_dystrophin_bench")
        packed = module("analyze_packed_layout")
        args = SimpleNamespace(diagnostic=False, expected_accessions=2, expected_accession_names=["P30281", "Q6UXY8"],
                               expected_lambdas=["0", "4"], expected_threads=[1, 16], expected_repeats=3, bootstrap=10000, seed=20260714)
        snapshot = abl.read_snapshot(root / "ablation/sparse-ablation-P30281-Q6UXY8-3x-2026-07-14.csv", args)
        require(len(snapshot.valid_pairs) == 24, "ablation expected 24 matched pairs")
        ablation = abl.summary_rows(snapshot, args)
        obs, verified = dyst.collect_observations(root / "july13", dyst.GROUPS)
        require(verified and len(obs) == 27, "July13 incomplete original manifest/grid")
        july13 = [dyst.observation_row(o) for o in obs]
        july23 = jul23(root / "july23", dyst)
        pk = packed.analyze_archive(root / "packed")
        pk["archive"] = "raw/packed (original byte streams restored in temporary directory)"
        source = (root / "inputs/NP_000100.3_dystrophin_Dp427c.fasta").read_text()
        protein = "".join(line.strip() for line in source.splitlines() if not line.startswith(">"))
        require(len(protein) == 3677 and sha(protein.encode()) == "c0b35a4bba1c1cece0cb576cbf1302a2047b94bbe218844e692cfcde480ff2f6", "dystrophin input identity mismatch")
        historical = (root / "july13/historical-dystseq.fasta").read_text()
        require(protein == "".join(line.strip() for line in historical.splitlines() if not line.startswith(">")), "different July13/July23 input sequence")
        old_manifest = dyst.parse_manifest(root / "july13")
        historical_readme = (root / "july13/README.md").read_text()
        original_builds = {
            "ablation": {key: snapshot.valid_pairs[0].sparse.values[key] for key in ["harness_sha256", "codon_table_sha256", "collector_sha256", "command_prefix_sha256", "gnu_time_sha256"]},
            "ablation_source_record": "raw/ablation/sparse-ablation-final-square-build.md.gz",
            "july13": {
                "sparse_binary_sha256": old_manifest["lineardesign-clean-881bf4513e51"],
                "local_dense_binary_sha256": old_manifest["LinearDesign_2D-local-openmp-50040b7e7956"],
                "local_dense_commit_as_recorded": re.search(r"source commit `([a-f0-9]{40})`", historical_readme).group(1),
                "sparse_source_boundary": "Historical v0039 label and binary hash retained; no complete frozen sparse-source tree is supplied by this archive.",
            },
            "july23": {"recorded_hashes": july23["recorded_dependency_sha256"],
                       "source_boundary": "The sparse binary identity differs from July13; a binary hash does not establish complete source identity."},
            "packed": pk["frozen"],
            "verification_boundary": "Original identity records are authenticated as retained bytes. Omitted historical binaries and solver source trees are not independently rebuilt or rehashed by this package.",
        }
        result = {"schema": "historical-summary-v1", "status": "pass", "input_manifest_sha256": sha((ROOT / "inputs.json").read_bytes()),
                  "analysis_script_sha256": sha(Path(__file__).read_bytes()), "archived_files_verified": len(manifest["files"]),
                  "ablation": ablation, "july13": july13, "july23": july23, "packed": pk, "original_build_identities": original_builds,
                  "scope": "Reanalysis of original records, not new execution. MFE/CAI values are printed solver outputs; independent archived validation reports are retained but not rerun."}
        return result


def outputs(summary):
    rows = {r["run_id"]: r for r in summary["july13"]}
    q6 = [r for r in summary["ablation"] if r["accession"] == "Q6UXY8"]
    noisy = next(r for r in q6 if r["lambda"] == "0" and r["threads_requested"] == 16)
    repeats = [float(r["wall_seconds"]) for r in rows.values() if r["group"] == "rep"]
    parallel = float(rows["clean_j1"]["wall_seconds"]) / float(rows["clean_j16"]["wall_seconds"])
    pk = summary["packed"]
    square = pk["layouts"]["B-square"]
    packed = pk["layouts"]["D-packed"]
    vals = {
        "HistAblationPairs": (sum(r["pairs"] for r in summary["ablation"]), 0),
        "HistAblationSolveLow": (min(r["solve_speedup_median"] for r in q6), 3),
        "HistAblationSolveHigh": (max(r["solve_speedup_median"] for r in q6), 3),
        "HistAblationNoisyLow": (noisy["solve_speedup_min"], 3), "HistAblationNoisyHigh": (noisy["solve_speedup_max"], 3),
        "HistAblationRSSMiB": (statistics.median(r["rss_dense_minus_sparse_median_kib"] / 1024 for r in q6), 0),
        "HistDystrophinAA": (len(rows) and int(rows["clean_j1"]["protein_aa"]), 0), "HistDystrophinNT": (int(rows["clean_j1"]["rna_nt"]), 0),
        "HistDystrophinNodes": (int(pk["cost_metadata"]["lattice_nodes"]), 0), "HistDystrophinMFE": (float(rows["clean_j1"]["mfe_kcal"]), 2),
        "HistParallelGain": (parallel, 2), "HistParallelEfficiencyPercent": (100 * parallel / 16, 0),
        "HistJulyThirteenRatio": (float(rows["ld_j16"]["wall_seconds"]) / float(rows["clean_j16"]["wall_seconds"]), 2),
        "HistJulyTwentyThreeRatio": (summary["july23"]["dense_over_sparse_wall"], 2),
        "HistRepeatMin": (min(repeats), 2), "HistRepeatMax": (max(repeats), 2),
        "HistPackedSolveIncreasePercent": (100 * (packed["median_solve_seconds"] / square["median_solve_seconds"] - 1), 2),
        "HistPackedRSSReductionPercent": (100 * (1 - packed["median_max_rss_kib"] / square["median_max_rss_kib"]), 2),
        "HistPackedRSSGiB": (packed["median_max_rss_kib"] / 1048576, 3),
    }
    claims = "% Generated by evaluation/historical/analyze.py; do not edit.\n" + "".join(f"\\newcommand{{\\{key}}}{{{value:.{places}f}}}\n" for key, (value, places) in vals.items())
    table = ["% Generated historical table body; protocols remain separate."]
    for protocol, label, run_id in [("July 13, node 0 ladder", "Historical dense", "ld_j1"), ("", "Sparse square", "clean_j1"), ("", "Historical dense", "ld_j16"), ("", "Sparse square", "clean_j16")]:
        r = rows[run_id]
        table.append(f"{protocol} & {label} & {r['threads']} & {float(r['wall_seconds']):,.2f} & {int(r['max_rss_kib']) / 1048576:.2f}\\\\")
    for i, r in enumerate(summary["july23"]["runs"]):
        table.append(f"{'July 23, node 1 pair' if i == 0 else ''} & {'Historical dense' if i == 0 else 'Sparse square'} & 16 & {r['wall_seconds']:,.2f} & {r['max_rss_gib']:.2f}\\\\")
    claims_json = {key: value for key, (value, places) in vals.items()}
    return {"summary.json": json.dumps(summary, indent=2, sort_keys=True) + "\n", "claims.json": json.dumps(claims_json, indent=2, sort_keys=True) + "\n",
            "claims.tex": claims, "dystrophin-table.tex": "\n".join(table) + "\n",
            "historical-results.tex": claims + "\\newcommand{\\HistDystrophinTableRows}{%\n" + "\n".join(table) + "\n}\n"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="compare all generated files without changing them")
    ap.add_argument("--output", type=Path, help="write regenerated files in this directory")
    args = ap.parse_args()
    try:
        for name, content in outputs(analyze()).items():
            dest = (args.output or ROOT) / name
            if args.check:
                require(dest.read_text() == content, f"generated output drift: {dest}")
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_text(content)
        print("PASS: historical input/source hashes, 24 ablation pairs, 27 July13 runs, 2 July23 runs, 5 packed runs, and generated claims")
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
