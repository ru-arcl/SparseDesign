#!/usr/bin/env python3
"""Read-only similarity clustering of the frozen 2,000-protein panel.

This module never invokes a sequence-design solver. It groups archived proteins
solely for sensitivity analysis of already recorded candidate-count statistics.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import zipfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def prepare(root, out):
    campaign = root / "research/results/sparse-balanced-codon-panel-2026-07-16"
    plan = json.loads((campaign / "plan-final.json").read_text())
    manifest = root / "research/data/uniprot-sprot-2026_02-sparse-balanced-2000.manifest.csv"
    archive = root / "research/data/uniprot-sprot-2026_02-sparse-balanced-2000.zip"
    if sha(manifest.read_bytes()) != plan["inputs"]["manifest_csv"]["sha256"] or sha(archive.read_bytes()) != plan["inputs"]["source"]["sha256"]:
        raise ValueError("Archived panel manifest/FASTA zip hash mismatch")
    rows = list(csv.DictReader(io.StringIO(manifest.read_text())))
    if len(rows) != 2000 or len({r["accession"] for r in rows}) != 2000:
        raise ValueError("Expected 2,000 unique protein accessions")
    fasta, index = [], []
    with zipfile.ZipFile(archive) as z:
        for row in sorted(rows, key=lambda r:r["accession"]):
            accession = row["accession"]
            if not re.fullmatch("[A-Z0-9]+", accession):
                raise ValueError("Unsafe FASTA identifier")
            raw = z.read(row["fasta_path"])
            if sha(raw) != row["fasta_sha256"]:
                raise ValueError("Individual frozen FASTA hash mismatch")
            lines = raw.decode().splitlines()
            if sum(line.startswith(">") for line in lines) != 1:
                raise ValueError("Expected a single archived FASTA record")
            sequence = "".join(line for line in lines if not line.startswith(">"))
            if sha(sequence.encode()) != row["sequence_sha256"] or len(sequence) != int(row["aa_length"]):
                raise ValueError("Archived protein sequence mismatch")
            fasta.append(f">{accession}\n{sequence}\n")
            index.append({k:row[k] for k in ("accession", "sequence_sha256", "aa_length", "length_bin")})
    out.mkdir(parents=True, exist_ok=True)
    raw_fasta = "".join(fasta).encode()
    (out / "panel.fasta.gz").write_bytes(gzip.compress(raw_fasta, mtime=0))
    with (out / "sequence-index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index[0]))
        writer.writeheader(); writer.writerows(index)
    provenance = {"source_zip_sha256": sha(archive.read_bytes()), "source_manifest_sha256": sha(manifest.read_bytes()), "fasta_sha256": sha(raw_fasta), "sequence_index_sha256": sha((out / "sequence-index.csv").read_bytes()), "panel_records": 2000, "fasta_order": "lexicographic accession", "purpose": "Read-only similarity clustering; no design optimization or new candidate measurements."}
    (out / "input-provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


def cluster(out, binary, lock, threads):
    provenance = json.loads((out / "input-provenance.json").read_text())
    raw = gzip.decompress((out / "panel.fasta.gz").read_bytes())
    if sha(raw) != provenance["fasta_sha256"]:
        raise ValueError("Clustering input drift")
    if (out / "clusters.tsv").exists():
        raise ValueError("Clustering output already exists; use a fresh output directory")
    binary = binary.resolve()
    version = subprocess.check_output([str(binary), "version"], text=True).strip()
    command_options = ["--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "2", "--single-step-clustering", "1", "--alignment-mode", "3", "--seq-id-mode", "0", "-e", "0.001", "-s", "7.5", "--max-seqs", "3000", "--threads", str(threads), "--split-memory-limit", "4G", "--dbtype", "1", "--createdb-mode", "0", "--remove-tmp-files", "1"]
    lock.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(version=version, binary_sha256=sha(binary.read_bytes()), script_sha256=sha(Path(__file__).read_bytes()), source_url="https://github.com/soedinglab/MMseqs2", documentation="https://github.com/soedinglab/MMseqs2/wiki", options=command_options, identity_definition="Identical aligned residues divided by alignment length (--alignment-mode 3 --seq-id-mode 0)", coverage="At least 0.8 of both query and target (--cov-mode 0)", clustering="One-step greedy sequence-length clustering (--cluster-mode 2); operational similarity groups, not annotated biological families.", sensitivity_limitations="Heuristic prefilter sensitivity 7.5; default low-complexity masking and composition-bias correction; E-value <=0.001 also required. Missing homology and order/algorithm effects remain possible.", input_fasta_sha256=sha(raw))
    (out / "executed-cluster-driver.py").write_bytes(Path(__file__).read_bytes())
    with lock.open("a") as guard:
        print("Waiting for shared compute lock for read-only MMseqs2 clustering", flush=True)
        fcntl.flock(guard, fcntl.LOCK_EX)
        print("Acquired clustering compute lock", flush=True)
        with tempfile.TemporaryDirectory(prefix="sparsedesign-density-cluster-") as temp:
            temp = Path(temp); fasta = temp / "panel.fasta"; fasta.write_bytes(raw); fasta.chmod(0o444)
            command = [str(binary), "easy-cluster", str(fasta), str(temp / "cluster"), str(temp / "work"), *command_options]
            metadata["command"] = command
            started = time.monotonic()
            with (out / "mmseqs.log").open("w") as log:
                result = subprocess.run(["/usr/bin/time", "-v", "-o", str(out / "time.txt"), *command], stdout=log, stderr=subprocess.STDOUT)
            metadata.update(returncode=result.returncode, wall_seconds=time.monotonic()-started)
            if result.returncode == 0:
                shutil.copyfile(temp / "cluster_cluster.tsv", out / "clusters.tsv")
                metadata["clusters_tsv_sha256"] = sha((out / "clusters.tsv").read_bytes())
            (out / "run-provenance.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
            result.check_returncode()
    rows = [line.split("\t") for line in (out / "clusters.tsv").read_text().splitlines()]
    expected = {row["accession"] for row in csv.DictReader((out / "sequence-index.csv").open())}
    if len(rows) != len(expected) or {x[1] for x in rows} != expected or not {x[0] for x in rows} <= expected:
        raise ValueError("MMseqs2 result does not partition the frozen panel")
    print(json.dumps({"sequences":len(rows), "clusters":len({r[0] for r in rows})}), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--prepare", action="store_true")
    p.add_argument("--mmseqs", type=Path, required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--lock", type=Path, required=True)
    a = p.parse_args(); out = a.out.resolve()
    if a.prepare:
        prepare(a.root.resolve(), out)
    cluster(out, a.mmseqs, a.lock.resolve(), a.threads)


if __name__ == "__main__":
    main()
