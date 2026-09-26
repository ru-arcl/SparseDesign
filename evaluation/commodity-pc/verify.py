#!/usr/bin/env python3
"""Verify and summarize the frozen ten-run i9 Dp427c subset without a solver."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics


ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def stats(values):
    values = sorted(values)
    require(len(values) == 5, "Expected five observations per layout or pair")
    return {"n": 5, "min": values[0], "q25": values[1],
            "median": statistics.median(values), "q75": values[3],
            "max": values[-1]}


def time_field(text, label):
    match = re.search(r"^\s*" + re.escape(label) + r":\s*(.+)$", text, re.M)
    require(match is not None, "Missing GNU time field: " + label)
    return match.group(1).strip()


def seconds(value):
    result = 0.0
    for field in value.split(":"):
        result = result * 60 + float(field)
    return result


def vmstat(snapshot, key):
    values = dict(line.split() for line in snapshot["vmstat"].splitlines())
    return int(values[key])


def verify():
    manifest = {}
    for line in (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, name = line.split("  ", 1)
        target = (ROOT / name).resolve()
        require(target.is_relative_to(ROOT), "Manifest path escapes dataset")
        require(name not in manifest, "Duplicate manifest entry: " + name)
        require(digest(target) == expected, "SHA-256 mismatch: " + name)
        manifest[name] = expected
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
              if p.is_file() and p.name != "SHA256SUMS"
              and "__pycache__" not in p.parts}
    require(actual == set(manifest), "Manifest does not cover exactly the dataset files")

    origin = read_json(ROOT / "SUBSET_ORIGIN.json")
    for name, info in origin["copied_files"].items():
        require(digest(ROOT / name) == info["sha256"],
                "Original file digest mismatch: " + name)
    protocol = read_json(ROOT / "protocol.json")
    build = read_json(ROOT / "build.json")
    require(digest(ROOT / "preflight-protocol.json") == protocol["preflight_protocol_sha256"],
            "Preflight protocol identity mismatch")
    require(digest(ROOT / "inputs/panel.json") == protocol["panel_sha256"],
            "Original panel identity mismatch")
    for name, expected in build["source_sha256"].items():
        require(digest(ROOT / "source" / name) == expected,
                "Frozen source mismatch: " + name)

    tasks = {t["id"]: t for t in protocol["tasks"] if t["phase"] == "workstation"}
    require(len(tasks) == 10, "Expected ten registered workstation tasks")
    records = sorted((ROOT / "timings/runs").glob("*/record.json"))
    require({p.parent.name for p in records} == set(tasks), "Run set differs from registration")
    observations = {"packed": {}, "square": {}}
    for path in records:
        record = read_json(path)
        task = record["task"]
        require(task == tasks[path.parent.name], "Task registration mismatch")
        require(record["status"] == "ok" and record["exit_code"] == 0,
                "Run did not complete successfully: " + task["id"])
        require(task["accession"] == "NP_000100.3" and task["aa_length"] == 3677
                and task["lambda"] == 0 and task["threads"] == 16
                and task["mode"] == "cli", "Unexpected workstation condition")
        require(digest(ROOT / task["fasta"]) == task["fasta_sha256"], "Input identity mismatch")
        require(record["protocol_sha256"] == digest(ROOT / "protocol.json"),
                "Recorded protocol identity mismatch")
        require(record["collector_sha256"] == digest(ROOT / "source/publication_campaign_final.py"),
                "Recorded collector identity mismatch")
        require(record["codon_table_sha256"] == digest(ROOT / "source/data/codon_usage_freq_table_human.csv"),
                "Recorded codon table identity mismatch")
        require(record["binary_sha256"] == build["binary_sha256"]["sparsedesign-" + task["layout"]],
                "Recorded binary identity mismatch")
        for name, expected in record["output_sha256"].items():
            require(digest(path.parent / name) == expected, "Raw output mismatch: " + task["id"])
        raw_time = (path.parent / "time.txt").read_text(encoding="utf-8")
        elapsed_label = "Elapsed (wall clock) time (h:mm:ss or m:ss)"
        rss_label = "Maximum resident set size (kbytes)"
        for label in [elapsed_label, rss_label, "Exit status"]:
            require(time_field(raw_time, label) == record["gnu_time"][label],
                    "Parsed and raw GNU time differ: " + task["id"])
        require(time_field(raw_time, "Exit status") == "0", "GNU time reports failure")
        require(record["metrics"]["objective_delta"] == 0, "Objective consistency check failed")
        observations[task["layout"]][task["repeat"]] = {
            "wall": seconds(time_field(raw_time, elapsed_label)),
            "rss": int(time_field(raw_time, rss_label)) / (1024 ** 2),
            "swap_in": vmstat(record["after"], "pswpin") - vmstat(record["before"], "pswpin"),
            "swap_out": vmstat(record["after"], "pswpout") - vmstat(record["before"], "pswpout"),
        }

    workstation = []
    for layout, rows in observations.items():
        require(set(rows) == set(range(1, 6)), "Expected repetitions 1 through 5")
        ordered = [rows[rep] for rep in range(1, 6)]
        workstation.append({
            "layout": layout, "expected": 5, "complete": True, "all_successful": True,
            "outcomes": {"missing": 0, "ok": 5}, "successful_repeats": list(range(1, 6)),
            "wall_seconds": stats([r["wall"] for r in ordered]),
            "rss_gib": stats([r["rss"] for r in ordered]),
            "host_swap_in_pages": [r["swap_in"] for r in ordered],
            "host_swap_out_pages": [r["swap_out"] for r in ordered],
        })
    pairs = [{"repeat": rep,
              "packed_over_square_peak_rss": observations["packed"][rep]["rss"] / observations["square"][rep]["rss"],
              "packed_over_square_wall": observations["packed"][rep]["wall"] / observations["square"][rep]["wall"]}
             for rep in range(1, 6)]
    paired = {"complete_pairs": 5, "expected_pairs": 5, "observations": pairs}
    for key in ["packed_over_square_peak_rss", "packed_over_square_wall"]:
        paired[key] = stats([p[key] for p in pairs])
    summary = {"workstation": workstation, "workstation_paired": paired}
    require(summary == read_json(ROOT / "analysis/summary.json"),
            "Recomputed summary differs from recorded workstation summary")
    with (ROOT / "analysis/workstation.csv").open(encoding="utf-8", newline="") as stream:
        csv_rows = list(csv.DictReader(stream))
    require(len(csv_rows) == 2, "Expected two CSV layout rows")
    for row, expected in zip(csv_rows, workstation):
        for key, value in expected.items():
            parsed = row[key] if isinstance(value, str) else (
                row[key] == "True" if isinstance(value, bool) else json.loads(row[key]))
            require(parsed == value, "CSV summary mismatch: " + key)
    return summary, len(manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the recomputed JSON summary outside this dataset")
    args = parser.parse_args()
    summary, count = verify()
    if args.output:
        require(not args.output.resolve().is_relative_to(ROOT), "Choose an output outside the frozen dataset")
        args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("Verified {} file hashes and all 10 workstation records.".format(count))
    for row in summary["workstation"]:
        print("{}: median {:.2f} s, {:.2f} GiB peak RSS (n=5)".format(
            row["layout"], row["wall_seconds"]["median"], row["rss_gib"]["median"]))


if __name__ == "__main__":
    main()
