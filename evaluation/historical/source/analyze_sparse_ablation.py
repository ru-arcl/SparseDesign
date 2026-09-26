#!/usr/bin/env python3
r"""Validate and summarize paired sparse/dense ablation measurements.

The input must be a CSV produced by ``collect_sparse_ablation.py``.  The
default is deliberately strict: every requested condition must have a
contiguous repetition grid, every block must contain one adjacent sparse and
dense run, every objective must match within the recorded tolerance, and all
successful measurements and provenance fields must be well formed and
consistent.  The collector's deterministic block IDs, mode orientations, and
randomized block ordering are reconstructed and checked.

Because a CSV does not encode dimensions that never started, final reports
should pass the intended dimensions explicitly.  For example::

  python3 research/analyze_sparse_ablation.py results.csv \
      --expected-accessions 2 \
      --expected-lambda 0 --expected-lambda 4 \
      --expected-threads 1 --expected-threads 16 \
      --expected-repeats 5 \
      --output results-summary.md \
      --summary-csv results-summary.csv

Use ``--diagnostic`` only for an in-progress or failed collection.  It permits
failed, mismatched, truncated-final, and missing blocks, reports all such
issues, and excludes them from statistics.  It does not relax schema,
provenance, duplicate-row, or successful-measurement validation.

Wall time means the collector's high-resolution external wall clock.  A
speedup is paired dense/sparse, so values above one favor the sparse kernel.
Confidence intervals are deterministic paired-bootstrap intervals for the
median of the within-block speedup ratios; they quantify repeat-level
resampling uncertainty for each condition, not uncertainty across proteins.
They are suppressed for conditions with fewer than three valid pairs.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
import statistics
import sys
from typing import Iterable, Sequence


SCHEDULE_VERSION = "complete-block-v1"

# Keep this explicit: accepting an adjacent but incompatible collector schema
# would make a publication report appear more reproducible than it is.
CSV_FIELDS = (
    "run_started_utc",
    "schedule_version",
    "schedule_seed",
    "schedule_index",
    "block_id",
    "block_position",
    "repeat",
    "mode",
    "pair_status",
    "pair_error",
    "objective_abs_diff",
    "objective_allowed_diff",
    "accession",
    "description",
    "input_order",
    "input_path",
    "input_file_sha256",
    "sequence_sha256",
    "protein_aa",
    "rna_nt",
    "lattice_nodes",
    "lambda",
    "threads_requested",
    "threads_effective",
    "objective_units",
    "solve_seconds",
    "external_wall_seconds",
    "host_load1_before",
    "host_load5_before",
    "host_load15_before",
    "host_load1_after",
    "host_load5_after",
    "host_load15_after",
    "time_elapsed_seconds",
    "time_user_seconds",
    "time_system_seconds",
    "time_cpu_percent",
    "max_rss_kib",
    "harness_path",
    "harness_sha256",
    "codon_table_path",
    "codon_table_sha256",
    "gnu_time_path",
    "gnu_time_sha256",
    "command_prefix_json",
    "command_prefix_executable",
    "command_prefix_sha256",
    "command_json",
    "working_directory",
    "omp_env_json",
    "environment_overrides_json",
    "objective_abs_tolerance",
    "objective_rel_tolerance",
    "timeout_seconds",
    "collector_sha256",
    "host",
    "platform",
    "python_version",
    "returncode",
    "status",
    "error",
    "harness_stdout",
    "harness_stderr",
    "gnu_time_output",
)

PROVENANCE_FIELDS = (
    "run_started_utc",
    "schedule_version",
    "schedule_seed",
    "harness_path",
    "harness_sha256",
    "codon_table_path",
    "codon_table_sha256",
    "gnu_time_path",
    "gnu_time_sha256",
    "command_prefix_json",
    "command_prefix_executable",
    "command_prefix_sha256",
    "working_directory",
    "omp_env_json",
    "environment_overrides_json",
    "objective_abs_tolerance",
    "objective_rel_tolerance",
    "timeout_seconds",
    "collector_sha256",
    "host",
    "platform",
    "python_version",
)

INPUT_FIELDS = (
    "description",
    "input_order",
    "input_path",
    "input_file_sha256",
    "sequence_sha256",
    "protein_aa",
)

PAIR_FIELDS = (
    "accession",
    "description",
    "input_order",
    "input_path",
    "input_file_sha256",
    "sequence_sha256",
    "protein_aa",
    "lambda",
    "threads_requested",
    "repeat",
)

LOAD_FIELDS = (
    "host_load1_before",
    "host_load5_before",
    "host_load15_before",
    "host_load1_after",
    "host_load5_after",
    "host_load15_after",
)

OMP_ENV_KEYS = (
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "OMP_DYNAMIC",
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OMP_WAIT_POLICY",
    "GOMP_CPU_AFFINITY",
    "KMP_AFFINITY",
)

HARNESS_FIELDS = (
    "mode",
    "lambda",
    "threads_requested",
    "threads_effective",
    "protein_aa",
    "rna_nt",
    "lattice_nodes",
    "objective_units",
    "solve_seconds",
)

HEX20 = re.compile(r"[0-9a-f]{20}")
HEX64 = re.compile(r"[0-9a-f]{64}")


class ValidationError(ValueError):
    """The input cannot support a trustworthy ablation report."""


@dataclass(frozen=True)
class Measurement:
    threads_effective: int
    objective: Decimal
    solve_seconds: float
    external_wall_seconds: float
    time_elapsed_seconds: float
    time_user_seconds: float
    time_system_seconds: float
    time_cpu_percent: float
    max_rss_kib: int
    loads: tuple[float, float, float, float, float, float]


@dataclass(frozen=True)
class RunRow:
    row_number: int
    values: dict[str, str]
    schedule_index: int
    block_id: str
    block_position: int
    repetition: int
    mode: str
    pair_status: str
    accession: str
    sequence_sha256: str
    protein_aa: int
    lambda_text: str
    lambda_value: float
    threads_requested: int
    status: str
    measurement: Measurement | None

    @property
    def logical_key(self) -> tuple[str, str, str, int, int]:
        return (
            self.accession,
            self.sequence_sha256,
            self.lambda_text,
            self.threads_requested,
            self.repetition,
        )

    @property
    def condition_key(self) -> tuple[str, str, int]:
        return (self.accession, self.lambda_text, self.threads_requested)


@dataclass(frozen=True)
class Pair:
    block_id: str
    sparse: RunRow
    dense: RunRow

    @property
    def condition_key(self) -> tuple[str, str, int]:
        return self.sparse.condition_key


@dataclass(frozen=True)
class Snapshot:
    rows: tuple[RunRow, ...]
    valid_pairs: tuple[Pair, ...]
    block_count: int
    excluded_blocks: int
    raw_rows: int
    skipped_truncated: int
    issues: tuple[str, ...]
    provenance: tuple[str, ...]
    csv_sha256: str
    expected_accessions: tuple[str, ...]
    expected_lambdas: tuple[str, ...]
    expected_threads: tuple[int, ...]
    expected_repeats: int
    dimensions_inferred: tuple[str, ...]


@dataclass(frozen=True)
class SpeedSummary:
    median: float
    q1: float
    q3: float
    minimum: float
    maximum: float
    ci_low: float | None
    ci_high: float | None
    bootstrap_valid: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("csv", type=Path, help="collector CSV to analyze read-only")
    parser.add_argument(
        "--output", type=Path, help="write Markdown here instead of standard output"
    )
    parser.add_argument(
        "--summary-csv", type=Path, help="also write one compact CSV row per condition"
    )
    parser.add_argument(
        "--diagnostic",
        action="store_true",
        help="report and exclude failed/incomplete blocks instead of refusing them",
    )
    parser.add_argument(
        "--expected-accessions",
        type=int,
        help="intended number of input accessions (recommended for final reports)",
    )
    parser.add_argument(
        "--expected-accession",
        dest="expected_accession_names",
        action="append",
        help="intended accession; repeat to validate the exact accession set",
    )
    parser.add_argument(
        "--expected-lambda",
        dest="expected_lambdas",
        action="append",
        help="intended lambda; repeat for every objective",
    )
    parser.add_argument(
        "--expected-threads",
        dest="expected_threads",
        action="append",
        type=int,
        help="intended requested thread count; repeat for every count",
    )
    parser.add_argument(
        "--expected-repeats",
        type=int,
        help="intended number of paired repetitions per condition",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=10000,
        help="paired-bootstrap replicates for median speedups (default: 10000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260714,
        help="base seed for deterministic bootstrap streams (default: 20260714)",
    )
    args = parser.parse_args()

    if args.expected_accessions is not None and args.expected_accessions < 1:
        parser.error("--expected-accessions must be positive")
    if args.expected_repeats is not None and args.expected_repeats < 1:
        parser.error("--expected-repeats must be positive")
    if args.bootstrap < 0:
        parser.error("--bootstrap must be non-negative")
    if args.expected_threads is not None:
        if any(value < 1 for value in args.expected_threads):
            parser.error("--expected-threads values must be positive")
        if len(set(args.expected_threads)) != len(args.expected_threads):
            parser.error("--expected-threads values must be distinct")
    if args.expected_accession_names is not None:
        if any(not value.strip() for value in args.expected_accession_names):
            parser.error("--expected-accession must not be empty")
        if len(set(args.expected_accession_names)) != len(args.expected_accession_names):
            parser.error("--expected-accession values must be distinct")
        if (
            args.expected_accessions is not None
            and len(args.expected_accession_names) != args.expected_accessions
        ):
            parser.error(
                "--expected-accessions disagrees with the number of "
                "--expected-accession values"
            )
    if args.expected_lambdas is not None:
        normalized: list[str] = []
        for text in args.expected_lambdas:
            normalized.append(canonical_lambda_text(text, "--expected-lambda"))
        if len(set(normalized)) != len(normalized):
            parser.error("--expected-lambda values must be distinct")
        args.expected_lambdas = normalized

    destinations = [path for path in (args.output, args.summary_csv) if path is not None]
    resolved_input = args.csv.resolve()
    resolved_destinations = [path.resolve() for path in destinations]
    if resolved_input in resolved_destinations:
        parser.error("output paths must not overwrite the input CSV")
    if len(resolved_destinations) != len(set(resolved_destinations)):
        parser.error("--output and --summary-csv must be different files")
    return args


def canonical_lambda_text(text: str, option: str = "lambda") -> str:
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(f"{option} is not a number: {text!r}") from error
    if (
        not math.isfinite(value)
        or value < 0.0
        or value > sys.float_info.max / 100.0
    ):
        raise ValueError(f"{option} must be finite and nonnegative: {text!r}")
    return "0" if value == 0.0 else format(value, ".17g")


def fail(row_number: int, message: str) -> ValidationError:
    return ValidationError(f"CSV row {row_number}: {message}")


def required_text(row: dict[str, str], name: str, row_number: int) -> str:
    value = row.get(name)
    if value is None or not value.strip():
        raise fail(row_number, f"{name} is empty")
    return value.strip()


def parse_integer(
    row: dict[str, str], name: str, row_number: int, minimum: int | None = None
) -> int:
    text = required_text(row, name, row_number)
    try:
        value = int(text, 10)
    except ValueError as error:
        raise fail(row_number, f"{name} is not an integer: {text!r}") from error
    if minimum is not None and value < minimum:
        raise fail(row_number, f"{name} must be at least {minimum}, got {value}")
    return value


def parse_float(
    row: dict[str, str],
    name: str,
    row_number: int,
    minimum: float | None = None,
    strictly_positive: bool = False,
) -> float:
    text = required_text(row, name, row_number)
    try:
        value = float(text)
    except ValueError as error:
        raise fail(row_number, f"{name} is not a number: {text!r}") from error
    if not math.isfinite(value):
        raise fail(row_number, f"{name} is not finite: {text!r}")
    if minimum is not None and value < minimum:
        raise fail(row_number, f"{name} must be at least {minimum:g}, got {text!r}")
    if strictly_positive and value <= 0.0:
        raise fail(row_number, f"{name} must be positive, got {text!r}")
    return value


def parse_decimal(
    row: dict[str, str], name: str, row_number: int, nonnegative: bool = False
) -> Decimal:
    text = required_text(row, name, row_number)
    try:
        value = Decimal(text)
    except InvalidOperation as error:
        raise fail(row_number, f"{name} is not a decimal: {text!r}") from error
    if not value.is_finite() or (nonnegative and value < 0):
        qualifier = "finite and nonnegative" if nonnegative else "finite"
        raise fail(row_number, f"{name} must be {qualifier}: {text!r}")
    return value


def parse_json(
    row: dict[str, str], name: str, row_number: int, expected_type: type
) -> object:
    text = required_text(row, name, row_number)
    try:
        value = json.loads(text)
    except json.JSONDecodeError as error:
        raise fail(row_number, f"{name} is not valid JSON: {error.msg}") from error
    if not isinstance(value, expected_type):
        raise fail(row_number, f"{name} must encode a {expected_type.__name__}")
    return value


def validate_sha(row: dict[str, str], name: str, row_number: int) -> str:
    value = required_text(row, name, row_number)
    if not HEX64.fullmatch(value):
        raise fail(row_number, f"{name} must be a lowercase SHA-256 digest")
    return value


def validate_provenance_row(row: dict[str, str], row_number: int) -> tuple[str, ...]:
    timestamp = required_text(row, "run_started_utc", row_number)
    try:
        parsed_timestamp = datetime.fromisoformat(timestamp)
    except ValueError as error:
        raise fail(row_number, "run_started_utc is not an ISO-8601 timestamp") from error
    if parsed_timestamp.tzinfo is None:
        raise fail(row_number, "run_started_utc must include a UTC offset")
    if parsed_timestamp.utcoffset() is None or parsed_timestamp.utcoffset().total_seconds() != 0:
        raise fail(row_number, "run_started_utc must use a UTC offset of +00:00")
    if required_text(row, "schedule_version", row_number) != SCHEDULE_VERSION:
        raise fail(
            row_number,
            f"unsupported schedule_version (expected {SCHEDULE_VERSION})",
        )
    parse_integer(row, "schedule_seed", row_number)
    for name in (
        "harness_sha256",
        "codon_table_sha256",
        "gnu_time_sha256",
        "collector_sha256",
    ):
        validate_sha(row, name, row_number)
    for name in (
        "harness_path",
        "codon_table_path",
        "gnu_time_path",
        "working_directory",
        "host",
        "platform",
        "python_version",
    ):
        required_text(row, name, row_number)

    prefix = parse_json(row, "command_prefix_json", row_number, list)
    if not all(isinstance(item, str) for item in prefix):
        raise fail(row_number, "command_prefix_json must contain only strings")
    prefix_executable = row.get("command_prefix_executable", "").strip()
    prefix_sha = row.get("command_prefix_sha256", "").strip()
    if prefix:
        if not prefix_executable:
            raise fail(row_number, "command_prefix_executable is empty for a nonempty prefix")
        if not HEX64.fullmatch(prefix_sha):
            raise fail(row_number, "command_prefix_sha256 must be set for a nonempty prefix")
    elif prefix_executable or prefix_sha:
        raise fail(row_number, "empty command prefix has executable/hash provenance")

    environment_values: dict[str, dict[str, str]] = {}
    for name in ("omp_env_json", "environment_overrides_json"):
        raw_values = parse_json(row, name, row_number, dict)
        if not all(
            isinstance(key, str) and isinstance(value, str)
            for key, value in raw_values.items()
        ):
            raise fail(row_number, f"{name} must map strings to strings")
        environment_values[name] = raw_values  # type: ignore[assignment]
    omp_values = environment_values["omp_env_json"]
    if set(omp_values) != set(OMP_ENV_KEYS):
        raise fail(row_number, "omp_env_json does not contain the collector's exact key set")
    overrides = environment_values["environment_overrides_json"]
    if "LC_ALL" not in overrides:
        raise fail(row_number, "environment_overrides_json does not record LC_ALL")
    for name in OMP_ENV_KEYS:
        if name in overrides and overrides[name] != omp_values[name]:
            raise fail(
                row_number,
                f"{name} differs between omp_env_json and environment_overrides_json",
            )

    parse_decimal(row, "objective_abs_tolerance", row_number, nonnegative=True)
    parse_decimal(row, "objective_rel_tolerance", row_number, nonnegative=True)
    timeout = row.get("timeout_seconds", "").strip()
    if timeout:
        value = parse_float(row, "timeout_seconds", row_number, strictly_positive=True)
        if not math.isfinite(value):  # Kept explicit for readability of the invariant.
            raise fail(row_number, "timeout_seconds must be finite")
    return tuple(row[name] for name in PROVENANCE_FIELDS)


def validate_command(
    row: dict[str, str],
    row_number: int,
    mode: str,
    schedule_index: int,
    input_order: int,
) -> None:
    command = parse_json(row, "command_json", row_number, list)
    if not command or not all(isinstance(token, str) for token in command):
        raise fail(row_number, "command_json must be a nonempty string argv")
    prefix = parse_json(row, "command_prefix_json", row_number, list)
    assert all(isinstance(token, str) for token in prefix)
    harness_index = 5 + len(prefix)
    if len(command) != harness_index + 10:
        raise fail(row_number, "command_json does not match the collector argv layout")
    if command[:3] != [required_text(row, "gnu_time_path", row_number), "-v", "-o"]:
        raise fail(row_number, "command_json has the wrong GNU time prefix")
    if command[4] != "--":
        raise fail(row_number, "command_json is missing GNU time's '--' separator")
    if command[5:harness_index] != prefix:
        raise fail(row_number, "command_json differs from command_prefix_json")

    expected_suffix = [
        "--mode",
        mode,
        "--lambda",
        required_text(row, "lambda", row_number),
        "--threads",
        required_text(row, "threads_requested", row_number),
        "--codon-table",
    ]
    if command[harness_index + 1 : harness_index + 8] != expected_suffix:
        raise fail(row_number, "command_json has the wrong harness options")

    time_path = Path(command[3])
    harness_path = Path(command[harness_index])
    codon_path = Path(command[harness_index + 8])
    input_path = Path(command[harness_index + 9])
    snapshot_paths = (time_path, harness_path, codon_path, input_path)
    if any(not path.is_absolute() for path in snapshot_paths):
        raise fail(row_number, "command_json snapshot paths must be absolute")
    if len({path.parent for path in snapshot_paths}) != 1:
        raise fail(row_number, "command_json snapshot paths do not share a temporary directory")
    expected_names = (
        f"time-{schedule_index:08d}.txt",
        "sparse_ablation.snapshot",
        "codon_table.snapshot.csv",
        f"input-{input_order:04d}.fasta",
    )
    if tuple(path.name for path in snapshot_paths) != expected_names:
        raise fail(row_number, "command_json snapshot names do not match the collector schedule")


def validate_harness_transcript(row: dict[str, str], row_number: int) -> None:
    lines = [line.strip() for line in row["harness_stdout"].splitlines() if line.strip()]
    if len(lines) != 1:
        raise fail(row_number, "harness_stdout must contain one nonempty result line")
    transcript: dict[str, str] = {}
    for token in lines[0].split():
        name, separator, value = token.partition("=")
        if not separator or not name or not value:
            raise fail(row_number, f"malformed harness_stdout token: {token!r}")
        if name in transcript:
            raise fail(row_number, f"duplicate harness_stdout field: {name}")
        transcript[name] = value
    missing = [name for name in HARNESS_FIELDS if name not in transcript]
    if missing:
        raise fail(
            row_number,
            "harness_stdout is missing fields: " + ", ".join(missing),
        )
    changed = [name for name in HARNESS_FIELDS if transcript[name] != row[name]]
    if changed:
        raise fail(
            row_number,
            "parsed fields differ from harness_stdout: " + ", ".join(changed),
        )


def parse_elapsed_transcript(text: str, row_number: int) -> Decimal:
    components = text.split(":")
    try:
        if len(components) == 2:
            hours = Decimal(0)
            minutes = Decimal(components[0])
            seconds = Decimal(components[1])
        elif len(components) == 3:
            hours = Decimal(components[0])
            minutes = Decimal(components[1])
            seconds = Decimal(components[2])
        else:
            raise fail(row_number, f"unrecognized GNU time elapsed value: {text!r}")
    except InvalidOperation as error:
        raise fail(row_number, f"unrecognized GNU time elapsed value: {text!r}") from error
    total = hours * 3600 + minutes * 60 + seconds
    if (
        not total.is_finite()
        or hours < 0
        or minutes < 0
        or seconds < 0
        or (len(components) == 3 and minutes >= 60)
        or seconds >= 60
    ):
        raise fail(row_number, f"invalid GNU time elapsed value: {text!r}")
    return total


def validate_gnu_time_transcript(row: dict[str, str], row_number: int) -> None:
    output = required_text(row, "gnu_time_output", row_number)
    values: dict[str, str] = {}
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if ": " not in line:
            continue
        name, value = line.split(": ", 1)
        if name in values:
            raise fail(row_number, f"duplicate GNU time field: {name}")
        values[name] = value.strip()

    transcript_fields = {
        "time_user_seconds": "User time (seconds)",
        "time_system_seconds": "System time (seconds)",
        "time_cpu_percent": "Percent of CPU this job got",
        "max_rss_kib": "Maximum resident set size (kbytes)",
        "returncode": "Exit status",
    }
    missing = [label for label in transcript_fields.values() if label not in values]
    elapsed_labels = [
        name for name in values if name.startswith("Elapsed (wall clock) time")
    ]
    if len(elapsed_labels) != 1:
        missing.append("one Elapsed (wall clock) time (...) field")
    if missing:
        raise fail(row_number, "GNU time transcript is missing: " + ", ".join(missing))

    for field, label in transcript_fields.items():
        if row[field].strip() != values[label]:
            raise fail(row_number, f"{field} differs from the GNU time transcript")
    elapsed = parse_elapsed_transcript(values[elapsed_labels[0]], row_number)
    recorded_elapsed = parse_decimal(row, "time_elapsed_seconds", row_number, nonnegative=True)
    if elapsed != recorded_elapsed:
        raise fail(
            row_number,
            f"time_elapsed_seconds={recorded_elapsed} but GNU time reports {elapsed}",
        )


def parse_measurement(row: dict[str, str], row_number: int, requested: int) -> Measurement:
    effective = parse_integer(row, "threads_effective", row_number, minimum=1)
    if effective > requested:
        raise fail(
            row_number,
            f"threads_effective={effective} exceeds threads_requested={requested}",
        )
    objective = parse_decimal(row, "objective_units", row_number)
    solve = parse_float(row, "solve_seconds", row_number, strictly_positive=True)
    external = parse_float(
        row, "external_wall_seconds", row_number, strictly_positive=True
    )
    elapsed = parse_float(row, "time_elapsed_seconds", row_number, minimum=0.0)
    user = parse_float(row, "time_user_seconds", row_number, minimum=0.0)
    system = parse_float(row, "time_system_seconds", row_number, minimum=0.0)
    cpu_text = required_text(row, "time_cpu_percent", row_number)
    if not cpu_text.endswith("%"):
        raise fail(row_number, "time_cpu_percent must end in '%'")
    cpu_row = dict(row)
    cpu_row["time_cpu_percent"] = cpu_text[:-1]
    cpu = parse_float(cpu_row, "time_cpu_percent", row_number, minimum=0.0)
    rss = parse_integer(row, "max_rss_kib", row_number, minimum=1)
    loads = tuple(parse_float(row, name, row_number, minimum=0.0) for name in LOAD_FIELDS)
    if external + 1e-6 < solve:
        raise fail(
            row_number,
            "external_wall_seconds is smaller than solve_seconds",
        )
    # GNU time prints elapsed time to hundredths, so allow one extra centisecond for rounding.
    if external + 0.011 < elapsed:
        raise fail(
            row_number,
            "external_wall_seconds is smaller than GNU time's elapsed measurement",
        )
    if not row.get("harness_stdout", "").strip():
        raise fail(row_number, "harness_stdout is empty for a successful run")
    if row.get("harness_stderr", "").strip():
        raise fail(row_number, "harness_stderr is nonempty for a successful run")
    validate_harness_transcript(row, row_number)
    validate_gnu_time_transcript(row, row_number)
    return Measurement(
        threads_effective=effective,
        objective=objective,
        solve_seconds=solve,
        external_wall_seconds=external,
        time_elapsed_seconds=elapsed,
        time_user_seconds=user,
        time_system_seconds=system,
        time_cpu_percent=cpu,
        max_rss_kib=rss,
        loads=loads,  # type: ignore[arg-type]
    )


def parse_run_row(row: dict[str, str], row_number: int) -> RunRow:
    schedule_index = parse_integer(row, "schedule_index", row_number, minimum=1)
    block_id = required_text(row, "block_id", row_number)
    if not HEX20.fullmatch(block_id):
        raise fail(row_number, "block_id must be 20 lowercase hexadecimal characters")
    position = parse_integer(row, "block_position", row_number, minimum=1)
    if position not in (1, 2):
        raise fail(row_number, "block_position must be 1 or 2")
    repetition = parse_integer(row, "repeat", row_number, minimum=1)
    mode = required_text(row, "mode", row_number)
    if mode not in {"sparse", "dense"}:
        raise fail(row_number, f"unknown mode: {mode!r}")
    pair_status = required_text(row, "pair_status", row_number)
    if pair_status not in {"match", "mismatch", "unresolved"}:
        raise fail(row_number, f"unknown pair_status: {pair_status!r}")

    accession = required_text(row, "accession", row_number)
    required_text(row, "description", row_number)
    input_order = parse_integer(row, "input_order", row_number, minimum=0)
    required_text(row, "input_path", row_number)
    validate_sha(row, "input_file_sha256", row_number)
    sequence_sha = validate_sha(row, "sequence_sha256", row_number)
    protein_aa = parse_integer(row, "protein_aa", row_number, minimum=1)
    lambda_text = required_text(row, "lambda", row_number)
    try:
        canonical = canonical_lambda_text(lambda_text)
    except ValueError as error:
        raise fail(row_number, str(error)) from error
    if lambda_text != canonical:
        raise fail(
            row_number,
            f"lambda is not in collector canonical form ({lambda_text!r} versus {canonical!r})",
        )
    requested = parse_integer(row, "threads_requested", row_number, minimum=1)
    validate_command(row, row_number, mode, schedule_index, input_order)

    status = required_text(row, "status", row_number)
    if status not in {"ok", "error", "timeout"}:
        raise fail(row_number, f"unknown status: {status!r}")
    measurement: Measurement | None = None
    if status == "ok":
        returncode = parse_integer(row, "returncode", row_number)
        if returncode != 0:
            raise fail(row_number, "successful run has a nonzero returncode")
        if row.get("error", "").strip():
            raise fail(row_number, "successful run has a nonempty error field")
        rna_nt = parse_integer(row, "rna_nt", row_number, minimum=1)
        if rna_nt != 3 * protein_aa:
            raise fail(row_number, f"rna_nt={rna_nt} does not equal 3*protein_aa")
        parse_integer(row, "lattice_nodes", row_number, minimum=2)
        measurement = parse_measurement(row, row_number, requested)
    else:
        if not row.get("error", "").strip():
            raise fail(row_number, f"status={status} has an empty error field")
        returncode = row.get("returncode", "").strip()
        if returncode:
            parse_integer(row, "returncode", row_number)
        # These are still expected for every attempted process, even when the
        # harness fails before producing its own measurements.
        parse_float(row, "external_wall_seconds", row_number, minimum=0.0)
        for name in LOAD_FIELDS:
            text = row.get(name, "").strip()
            if text:
                parse_float(row, name, row_number, minimum=0.0)

    return RunRow(
        row_number=row_number,
        values=row,
        schedule_index=schedule_index,
        block_id=block_id,
        block_position=position,
        repetition=repetition,
        mode=mode,
        pair_status=pair_status,
        accession=accession,
        sequence_sha256=sequence_sha,
        protein_aa=protein_aa,
        lambda_text=lambda_text,
        lambda_value=float(lambda_text),
        threads_requested=requested,
        status=status,
        measurement=measurement,
    )


def stable_hash(label: str, fields: Iterable[object]) -> str:
    hasher = hashlib.sha256()
    hasher.update(label.encode("ascii"))
    for field in fields:
        hasher.update(b"\0")
        hasher.update(str(field).encode("utf-8"))
    return hasher.hexdigest()


def expected_block_id(key: tuple[str, str, str, int, int]) -> str:
    return stable_hash("block-id-v1", key)[:20]


def expected_mode_order(row: RunRow, schedule_seed: int) -> tuple[str, str]:
    orientation = stable_hash(
        "mode-orientation-v1",
        (
            schedule_seed,
            row.accession,
            row.sequence_sha256,
            row.lambda_text,
            row.threads_requested,
        ),
    )
    dense_first = int(orientation[-1], 16) % 2
    first_is_dense = bool(dense_first ^ ((row.repetition - 1) % 2))
    return ("dense", "sparse") if first_is_dense else ("sparse", "dense")


def add_issue(issues: list[str], message: str) -> None:
    if message not in issues:
        issues.append(message)


def validate_pair_block(
    block_id: str,
    rows: Sequence[RunRow],
    schedule_seed: int,
    diagnostic: bool,
    issues: list[str],
) -> Pair | None:
    def incomplete(message: str) -> None:
        detail = f"block {block_id}: {message}"
        if diagnostic:
            add_issue(issues, detail)
        else:
            raise ValidationError(detail)

    if len(rows) != 2:
        incomplete(f"contains {len(rows)} rows, expected 2")
        return None
    by_position = {row.block_position: row for row in rows}
    if len(by_position) != 2:
        raise ValidationError(f"block {block_id}: duplicate block_position")
    ordered = [by_position[1], by_position[2]]
    if len({row.mode for row in ordered}) != 2:
        raise ValidationError(f"block {block_id}: duplicate mode")
    if {row.mode for row in ordered} != {"sparse", "dense"}:
        raise ValidationError(f"block {block_id}: does not contain sparse and dense")
    for field in PAIR_FIELDS:
        if ordered[0].values[field] != ordered[1].values[field]:
            raise ValidationError(f"block {block_id}: paired rows differ in {field}")
    if ordered[0].logical_key != ordered[1].logical_key:
        raise ValidationError(f"block {block_id}: paired rows have different logical keys")
    computed_id = expected_block_id(ordered[0].logical_key)
    if block_id != computed_id:
        raise ValidationError(
            f"block {block_id}: ID does not match collector identity ({computed_id})"
        )
    expected_modes = expected_mode_order(ordered[0], schedule_seed)
    observed_modes = tuple(row.mode for row in ordered)
    if observed_modes != expected_modes:
        raise ValidationError(
            f"block {block_id}: mode order {observed_modes} does not match {expected_modes}"
        )
    if ordered[1].schedule_index != ordered[0].schedule_index + 1:
        incomplete("paired rows are not adjacent in the schedule")

    pair_statuses = {row.pair_status for row in ordered}
    pair_errors = {row.values["pair_error"] for row in ordered}
    differences = {row.values["objective_abs_diff"] for row in ordered}
    allowances = {row.values["objective_allowed_diff"] for row in ordered}
    if len(pair_statuses) != 1 or len(pair_errors) != 1:
        raise ValidationError(f"block {block_id}: pair status/error differs between rows")
    if len(differences) != 1 or len(allowances) != 1:
        raise ValidationError(f"block {block_id}: objective comparison differs between rows")
    pair_status = ordered[0].pair_status

    if pair_status == "unresolved":
        if ordered[0].values["objective_abs_diff"].strip() or ordered[0].values[
            "objective_allowed_diff"
        ].strip():
            raise ValidationError(f"block {block_id}: unresolved pair has objective deltas")
        if not ordered[0].values["pair_error"].strip():
            raise ValidationError(f"block {block_id}: unresolved pair has no pair_error")
        incomplete(ordered[0].values["pair_error"].strip())
        return None

    if any(row.status != "ok" or row.measurement is None for row in ordered):
        raise ValidationError(f"block {block_id}: resolved pair contains a failed run")
    sparse = next(row for row in ordered if row.mode == "sparse")
    dense = next(row for row in ordered if row.mode == "dense")
    assert sparse.measurement is not None and dense.measurement is not None
    if sparse.measurement.threads_effective != dense.measurement.threads_effective:
        raise ValidationError(
            f"block {block_id}: effective thread counts differ "
            f"({sparse.measurement.threads_effective} versus "
            f"{dense.measurement.threads_effective})"
        )
    difference = abs(sparse.measurement.objective - dense.measurement.objective)
    absolute_tolerance = parse_decimal(
        sparse.values, "objective_abs_tolerance", sparse.row_number, nonnegative=True
    )
    relative_tolerance = parse_decimal(
        sparse.values, "objective_rel_tolerance", sparse.row_number, nonnegative=True
    )
    scale = max(abs(sparse.measurement.objective), abs(dense.measurement.objective))
    allowed = max(absolute_tolerance, relative_tolerance * scale)
    recorded_difference = parse_decimal(
        sparse.values, "objective_abs_diff", sparse.row_number, nonnegative=True
    )
    recorded_allowed = parse_decimal(
        sparse.values, "objective_allowed_diff", sparse.row_number, nonnegative=True
    )
    if recorded_difference != difference:
        raise ValidationError(
            f"block {block_id}: objective_abs_diff={recorded_difference} but computed {difference}"
        )
    if recorded_allowed != allowed:
        raise ValidationError(
            f"block {block_id}: objective_allowed_diff={recorded_allowed} but computed {allowed}"
        )
    matched = difference <= allowed
    expected_status = "match" if matched else "mismatch"
    if pair_status != expected_status:
        raise ValidationError(
            f"block {block_id}: pair_status={pair_status} but computed {expected_status}"
        )
    if matched and sparse.values["pair_error"].strip():
        raise ValidationError(f"block {block_id}: matched pair has a pair_error")
    if not matched and not sparse.values["pair_error"].strip():
        raise ValidationError(f"block {block_id}: mismatched pair has no pair_error")
    if not matched:
        incomplete(sparse.values["pair_error"].strip())
        return None
    return Pair(block_id=block_id, sparse=sparse, dense=dense)


def resolve_dimensions(
    rows: Sequence[RunRow], args: argparse.Namespace, issues: list[str]
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[int, ...], int, tuple[str, ...]]:
    observed_accessions = {row.accession for row in rows}
    observed_lambdas = {row.lambda_text for row in rows}
    observed_threads = {row.threads_requested for row in rows}
    observed_repeats = {row.repetition for row in rows}
    if not observed_accessions or not observed_lambdas or not observed_threads or not observed_repeats:
        raise ValidationError("CSV contains no complete data rows")

    inferred: list[str] = []
    if args.expected_accession_names is not None:
        accessions = tuple(args.expected_accession_names)
        unexpected = sorted(observed_accessions - set(accessions))
        if unexpected:
            raise ValidationError("unexpected accessions: " + ", ".join(unexpected))
        missing_accessions = sorted(set(accessions) - observed_accessions)
        if missing_accessions:
            add_issue(
                issues,
                "missing expected accessions: " + ", ".join(missing_accessions),
            )
    else:
        accessions = tuple(sorted(observed_accessions))
        inferred.append("accessions")
    if args.expected_accessions is not None:
        if len(observed_accessions) > args.expected_accessions:
            raise ValidationError(
                f"observed {len(observed_accessions)} accessions, more than expected "
                f"{args.expected_accessions}"
            )
        if len(observed_accessions) < args.expected_accessions:
            add_issue(
                issues,
                f"observed {len(observed_accessions)} of {args.expected_accessions} expected accessions",
            )

    if args.expected_lambdas is not None:
        lambdas = tuple(sorted(args.expected_lambdas, key=float))
        unexpected = sorted(observed_lambdas - set(lambdas), key=float)
        if unexpected:
            raise ValidationError("unexpected lambda values: " + ", ".join(unexpected))
    else:
        lambdas = tuple(sorted(observed_lambdas, key=float))
        inferred.append("lambda values")

    if args.expected_threads is not None:
        threads = tuple(sorted(args.expected_threads))
        unexpected_threads = sorted(observed_threads - set(threads))
        if unexpected_threads:
            raise ValidationError(
                "unexpected thread counts: " + ", ".join(map(str, unexpected_threads))
            )
    else:
        threads = tuple(sorted(observed_threads))
        inferred.append("thread counts")

    if args.expected_repeats is not None:
        repeats = args.expected_repeats
        if max(observed_repeats) > repeats:
            raise ValidationError(
                f"observed repetition {max(observed_repeats)}, above expected {repeats}"
            )
    else:
        repeats = max(observed_repeats)
        inferred.append("repeat count")
    invalid_repeats = sorted(value for value in observed_repeats if value > repeats)
    if invalid_repeats:
        raise ValidationError("unexpected repetitions: " + ", ".join(map(str, invalid_repeats)))
    return accessions, lambdas, threads, repeats, tuple(inferred)


def read_snapshot(path: Path, args: argparse.Namespace) -> Snapshot:
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ValidationError(f"cannot read {path}: {error}") from error
    if not raw:
        raise ValidationError(f"{path} is empty")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError(f"{path} is not UTF-8 text") from error
    try:
        table = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except csv.Error as error:
        raise ValidationError(f"cannot parse CSV: {error}") from error
    if not table:
        raise ValidationError(f"{path} has no CSV records")
    header = table[0]
    if len(header) != len(set(header)):
        raise ValidationError("CSV header contains duplicate column names")
    if tuple(header) != CSV_FIELDS:
        missing = [name for name in CSV_FIELDS if name not in header]
        extra = [name for name in header if name not in CSV_FIELDS]
        details: list[str] = []
        if missing:
            details.append("missing: " + ", ".join(missing))
        if extra:
            details.append("unexpected: " + ", ".join(extra))
        if not details:
            details.append("column order differs from the collector schema")
        raise ValidationError("incompatible CSV schema (" + "; ".join(details) + ")")

    issues: list[str] = []
    parsed_rows: list[RunRow] = []
    provenance: tuple[str, ...] | None = None
    skipped_truncated = 0
    nonempty_records = [index for index, values in enumerate(table[1:], start=1) if values]
    last_nonempty = nonempty_records[-1] if nonempty_records else -1
    for index, values in enumerate(table[1:], start=1):
        if not values:
            continue
        row_number = index + 1
        if len(values) != len(header):
            if args.diagnostic and index == last_nonempty:
                skipped_truncated += 1
                add_issue(
                    issues,
                    f"row {row_number}: truncated final CSV record "
                    f"({len(values)} of {len(header)} fields)",
                )
                continue
            raise fail(
                row_number,
                f"has {len(values)} fields but the header has {len(header)}",
            )
        row = dict(zip(header, values))
        row_provenance = validate_provenance_row(row, row_number)
        if provenance is None:
            provenance = row_provenance
        elif row_provenance != provenance:
            changed = [
                name
                for name, before, after in zip(PROVENANCE_FIELDS, provenance, row_provenance)
                if before != after
            ]
            raise fail(
                row_number,
                "run provenance differs from the first row: " + ", ".join(changed),
            )
        parsed_rows.append(parse_run_row(row, row_number))
    if not parsed_rows or provenance is None:
        raise ValidationError("CSV contains no complete data rows")

    schedule_indices: dict[int, int] = {}
    schedule_rows = sorted(parsed_rows, key=lambda row: row.schedule_index)
    for row in parsed_rows:
        previous = schedule_indices.setdefault(row.schedule_index, row.row_number)
        if previous != row.row_number:
            raise ValidationError(
                f"duplicate schedule_index={row.schedule_index} in CSV rows "
                f"{previous} and {row.row_number}"
            )
    expected_indices = list(range(1, len(schedule_rows) + 1))
    observed_indices = [row.schedule_index for row in schedule_rows]
    if observed_indices != expected_indices:
        message = "schedule_index values are not contiguous from 1"
        if args.diagnostic:
            add_issue(issues, message)
        else:
            raise ValidationError(message)
    file_order = [row.schedule_index for row in parsed_rows]
    if file_order != sorted(file_order):
        raise ValidationError("CSV rows are not in schedule_index order")

    input_properties: dict[str, tuple[str, ...]] = {}
    successful_input_properties: dict[str, tuple[str, str]] = {}
    order_to_accession: dict[int, str] = {}
    for row in parsed_rows:
        properties = tuple(row.values[name] for name in INPUT_FIELDS)
        previous = input_properties.setdefault(row.accession, properties)
        if previous != properties:
            raise ValidationError(f"{row.accession}: input provenance changes across rows")
        input_order = int(row.values["input_order"])
        previous_accession = order_to_accession.setdefault(input_order, row.accession)
        if previous_accession != row.accession:
            raise ValidationError(
                f"input_order={input_order} is shared by {previous_accession} and {row.accession}"
            )
        if row.status == "ok":
            successful = (row.values["rna_nt"], row.values["lattice_nodes"])
            previous_successful = successful_input_properties.setdefault(
                row.accession, successful
            )
            if previous_successful != successful:
                raise ValidationError(
                    f"{row.accession}: RNA length or lattice size changes across successful rows"
                )
    observed_orders = sorted(order_to_accession)
    if observed_orders != list(range(len(observed_orders))):
        raise ValidationError("input_order values are not contiguous from 0")

    block_rows: dict[str, list[RunRow]] = {}
    logical_to_block: dict[tuple[str, str, str, int, int], str] = {}
    for row in parsed_rows:
        block_rows.setdefault(row.block_id, []).append(row)
        prior = logical_to_block.setdefault(row.logical_key, row.block_id)
        if prior != row.block_id:
            raise ValidationError(
                f"logical block {row.logical_key} uses both {prior} and {row.block_id}"
            )
    schedule_seed = int(provenance[PROVENANCE_FIELDS.index("schedule_seed")])
    valid_pairs: list[Pair] = []
    for block_id, rows in block_rows.items():
        for row in rows:
            computed_id = expected_block_id(row.logical_key)
            if block_id != computed_id:
                raise ValidationError(
                    f"block {block_id}: ID does not match collector identity ({computed_id})"
                )
        pair = validate_pair_block(
            block_id, rows, schedule_seed, args.diagnostic, issues
        )
        if pair is not None:
            valid_pairs.append(pair)

    block_order = sorted(
        block_rows,
        key=lambda block_id: min(row.schedule_index for row in block_rows[block_id]),
    )
    random_order = sorted(
        block_rows,
        key=lambda block_id: (
            stable_hash(
                "block-order-v1",
                (schedule_seed, *block_rows[block_id][0].logical_key),
            ),
            block_id,
        ),
    )
    if block_order != random_order:
        raise ValidationError("block order does not match the seeded collector schedule")

    accessions, lambdas, threads, repeats, inferred = resolve_dimensions(
        parsed_rows, args, issues
    )
    observed_keys = set(logical_to_block)
    sequence_by_accession = {
        accession: properties[INPUT_FIELDS.index("sequence_sha256")]
        for accession, properties in input_properties.items()
    }
    missing_blocks: list[tuple[str, str, int, int]] = []
    for accession in accessions:
        sequence_sha = sequence_by_accession.get(accession)
        if sequence_sha is None:
            # Exact missing accession identities are only knowable when names
            # were supplied.  Count-only missing accessions are already an issue.
            continue
        for lambda_text in lambdas:
            for requested in threads:
                for repetition in range(1, repeats + 1):
                    key = (accession, sequence_sha, lambda_text, requested, repetition)
                    if key not in observed_keys:
                        missing_blocks.append((accession, lambda_text, requested, repetition))
    if missing_blocks:
        preview = ", ".join(
            f"{accession}/lambda={lambda_text}/j={requested}/r={repetition}"
            for accession, lambda_text, requested, repetition in missing_blocks[:8]
        )
        if len(missing_blocks) > 8:
            preview += f", and {len(missing_blocks) - 8} more"
        add_issue(issues, f"missing {len(missing_blocks)} scheduled blocks: {preview}")

    # An exact DP objective should be independent of mode, repetition, and
    # thread count.  Respect the explicitly recorded tolerances when checking
    # across otherwise valid pairs.
    objectives: dict[tuple[str, str], tuple[Decimal, RunRow]] = {}
    for pair in valid_pairs:
        for row in (pair.sparse, pair.dense):
            assert row.measurement is not None
            key = (row.accession, row.lambda_text)
            previous = objectives.setdefault(key, (row.measurement.objective, row))
            reference, reference_row = previous
            absolute = Decimal(row.values["objective_abs_tolerance"])
            relative = Decimal(row.values["objective_rel_tolerance"])
            allowed = max(absolute, relative * max(abs(reference), abs(row.measurement.objective)))
            if abs(reference - row.measurement.objective) > allowed:
                raise ValidationError(
                    f"objective changes across repetitions/threads for {key}: "
                    f"CSV rows {reference_row.row_number} and {row.row_number}"
                )

    if issues and not args.diagnostic:
        preview = "; ".join(issues[:8])
        if len(issues) > 8:
            preview += f"; and {len(issues) - 8} more"
        raise ValidationError(
            "incomplete or failed ablation collection (use --diagnostic only for "
            f"inspection): {preview}"
        )

    valid_pairs.sort(
        key=lambda pair: (
            pair.sparse.accession,
            pair.sparse.lambda_value,
            pair.sparse.threads_requested,
            pair.sparse.repetition,
        )
    )
    return Snapshot(
        rows=tuple(parsed_rows),
        valid_pairs=tuple(valid_pairs),
        block_count=len(block_rows),
        excluded_blocks=len(block_rows) - len(valid_pairs),
        raw_rows=sum(1 for values in table[1:] if values),
        skipped_truncated=skipped_truncated,
        issues=tuple(issues),
        provenance=provenance,
        csv_sha256=hashlib.sha256(raw).hexdigest(),
        expected_accessions=accessions,
        expected_lambdas=lambdas,
        expected_threads=threads,
        expected_repeats=repeats,
        dimensions_inferred=inferred,
    )


def quantile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("quantile requires at least one value")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("quantile probability must be in [0, 1]")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def paired_speed_summary(
    ratios: Sequence[float], replicates: int, seed: int, stream: str
) -> SpeedSummary:
    if not ratios or any(not math.isfinite(value) or value <= 0 for value in ratios):
        raise ValueError("paired speedups must be finite and positive")
    medians: list[float] = []
    if replicates and len(ratios) >= 3:
        digest = hashlib.sha256(f"{seed}:{stream}".encode("utf-8")).digest()
        rng = random.Random(int.from_bytes(digest[:8], "big"))
        count = len(ratios)
        for _ in range(replicates):
            sample = [ratios[rng.randrange(count)] for _ in range(count)]
            medians.append(float(statistics.median(sample)))
    return SpeedSummary(
        median=float(statistics.median(ratios)),
        q1=quantile(ratios, 0.25),
        q3=quantile(ratios, 0.75),
        minimum=min(ratios),
        maximum=max(ratios),
        ci_low=quantile(medians, 0.025) if medians else None,
        ci_high=quantile(medians, 0.975) if medians else None,
        bootstrap_valid=len(medians),
    )


def group_pairs(pairs: Sequence[Pair]) -> dict[tuple[str, str, int], list[Pair]]:
    result: dict[tuple[str, str, int], list[Pair]] = {}
    for pair in pairs:
        result.setdefault(pair.condition_key, []).append(pair)
    for values in result.values():
        values.sort(key=lambda pair: pair.sparse.repetition)
    return dict(
        sorted(result.items(), key=lambda item: (item[0][0], float(item[0][1]), item[0][2]))
    )


def pair_metrics(
    pairs: Sequence[Pair], bootstrap: int, seed: int, stream_prefix: str
) -> dict[str, object]:
    sparse_measurements = [pair.sparse.measurement for pair in pairs]
    dense_measurements = [pair.dense.measurement for pair in pairs]
    assert all(item is not None for item in sparse_measurements)
    assert all(item is not None for item in dense_measurements)
    sparse = [item for item in sparse_measurements if item is not None]
    dense = [item for item in dense_measurements if item is not None]
    wall_ratios = [
        dense_item.external_wall_seconds / sparse_item.external_wall_seconds
        for sparse_item, dense_item in zip(sparse, dense)
    ]
    solve_ratios = [
        dense_item.solve_seconds / sparse_item.solve_seconds
        for sparse_item, dense_item in zip(sparse, dense)
    ]
    rss_deltas = [
        dense_item.max_rss_kib - sparse_item.max_rss_kib
        for sparse_item, dense_item in zip(sparse, dense)
    ]
    loads = [value for item in (*sparse, *dense) for value in item.loads]
    load1 = [item.loads[index] for item in (*sparse, *dense) for index in (0, 3)]
    load5 = [item.loads[index] for item in (*sparse, *dense) for index in (1, 4)]
    load15 = [item.loads[index] for item in (*sparse, *dense) for index in (2, 5)]
    assert loads
    return {
        "sparse": sparse,
        "dense": dense,
        "wall": paired_speed_summary(wall_ratios, bootstrap, seed, stream_prefix + ":wall"),
        "solve": paired_speed_summary(solve_ratios, bootstrap, seed, stream_prefix + ":solve"),
        "rss_deltas": rss_deltas,
        "load1": load1,
        "load5": load5,
        "load15": load15,
    }


def md_escape(value: object) -> str:
    return str(value).replace("|", r"\|").replace("\n", " ")


def emit_table(lines: list[str], headers: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(md_escape(cell) for cell in row) + " |")


def fmt(value: float, significant: int = 4) -> str:
    return format(value, f".{significant}g")


def fmt_seconds(value: float) -> str:
    if value < 0.01:
        return f"{value:.6f}"
    if value < 10:
        return f"{value:.4f}"
    return f"{value:.3f}"


def fmt_range(values: Sequence[float], suffix: str = "") -> str:
    return f"{fmt(min(values))}–{fmt(max(values))}{suffix}"


def fmt_speed(summary: SpeedSummary) -> str:
    ci = (
        f"{fmt(summary.ci_low)}–{fmt(summary.ci_high)}"
        if summary.ci_low is not None and summary.ci_high is not None
        else "not computed"
    )
    return (
        f"{fmt(summary.median)}× "
        f"[IQR {fmt(summary.q1)}–{fmt(summary.q3)}; "
        f"range {fmt(summary.minimum)}–{fmt(summary.maximum)}; 95% CI {ci}]"
    )


def provenance_map(snapshot: Snapshot) -> dict[str, str]:
    return dict(zip(PROVENANCE_FIELDS, snapshot.provenance))


def build_report(path: Path, snapshot: Snapshot, args: argparse.Namespace) -> str:
    lines: list[str] = ["# Sparse/dense multiloop ablation analysis", ""]
    if args.diagnostic:
        state = "DIAGNOSTIC"
        detail = "issues present; invalid/incomplete blocks excluded" if snapshot.issues else "no issues detected"
    else:
        state = "PASS"
        detail = "strict schema, provenance, schedule, objective, and completeness checks passed"
    lines.append(f"**Validation: {state}.** {detail}.")
    lines.append("")
    lines.append(
        f"Input: `{md_escape(path)}` (SHA-256 `{snapshot.csv_sha256}`). "
        "Wall time is the collector's high-resolution external clock. Paired "
        "speedup is dense/sparse, so values above 1 favor sparse."
    )
    lines.append("")
    lines.append("## Coverage")
    lines.append("")
    coverage_rows = [
        ("CSV data records", f"{snapshot.raw_rows:,}"),
        ("Parsed records", f"{len(snapshot.rows):,}"),
        ("Scheduled blocks observed", f"{snapshot.block_count:,}"),
        ("Valid matched pairs analyzed", f"{len(snapshot.valid_pairs):,}"),
        ("Blocks excluded", f"{snapshot.excluded_blocks:,}"),
        ("Truncated final records skipped", f"{snapshot.skipped_truncated:,}"),
        ("Accessions", f"{len(snapshot.expected_accessions):,}"),
        ("Lambda values", ", ".join(snapshot.expected_lambdas)),
        ("Requested threads", ", ".join(map(str, snapshot.expected_threads))),
        ("Paired repetitions per condition", str(snapshot.expected_repeats)),
    ]
    emit_table(lines, ("Item", "Value"), coverage_rows)
    lines.append("")
    if snapshot.dimensions_inferred:
        lines.append(
            "Completeness caveat: the following dimensions were inferred from "
            "observed rows because no corresponding `--expected-*` option was "
            "provided: " + ", ".join(snapshot.dimensions_inferred) + "."
        )
        lines.append("")

    if snapshot.issues:
        lines.append("## Diagnostic issues")
        lines.append("")
        for issue in snapshot.issues[:100]:
            lines.append(f"- {md_escape(issue)}")
        if len(snapshot.issues) > 100:
            lines.append(f"- … and {len(snapshot.issues) - 100} more issues")
        lines.append("")

    grouped = group_pairs(snapshot.valid_pairs)
    metrics_by_key = {
        key: pair_metrics(
            pairs,
            args.bootstrap,
            args.seed,
            f"{key[0]}:{key[1]}:{key[2]}",
        )
        for key, pairs in grouped.items()
    }
    lines.append("## Paired performance")
    lines.append("")
    if not grouped:
        lines.append("No valid matched pairs are available for analysis.")
        lines.append("")
    else:
        performance_rows: list[tuple[object, ...]] = []
        for key, pairs in grouped.items():
            accession, lambda_text, requested = key
            values = metrics_by_key[key]
            sparse = values["sparse"]
            dense = values["dense"]
            wall = values["wall"]
            solve = values["solve"]
            assert isinstance(sparse, list) and isinstance(dense, list)
            assert isinstance(wall, SpeedSummary) and isinstance(solve, SpeedSummary)
            performance_rows.append(
                (
                    accession,
                    pairs[0].sparse.protein_aa,
                    lambda_text,
                    requested,
                    len(pairs),
                    fmt_seconds(statistics.median(item.external_wall_seconds for item in sparse)),
                    fmt_seconds(statistics.median(item.external_wall_seconds for item in dense)),
                    fmt_speed(wall),
                    fmt_seconds(statistics.median(item.solve_seconds for item in sparse)),
                    fmt_seconds(statistics.median(item.solve_seconds for item in dense)),
                    fmt_speed(solve),
                )
            )
        emit_table(
            lines,
            (
                "Accession",
                "AA",
                "λ",
                "j",
                "N",
                "Sparse wall (s), median",
                "Dense wall (s), median",
                "Paired wall speedup",
                "Sparse solve (s), median",
                "Dense solve (s), median",
                "Paired solve speedup",
            ),
            performance_rows,
        )
        lines.append("")

    lines.append("## Memory, threads, and host load")
    lines.append("")
    if not grouped:
        lines.append("No valid matched pairs are available for analysis.")
        lines.append("")
    else:
        resource_rows: list[tuple[object, ...]] = []
        for key, pairs in grouped.items():
            accession, lambda_text, requested = key
            values = metrics_by_key[key]
            sparse = values["sparse"]
            dense = values["dense"]
            deltas = values["rss_deltas"]
            load1 = values["load1"]
            load5 = values["load5"]
            load15 = values["load15"]
            assert isinstance(sparse, list) and isinstance(dense, list)
            assert isinstance(deltas, list)
            assert isinstance(load1, list) and isinstance(load5, list) and isinstance(load15, list)
            sparse_rss = [item.max_rss_kib / 1024.0 for item in sparse]
            dense_rss = [item.max_rss_kib / 1024.0 for item in dense]
            delta_mib = [value / 1024.0 for value in deltas]
            sparse_threads = [float(item.threads_effective) for item in sparse]
            dense_threads = [float(item.threads_effective) for item in dense]
            resource_rows.append(
                (
                    accession,
                    lambda_text,
                    requested,
                    len(pairs),
                    f"{fmt(statistics.median(sparse_rss))} / {fmt(max(sparse_rss))}",
                    f"{fmt(statistics.median(dense_rss))} / {fmt(max(dense_rss))}",
                    f"{fmt(statistics.median(delta_mib))} [{fmt(min(delta_mib))}–{fmt(max(delta_mib))}]",
                    f"{fmt(statistics.median(sparse_threads))} [{fmt_range(sparse_threads)}]",
                    f"{fmt(statistics.median(dense_threads))} [{fmt_range(dense_threads)}]",
                    fmt_range(load1),
                    fmt_range(load5),
                    fmt_range(load15),
                )
            )
        emit_table(
            lines,
            (
                "Accession",
                "λ",
                "j",
                "N",
                "Sparse RSS MiB, median / max",
                "Dense RSS MiB, median / max",
                "Dense−sparse RSS MiB, median [range]",
                "Sparse effective j, median [range]",
                "Dense effective j, median [range]",
                "Host load1 range",
                "Host load5 range",
                "Host load15 range",
            ),
            resource_rows,
        )
        lines.append("")

    lines.append("## Run provenance")
    lines.append("")
    provenance = provenance_map(snapshot)
    provenance_rows = [(name, provenance[name]) for name in PROVENANCE_FIELDS]
    provenance_rows.extend(
        [
            ("bootstrap_replicates", str(args.bootstrap)),
            ("bootstrap_seed", str(args.seed)),
        ]
    )
    emit_table(lines, ("Field", "Value"), provenance_rows)
    lines.append("")
    lines.append(
        "RSS is GNU time's maximum resident set size. Host-load ranges pool "
        "before/after measurements from both modes within each condition. "
        "Bootstrap resampling keeps each dense/sparse block paired."
    )
    lines.append("")
    return "\n".join(lines)


SUMMARY_FIELDS = (
    "source_csv_sha256",
    "accession",
    "protein_aa",
    "lambda",
    "threads_requested",
    "pairs",
    "sparse_wall_median_seconds",
    "dense_wall_median_seconds",
    "wall_speedup_median",
    "wall_speedup_q1",
    "wall_speedup_q3",
    "wall_speedup_min",
    "wall_speedup_max",
    "wall_speedup_ci95_low",
    "wall_speedup_ci95_high",
    "sparse_solve_median_seconds",
    "dense_solve_median_seconds",
    "solve_speedup_median",
    "solve_speedup_q1",
    "solve_speedup_q3",
    "solve_speedup_min",
    "solve_speedup_max",
    "solve_speedup_ci95_low",
    "solve_speedup_ci95_high",
    "sparse_rss_median_kib",
    "sparse_rss_max_kib",
    "dense_rss_median_kib",
    "dense_rss_max_kib",
    "rss_dense_minus_sparse_median_kib",
    "rss_dense_minus_sparse_min_kib",
    "rss_dense_minus_sparse_max_kib",
    "sparse_effective_threads_median",
    "sparse_effective_threads_min",
    "sparse_effective_threads_max",
    "dense_effective_threads_median",
    "dense_effective_threads_min",
    "dense_effective_threads_max",
    "host_load1_min",
    "host_load1_max",
    "host_load5_min",
    "host_load5_max",
    "host_load15_min",
    "host_load15_max",
    "bootstrap_replicates",
    "bootstrap_seed",
)


def summary_rows(snapshot: Snapshot, args: argparse.Namespace) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, pairs in group_pairs(snapshot.valid_pairs).items():
        accession, lambda_text, requested = key
        values = pair_metrics(
            pairs,
            args.bootstrap,
            args.seed,
            f"{accession}:{lambda_text}:{requested}",
        )
        sparse = values["sparse"]
        dense = values["dense"]
        wall = values["wall"]
        solve = values["solve"]
        deltas = values["rss_deltas"]
        load1 = values["load1"]
        load5 = values["load5"]
        load15 = values["load15"]
        assert isinstance(sparse, list) and isinstance(dense, list)
        assert isinstance(wall, SpeedSummary) and isinstance(solve, SpeedSummary)
        assert isinstance(deltas, list)
        assert isinstance(load1, list) and isinstance(load5, list) and isinstance(load15, list)
        sparse_rss = [item.max_rss_kib for item in sparse]
        dense_rss = [item.max_rss_kib for item in dense]
        sparse_threads = [item.threads_effective for item in sparse]
        dense_threads = [item.threads_effective for item in dense]
        rows.append(
            {
                "source_csv_sha256": snapshot.csv_sha256,
                "accession": accession,
                "protein_aa": pairs[0].sparse.protein_aa,
                "lambda": lambda_text,
                "threads_requested": requested,
                "pairs": len(pairs),
                "sparse_wall_median_seconds": statistics.median(
                    item.external_wall_seconds for item in sparse
                ),
                "dense_wall_median_seconds": statistics.median(
                    item.external_wall_seconds for item in dense
                ),
                "wall_speedup_median": wall.median,
                "wall_speedup_q1": wall.q1,
                "wall_speedup_q3": wall.q3,
                "wall_speedup_min": wall.minimum,
                "wall_speedup_max": wall.maximum,
                "wall_speedup_ci95_low": "" if wall.ci_low is None else wall.ci_low,
                "wall_speedup_ci95_high": "" if wall.ci_high is None else wall.ci_high,
                "sparse_solve_median_seconds": statistics.median(
                    item.solve_seconds for item in sparse
                ),
                "dense_solve_median_seconds": statistics.median(
                    item.solve_seconds for item in dense
                ),
                "solve_speedup_median": solve.median,
                "solve_speedup_q1": solve.q1,
                "solve_speedup_q3": solve.q3,
                "solve_speedup_min": solve.minimum,
                "solve_speedup_max": solve.maximum,
                "solve_speedup_ci95_low": "" if solve.ci_low is None else solve.ci_low,
                "solve_speedup_ci95_high": "" if solve.ci_high is None else solve.ci_high,
                "sparse_rss_median_kib": statistics.median(sparse_rss),
                "sparse_rss_max_kib": max(sparse_rss),
                "dense_rss_median_kib": statistics.median(dense_rss),
                "dense_rss_max_kib": max(dense_rss),
                "rss_dense_minus_sparse_median_kib": statistics.median(deltas),
                "rss_dense_minus_sparse_min_kib": min(deltas),
                "rss_dense_minus_sparse_max_kib": max(deltas),
                "sparse_effective_threads_median": statistics.median(sparse_threads),
                "sparse_effective_threads_min": min(sparse_threads),
                "sparse_effective_threads_max": max(sparse_threads),
                "dense_effective_threads_median": statistics.median(dense_threads),
                "dense_effective_threads_min": min(dense_threads),
                "dense_effective_threads_max": max(dense_threads),
                "host_load1_min": min(load1),
                "host_load1_max": max(load1),
                "host_load5_min": min(load5),
                "host_load5_max": max(load5),
                "host_load15_min": min(load15),
                "host_load15_max": max(load15),
                "bootstrap_replicates": args.bootstrap,
                "bootstrap_seed": args.seed,
            }
        )
    return rows


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")


def write_summary_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    snapshot = read_snapshot(args.csv, args)
    report = build_report(args.csv, snapshot, args)
    if args.output is None:
        sys.stdout.write(report)
    else:
        write_text(args.output, report)
    if args.summary_csv is not None:
        write_summary_csv(args.summary_csv, summary_rows(snapshot, args))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValidationError, csv.Error, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
