#!/usr/bin/env python3
r"""Strictly validate and summarize the packed-layout dystrophin experiment.

The expected archive contains four cost-only runs in the counterbalanced order
``B, D, D, B`` followed by one production-interface D run.  This analyzer is
read-only with respect to that archive.  It verifies the frozen build-record
hash chain, the final per-run checksum manifest, empty stderr transcripts,
run chronology, commands and binding environment, GNU ``time -v`` exit/RSS/
wall fields, and exact equality of the cost-run metadata and objective.

The Markdown report summarizes per-layout medians and clearly named ratios.
For the full CLI run it records transcript and design-field metadata, but it
does not independently validate translation, folding energy, structure
optimality, or CAI.  Use ``--output`` for Markdown and ``--summary-json`` for
the machine-readable form.  No report is written unless all checks pass.

Examples::

  python3 research/analyze_packed_layout.py \
      --output /tmp/packed-layout.md \
      --summary-json /tmp/packed-layout.json

  python3 research/analyze_packed_layout.py --self-test
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import stat
import statistics
import sys
import tempfile
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = (
    REPO_ROOT / "research" / "results" / "packed-layout-2026-07-14"
)

CHECKSUM_LINE = re.compile(
    r"^(?P<digest>[0-9a-fA-F]{64}) (?P<mode>[ *])(?P<name>.+)$"
)
ISO_UTC = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[^\s]+(?:Z|[+-][0-9]{2}:[0-9]{2})")
UPTIME_LOAD = re.compile(
    r"load average: (?P<one>[0-9]+(?:\.[0-9]+)?), "
    r"(?P<five>[0-9]+(?:\.[0-9]+)?), (?P<fifteen>[0-9]+(?:\.[0-9]+)?)"
)

# These values anchor this analyzer to the immutable build snapshot used by
# the 2026-07-14 experiment.  The temporary binaries themselves are not
# retained in the result directory; their hashes are cross-checked against
# every run's pre-execution metadata instead.
EXPECTED_FROZEN_HASHES = {
    "SOURCE_SHA256SUMS": "f63f19747ce97d7ca95a7d4dc93ae5998cafcbdda8d072190f7489cf53128d44",
    "BUILD_FLAGS": "f665c6fdbe5aa9ac96ebc12b3531c195d49eeec102b9f3bce900f7c939663b45",
    "BIN_SHA256SUMS": "762c8cc1d2e0234f070b3e18aa304f1baaa50d273f22242cca9508cbe6e9b2f3",
    "VALIDATION": "89bc2e06e23ab4b456347195a05cb77cccbf4bd1b9e8e0c3ebb73fe28d7ff7b9",
}
EXPECTED_FROZEN_MANIFEST_SHA256 = (
    "3684991e875f7ad0e0a839b56b9769ff6a7a744145d7c084606c13f50ef611c1"
)
EXPECTED_RUN_SCRIPT_SHA256 = (
    "97c707d56e490dfac0afd86d4dbb8bb7276507076c8a3ffac6c131e4a186fc45"
)
EXPECTED_BINDER_SHA256 = (
    "0bbada9fb6e1dda5db30ae218ad0685a371578d80c3c3fae16cac1d73aac8177"
)
EXPECTED_BINDER_NAME = "numabind0-epyc-0bbada9fb6e1"
EXPECTED_INPUT_SHA256 = (
    "8bddb035f39b0dc2bb6c0e49656c170ad0da417bc814bd0f568922796d58cbf4"
)
EXPECTED_INPUT_NAME = "NP_000100.3_dystrophin_Dp427c.fasta"
EXPECTED_INPUT_HEADER = ">NP_000100.3 dystrophin isoform Dp427c [Homo sapiens]"
EXPECTED_PROTEIN_AA = 3677
EXPECTED_RNA_NT = 3 * EXPECTED_PROTEIN_AA
EXPECTED_LATTICE_NODES = 12824
EXPECTED_OBJECTIVE_UNITS = Decimal("-716140")
EXPECTED_CODON_SHA256 = (
    "e6c4a49036fe21f3266e7e53130011bd7725c62ebeabbe442f7f8f1cdf5890c4"
)
EXPECTED_CODON_NAME = "codon_usage_freq_table_human.csv"

EXPECTED_SOURCE_NAMES = frozenset(
    {
        "lineardesign-clean/src/main.cc",
        "lineardesign-clean/src/fold_turner.cc",
        "lineardesign-clean/src/fold_turner.h",
        "lineardesign-clean/src/energy.cc",
        "lineardesign-clean/src/energy.h",
        "lineardesign-clean/src/dfa.cc",
        "lineardesign-clean/src/dfa.h",
        "lineardesign-clean/src/codon_table.cc",
        "lineardesign-clean/src/codon_table.h",
        "research/tools/sparse_ablation.cc",
        "lineardesign-clean/data/codon_usage_freq_table_human.csv",
    }
)
EXPECTED_BINARY_NAMES = frozenset(
    {f"lineardesign-{letter}" for letter in "ABCD"}
    | {f"cost-{letter}" for letter in "ABCD"}
)
OMP_ASSIGNMENTS = (
    "OMP_PROC_BIND=close",
    "OMP_PLACES=cores",
    "OMP_DYNAMIC=false",
    "OMP_WAIT_POLICY=active",
    "OMP_THREAD_LIMIT=16",
)
COST_FIELDS = (
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


class ValidationError(ValueError):
    """The archive cannot support a trustworthy packed-layout report."""


@dataclass(frozen=True)
class RunSpec:
    run_id: str
    layout: str
    binary_name: str
    full_cli: bool = False

    @property
    def files(self) -> tuple[str, ...]:
        return tuple(f"{self.run_id}.{suffix}" for suffix in ("out", "err", "time", "meta"))


RUN_SPECS = (
    RunSpec("dmd-01-B-square", "B-square", "cost-B"),
    RunSpec("dmd-02-D-packed", "D-packed", "cost-D"),
    RunSpec("dmd-03-D-packed", "D-packed", "cost-D"),
    RunSpec("dmd-04-B-square", "B-square", "cost-B"),
    RunSpec("dmd-05-D-packed-full", "D-packed", "lineardesign-D", True),
)
EXPECTED_RUN_FILES = frozenset(name for spec in RUN_SPECS for name in spec.files)


@dataclass(frozen=True)
class FrozenBuild:
    file_hashes: dict[str, str]
    binary_hashes: dict[str, str]
    source_hashes: dict[str, str]
    dmd_input_sha256: str


@dataclass(frozen=True)
class Meta:
    started_utc: datetime
    ended_utc: datetime
    loads_before: tuple[float, float, float]
    loads_after: tuple[float, float, float]
    paths: dict[str, str]
    hashes: dict[str, str]


@dataclass(frozen=True)
class Timing:
    command: str
    user_seconds: Decimal
    system_seconds: Decimal
    cpu_percent: Decimal
    elapsed_text: str
    wall_seconds: Decimal
    max_rss_kib: int
    exit_status: int


@dataclass(frozen=True)
class CostResult:
    fields: dict[str, str]
    solve_seconds: Decimal


@dataclass(frozen=True)
class CliResult:
    header: str
    rna_nt: int
    base_pairs: int
    mfe_kcal: Decimal
    cai: Decimal
    sequence_sha256: str
    structure_sha256: str


@dataclass(frozen=True)
class Observation:
    spec: RunSpec
    meta: Meta
    timing: Timing
    output_bytes: int
    output_sha256: str
    cost: CostResult | None
    cli: CliResult | None


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def read_regular(path: Path, label: str, *, allow_empty: bool = False) -> bytes:
    try:
        status = path.lstat()
    except OSError as error:
        raise ValidationError(f"missing or unreadable {label}: {path}: {error}") from error
    if path.is_symlink() or not stat.S_ISREG(status.st_mode):
        raise ValidationError(f"missing or non-regular {label}: {path}")
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ValidationError(f"cannot read {label} {path}: {error}") from error
    if not allow_empty and not raw:
        raise ValidationError(f"{label} is empty: {path}")
    return raw


def decode_text(raw: bytes, label: str, *, require_final_newline: bool = True) -> str:
    try:
        text = raw.decode("utf-8")
    except UnicodeError as error:
        raise ValidationError(f"{label} is not valid UTF-8") from error
    if "\x00" in text:
        raise ValidationError(f"{label} contains a NUL byte")
    if require_final_newline and not text.endswith("\n"):
        raise ValidationError(f"{label} has no final newline (possibly truncated)")
    return text


def parse_checksum_manifest(raw: bytes, label: str) -> dict[str, str]:
    text = decode_text(raw, label)
    entries: dict[str, str] = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line:
            raise ValidationError(f"{label} line {line_number} is blank")
        match = CHECKSUM_LINE.fullmatch(line)
        if match is None:
            raise ValidationError(f"{label} line {line_number} is malformed")
        name = match.group("name")
        if not name or name in entries:
            raise ValidationError(f"{label} has a duplicate/empty entry for {name!r}")
        entries[name] = match.group("digest").lower()
    if not entries:
        raise ValidationError(f"{label} has no checksum entries")
    return entries


def decimal_value(
    text: str, label: str, *, nonnegative: bool = False, positive: bool = False
) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation as error:
        raise ValidationError(f"{label} is not numeric: {text!r}") from error
    if not value.is_finite():
        raise ValidationError(f"{label} must be finite")
    if positive and value <= 0:
        raise ValidationError(f"{label} must be positive")
    if nonnegative and value < 0:
        raise ValidationError(f"{label} must be nonnegative")
    return value


def integer_value(text: str, label: str, *, minimum: int = 0) -> int:
    if re.fullmatch(r"[0-9]+", text) is None:
        raise ValidationError(f"{label} is not an unsigned integer: {text!r}")
    value = int(text)
    if value < minimum:
        raise ValidationError(f"{label} must be at least {minimum}")
    return value


def canonical_decimal(value: Decimal) -> str:
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"


def parse_frozen_build(
    root: Path,
    *,
    expected_hashes: Mapping[str, str] = EXPECTED_FROZEN_HASHES,
    expected_manifest_sha256: str = EXPECTED_FROZEN_MANIFEST_SHA256,
    expected_script_sha256: str = EXPECTED_RUN_SCRIPT_SHA256,
) -> FrozenBuild:
    script_raw = read_regular(root / "run-dystrophin-j16.sh", "run script")
    observed_script = sha256_bytes(script_raw)
    if observed_script != expected_script_sha256:
        raise ValidationError(
            f"run script SHA-256 changed: {observed_script} != {expected_script_sha256}"
        )

    manifest_path = root / "frozen-MANIFEST_SHA256SUMS"
    manifest_raw = read_regular(manifest_path, "frozen manifest")
    observed_manifest_hash = sha256_bytes(manifest_raw)
    if observed_manifest_hash != expected_manifest_sha256:
        raise ValidationError(
            "frozen manifest SHA-256 changed: "
            f"{observed_manifest_hash} != {expected_manifest_sha256}"
        )
    manifest = parse_checksum_manifest(manifest_raw, manifest_path.name)
    if manifest != dict(expected_hashes):
        raise ValidationError("frozen manifest entries differ from the anchored build snapshot")

    retained: dict[str, bytes] = {}
    for original_name, expected_hash in expected_hashes.items():
        path = root / f"frozen-{original_name}"
        raw = read_regular(path, f"frozen {original_name}")
        observed = sha256_bytes(raw)
        if observed != expected_hash:
            raise ValidationError(
                f"{path.name} SHA-256 changed: {observed} != {expected_hash}"
            )
        retained[original_name] = raw

    source = parse_checksum_manifest(retained["SOURCE_SHA256SUMS"], "frozen-SOURCE_SHA256SUMS")
    if frozenset(source) != EXPECTED_SOURCE_NAMES:
        missing = sorted(EXPECTED_SOURCE_NAMES - frozenset(source))
        extra = sorted(frozenset(source) - EXPECTED_SOURCE_NAMES)
        raise ValidationError(f"frozen source manifest grid changed; missing={missing}, extra={extra}")
    codon_source = "lineardesign-clean/data/codon_usage_freq_table_human.csv"
    if source[codon_source] != EXPECTED_CODON_SHA256:
        raise ValidationError("frozen source manifest has the wrong codon-table SHA-256")

    binaries = parse_checksum_manifest(retained["BIN_SHA256SUMS"], "frozen-BIN_SHA256SUMS")
    if frozenset(binaries) != EXPECTED_BINARY_NAMES:
        missing = sorted(EXPECTED_BINARY_NAMES - frozenset(binaries))
        extra = sorted(frozenset(binaries) - EXPECTED_BINARY_NAMES)
        raise ValidationError(f"frozen binary manifest grid changed; missing={missing}, extra={extra}")
    if len(set(binaries.values())) != len(binaries):
        raise ValidationError("frozen binary manifest contains duplicate binary hashes")

    build_text = decode_text(retained["BUILD_FLAGS"], "frozen-BUILD_FLAGS")
    required_build_lines = (
        "common: -std=c++11 -O3 -flto -Wall -fopenmp",
        "A: -DLDCLEAN_SQUARE_MEMOS -DLDCLEAN_NO_SAME_LAYER_PRESEED",
        "B: -DLDCLEAN_SQUARE_MEMOS",
        "C: -DLDCLEAN_PACKED_RANK_MAP",
        "D: (no layout macro; position-major packed default)",
        "CLI sources:",
        "cost source:",
    )
    for marker in required_build_lines:
        if marker not in build_text:
            raise ValidationError(f"frozen-BUILD_FLAGS is missing {marker!r}")

    validation_text = decode_text(retained["VALIDATION"], "frozen-VALIDATION")
    required_validation_markers = (
        "Normal default-D suite: PASS",
        "ASan default-D suite: PASS",
        "UBSan default-D suite: PASS",
        "A/B/C/D CLI output SHA-256:",
        "A/B/C/D cost canonical-metadata SHA-256:",
    )
    for marker in required_validation_markers:
        if marker not in validation_text:
            raise ValidationError(f"frozen-VALIDATION is missing {marker!r}")
    dmd_matches = re.findall(r"^DMD input SHA-256: ([0-9a-f]{64})$", validation_text, re.MULTILINE)
    if len(dmd_matches) != 1:
        raise ValidationError("frozen-VALIDATION must contain exactly one DMD input hash")
    if dmd_matches[0] != EXPECTED_INPUT_SHA256:
        raise ValidationError("frozen-VALIDATION has the wrong DMD input SHA-256")

    return FrozenBuild(dict(expected_hashes), binaries, source, dmd_matches[0])


def parse_utc_timestamp(text: str, label: str) -> datetime:
    if ISO_UTC.fullmatch(text) is None:
        raise ValidationError(f"{label} is not an ISO-8601 timestamp: {text!r}")
    normalized = text[:-1] + "+00:00" if text.endswith("Z") else text
    try:
        value = datetime.fromisoformat(normalized)
    except ValueError as error:
        raise ValidationError(f"{label} is not an ISO-8601 timestamp: {text!r}") from error
    if value.tzinfo is None or value.utcoffset() != timedelta(0):
        raise ValidationError(f"{label} must carry an explicit UTC offset")
    return value


def parse_load(line: str, label: str) -> tuple[float, float, float]:
    match = UPTIME_LOAD.search(line)
    if match is None:
        raise ValidationError(f"{label} has no parseable load averages")
    values = tuple(float(match.group(name)) for name in ("one", "five", "fifteen"))
    if not all(math.isfinite(value) and value >= 0 for value in values):
        raise ValidationError(f"{label} has invalid load averages")
    return values  # type: ignore[return-value]


def parse_meta(raw: bytes, spec: RunSpec, frozen: FrozenBuild) -> Meta:
    label = f"{spec.run_id}.meta"
    text = decode_text(raw, label)
    lines = text.splitlines()
    if len(lines) < 12:
        raise ValidationError(f"{label} is too short (possibly truncated)")
    if ISO_UTC.fullmatch(lines[0]) is None or ISO_UTC.fullmatch(lines[-2]) is None:
        raise ValidationError(f"{label} does not have before/after timestamps at its boundaries")
    started = parse_utc_timestamp(lines[0], f"{label} start")
    ended = parse_utc_timestamp(lines[-2], f"{label} end")
    if ended <= started:
        raise ValidationError(f"{label} end time does not follow its start time")
    before_load = parse_load(lines[1], f"{label} before uptime")
    after_load = parse_load(lines[-1], f"{label} after uptime")

    affinity = [line.strip() for line in lines if "current affinity list:" in line]
    if len(affinity) != 1 or re.fullmatch(
        r"pid [0-9]+'s current affinity list: 0-15", affinity[0]
    ) is None:
        raise ValidationError(f"{label} does not record the expected 0-15 CPU affinity")
    if sum(line.startswith("procs ") for line in lines) != 1:
        raise ValidationError(f"{label} must contain exactly one vmstat transcript")
    columns = [index for index, line in enumerate(lines) if line.strip().startswith("r  b")]
    if len(columns) != 1:
        raise ValidationError(f"{label} has an incomplete vmstat header")
    column_index = columns[0]
    samples = lines[column_index + 1 : column_index + 3]
    if len(samples) != 2 or any(
        len(line.split()) != 17 or any(not value.isdigit() for value in line.split())
        for line in samples
    ):
        raise ValidationError(f"{label} must contain exactly two complete vmstat samples")
    if column_index + 3 >= len(lines) or ISO_UTC.fullmatch(lines[column_index + 3]) is None:
        raise ValidationError(f"{label} has data between vmstat and its ending timestamp")

    checksum_rows: list[tuple[str, str]] = []
    for line in lines:
        match = CHECKSUM_LINE.fullmatch(line)
        if match is not None:
            checksum_rows.append((match.group("name"), match.group("digest").lower()))
    if len(checksum_rows) != 4:
        raise ValidationError(f"{label} must contain exactly four dependency hashes")

    expected_roles = {
        "binary": (spec.binary_name, frozen.binary_hashes[spec.binary_name]),
        "binder": (EXPECTED_BINDER_NAME, EXPECTED_BINDER_SHA256),
        "input": (EXPECTED_INPUT_NAME, frozen.dmd_input_sha256),
        "codon": (EXPECTED_CODON_NAME, EXPECTED_CODON_SHA256),
    }
    paths: dict[str, str] = {}
    hashes: dict[str, str] = {}
    unused = list(checksum_rows)
    for role, (basename, expected_hash) in expected_roles.items():
        matches = [(path, digest) for path, digest in unused if Path(path).name == basename]
        if len(matches) != 1:
            raise ValidationError(f"{label} must hash exactly one {role} named {basename}")
        path, digest = matches[0]
        if digest != expected_hash:
            raise ValidationError(
                f"{label} {role} SHA-256 {digest} differs from frozen {expected_hash}"
            )
        unused.remove(matches[0])
        paths[role] = path
        hashes[role] = digest
    if unused:
        raise ValidationError(f"{label} contains unrecognized dependency hashes")
    return Meta(started, ended, before_load, after_load, paths, hashes)


def elapsed_seconds(text: str, label: str) -> Decimal:
    parts = text.split(":")
    if len(parts) not in (2, 3) or any(not part for part in parts):
        raise ValidationError(f"{label} is malformed: {text!r}")
    if any(re.fullmatch(r"[0-9]+", part) is None for part in parts[:-1]):
        raise ValidationError(f"{label} is malformed: {text!r}")
    seconds = decimal_value(parts[-1], label, nonnegative=True)
    if seconds >= 60:
        raise ValidationError(f"{label} seconds field must be below 60")
    if len(parts) == 2:
        hours = Decimal(0)
        minutes = Decimal(parts[0])
    else:
        hours = Decimal(parts[0])
        minutes = Decimal(parts[1])
        if minutes >= 60:
            raise ValidationError(f"{label} minutes field must be below 60")
    return hours * 3600 + minutes * 60 + seconds


def unique_match(text: str, pattern: str, label: str, filename: str) -> str:
    values = re.findall(pattern, text, re.MULTILINE)
    if len(values) != 1:
        raise ValidationError(
            f"{filename}: expected exactly one GNU-time {label}, observed {len(values)}"
        )
    value = values[0].strip()
    if not value:
        raise ValidationError(f"{filename}: GNU-time {label} is empty")
    return value


def parse_timing(raw: bytes, filename: str) -> Timing:
    text = decode_text(raw, filename)
    command = unique_match(
        text, r"^[ \t]*Command being timed:[ \t]*(.+)$", "command", filename
    )
    user_text = unique_match(
        text, r"^[ \t]*User time \(seconds\):[ \t]*(\S+)[ \t]*$", "user", filename
    )
    system_text = unique_match(
        text, r"^[ \t]*System time \(seconds\):[ \t]*(\S+)[ \t]*$", "system", filename
    )
    cpu_text = unique_match(
        text, r"^[ \t]*Percent of CPU this job got:[ \t]*(\S+)[ \t]*$", "CPU", filename
    )
    elapsed_text = unique_match(
        text,
        r"^[ \t]*Elapsed \(wall clock\) time \([^\n]*\):[ \t]*(\S+)[ \t]*$",
        "elapsed",
        filename,
    )
    rss_text = unique_match(
        text,
        r"^[ \t]*Maximum resident set size \(kbytes\):[ \t]*(\S+)[ \t]*$",
        "maximum RSS",
        filename,
    )
    exit_text = unique_match(
        text, r"^[ \t]*Exit status:[ \t]*(\S+)[ \t]*$", "exit status", filename
    )
    if not text.splitlines()[-1].strip().startswith("Exit status:"):
        raise ValidationError(f"{filename}: exit status is not the final field (possibly truncated)")

    user = decimal_value(user_text, f"{filename} user seconds", nonnegative=True)
    system = decimal_value(system_text, f"{filename} system seconds", nonnegative=True)
    if re.fullmatch(r"[0-9]+%", cpu_text) is None:
        raise ValidationError(f"{filename}: malformed CPU percentage {cpu_text!r}")
    cpu = decimal_value(cpu_text[:-1], f"{filename} CPU percentage", nonnegative=True)
    wall = elapsed_seconds(elapsed_text, f"{filename} elapsed time")
    if wall <= 0:
        raise ValidationError(f"{filename}: elapsed wall time must be positive")
    rss = integer_value(rss_text, f"{filename} maximum RSS", minimum=1)
    exit_status = integer_value(exit_text, f"{filename} exit status")
    if exit_status != 0:
        raise ValidationError(f"{filename}: exit status is {exit_status}, expected 0")
    calculated_cpu = (user + system) * 100 / wall
    # GNU time reports integer CPU percent while wall/user/system are rounded.
    if abs(int(calculated_cpu) - int(cpu)) > 1:
        raise ValidationError(
            f"{filename}: CPU percentage {cpu_text} is inconsistent with user+system/wall "
            f"({canonical_decimal(calculated_cpu)}%)"
        )
    return Timing(command, user, system, cpu, elapsed_text, wall, rss, exit_status)


def command_arguments(command: str, filename: str) -> list[str]:
    rendered = command.strip()
    if len(rendered) >= 2 and rendered[0] == rendered[-1] == '"':
        rendered = rendered[1:-1]
    try:
        arguments = shlex.split(rendered)
    except ValueError as error:
        raise ValidationError(f"{filename}: malformed command rendering") from error
    if not arguments:
        raise ValidationError(f"{filename}: empty command rendering")
    return arguments


def validate_command(spec: RunSpec, timing: Timing, meta: Meta) -> None:
    filename = f"{spec.run_id}.time"
    arguments = command_arguments(timing.command, filename)
    prefix_length = 1 + len(OMP_ASSIGNMENTS)
    if len(arguments) <= prefix_length or arguments[0] != "env":
        raise ValidationError(f"{filename}: command does not begin with env")
    if tuple(arguments[1:prefix_length]) != OMP_ASSIGNMENTS:
        raise ValidationError(f"{filename}: OpenMP environment/order differs from the protocol")
    child = arguments[prefix_length:]
    if len(child) < 2:
        raise ValidationError(f"{filename}: command has no binder/binary pair")
    if child[0] != meta.paths["binder"] or child[1] != meta.paths["binary"]:
        raise ValidationError(f"{filename}: binder/binary paths differ from pre-run metadata")
    solver_args = child[2:]
    if spec.full_cli:
        expected = ["-l", "0", "-j", "16", "-c", meta.paths["codon"]]
    else:
        expected = [
            "--mode",
            "sparse",
            "-l",
            "0",
            "-j",
            "16",
            "-c",
            meta.paths["codon"],
            meta.paths["input"],
        ]
    if solver_args != expected:
        raise ValidationError(f"{filename}: solver arguments differ from the frozen protocol")


def parse_cost_output(raw: bytes, filename: str, timing: Timing) -> CostResult:
    text = decode_text(raw, filename)
    lines = text.splitlines()
    if len(lines) != 1 or not lines[0].strip():
        raise ValidationError(f"{filename}: expected exactly one nonempty result line")
    tokens = lines[0].split()
    parsed: list[tuple[str, str]] = []
    seen: set[str] = set()
    for token in tokens:
        name, separator, value = token.partition("=")
        if not separator or not name or not value or name in seen:
            raise ValidationError(f"{filename}: malformed/duplicate result token {token!r}")
        seen.add(name)
        parsed.append((name, value))
    if tuple(name for name, _ in parsed) != COST_FIELDS:
        raise ValidationError(f"{filename}: cost result schema/order differs from the harness")
    fields = dict(parsed)
    if fields["mode"] != "sparse" or decimal_value(fields["lambda"], f"{filename} lambda") != 0:
        raise ValidationError(f"{filename}: expected sparse mode at lambda zero")
    if integer_value(fields["threads_requested"], f"{filename} requested threads", minimum=1) != 16:
        raise ValidationError(f"{filename}: requested thread count is not 16")
    if integer_value(fields["threads_effective"], f"{filename} effective threads", minimum=1) != 16:
        raise ValidationError(f"{filename}: effective thread count is not 16")
    protein = integer_value(fields["protein_aa"], f"{filename} protein length", minimum=1)
    rna = integer_value(fields["rna_nt"], f"{filename} RNA length", minimum=1)
    nodes = integer_value(fields["lattice_nodes"], f"{filename} lattice nodes", minimum=1)
    if protein != EXPECTED_PROTEIN_AA or rna != EXPECTED_RNA_NT or rna != 3 * protein:
        raise ValidationError(f"{filename}: DMD protein/RNA length metadata is inconsistent")
    if nodes != EXPECTED_LATTICE_NODES:
        raise ValidationError(
            f"{filename}: lattice node count is {nodes}, expected {EXPECTED_LATTICE_NODES}"
        )
    objective = decimal_value(fields["objective_units"], f"{filename} objective")
    if objective != EXPECTED_OBJECTIVE_UNITS:
        raise ValidationError(
            f"{filename}: objective is {objective}, expected {EXPECTED_OBJECTIVE_UNITS}"
        )
    solve = decimal_value(fields["solve_seconds"], f"{filename} solve seconds", positive=True)
    if solve > timing.wall_seconds + Decimal("0.02"):
        raise ValidationError(f"{filename}: solve_seconds exceeds GNU wall time")
    return CostResult(fields, solve)


def validate_dot_bracket(structure: str, filename: str) -> int:
    depth = 0
    pairs = 0
    for index, symbol in enumerate(structure):
        if symbol == "(":
            depth += 1
            pairs += 1
        elif symbol == ")":
            depth -= 1
            if depth < 0:
                raise ValidationError(f"{filename}: structure closes before opening at {index}")
        elif symbol != ".":
            raise ValidationError(f"{filename}: unsupported structure character {symbol!r}")
    if depth:
        raise ValidationError(f"{filename}: structure has {depth} unclosed base pairs")
    return pairs


def parse_cli_output(raw: bytes, filename: str, objective: Decimal) -> CliResult:
    text = decode_text(raw, filename)
    lines = text.splitlines()
    if len(lines) != 4:
        raise ValidationError(f"{filename}: expected exactly four CLI output lines")
    header = lines[0]
    if header != EXPECTED_INPUT_HEADER:
        raise ValidationError(f"{filename}: output header differs from the frozen DMD input")
    sequence_match = re.fullmatch(r"mRNA sequence:[ \t]+([ACGU]+)[ \t]*", lines[1])
    structure_match = re.fullmatch(r"mRNA structure:[ \t]+([().]+)[ \t]*", lines[2])
    summary_match = re.fullmatch(
        r"mRNA folding free energy: ([+-]?[0-9]+\.[0-9]{2}) kcal/mol; "
        r"mRNA CAI: ([0-9]+\.[0-9]{3})",
        lines[3],
    )
    if sequence_match is None or structure_match is None or summary_match is None:
        raise ValidationError(f"{filename}: CLI transcript schema is malformed")
    sequence = sequence_match.group(1)
    structure = structure_match.group(1)
    if len(sequence) != EXPECTED_RNA_NT or len(structure) != len(sequence):
        raise ValidationError(f"{filename}: sequence/structure lengths are inconsistent")
    base_pairs = validate_dot_bracket(structure, filename)
    mfe = decimal_value(summary_match.group(1), f"{filename} printed MFE")
    cai = decimal_value(summary_match.group(2), f"{filename} printed CAI", nonnegative=True)
    if cai > 1:
        raise ValidationError(f"{filename}: printed CAI is above one")
    if mfe * 100 != objective:
        raise ValidationError(
            f"{filename}: printed lambda-zero MFE does not equal the cost-run objective"
        )
    return CliResult(
        header,
        len(sequence),
        base_pairs,
        mfe,
        cai,
        sha256_bytes(sequence.encode("ascii")),
        sha256_bytes(structure.encode("ascii")),
    )


def validate_run_manifest(
    root: Path,
) -> tuple[dict[str, str], dict[str, bytes], str]:
    manifest_path = root / "SHA256SUMS"
    manifest_raw = read_regular(manifest_path, "final run checksum manifest")
    manifest_text = decode_text(manifest_raw, manifest_path.name)
    manifest: dict[str, str] = {}
    resolved_root = root.resolve()
    for line_number, line in enumerate(manifest_text.splitlines(), 1):
        match = CHECKSUM_LINE.fullmatch(line)
        if match is None:
            raise ValidationError(f"{manifest_path.name} line {line_number} is malformed")
        recorded = Path(match.group("name"))
        recorded_absolute = recorded if recorded.is_absolute() else root / recorded
        # Original logs retain their historical absolute location. Resolve only
        # the expected direct filename inside this archive, never that old path.
        historical_parent = Path("/common/home/jy512/dev/rna/research/results/packed-layout-2026-07-14")
        if (".." in recorded.parts or
            (recorded_absolute.parent.resolve() != resolved_root and
             recorded.parent != historical_parent)):
            raise ValidationError(
                f"{manifest_path.name} line {line_number} does not name a direct archive child"
            )
        name = recorded.name
        if not name or name in manifest:
            raise ValidationError(
                f"{manifest_path.name} has a duplicate/empty entry for {name!r}"
            )
        manifest[name] = match.group("digest").lower()
    if frozenset(manifest) != EXPECTED_RUN_FILES:
        missing = sorted(EXPECTED_RUN_FILES - frozenset(manifest))
        extra = sorted(frozenset(manifest) - EXPECTED_RUN_FILES)
        raise ValidationError(f"final run manifest grid differs; missing={missing}, extra={extra}")
    discovered = {path.name for path in root.glob("dmd-*")}
    if discovered != EXPECTED_RUN_FILES:
        missing = sorted(EXPECTED_RUN_FILES - discovered)
        extra = sorted(discovered - EXPECTED_RUN_FILES)
        raise ValidationError(f"retained dmd-* file grid differs; missing={missing}, extra={extra}")
    retained: dict[str, bytes] = {}
    for name, expected in manifest.items():
        raw = read_regular(root / name, f"retained run file {name}", allow_empty=name.endswith(".err"))
        observed = sha256_bytes(raw)
        if observed != expected:
            raise ValidationError(f"{name} SHA-256 changed: {observed} != {expected}")
        retained[name] = raw
    return manifest, retained, sha256_bytes(manifest_raw)


def analyze_archive(
    root: Path,
    *,
    expected_hashes: Mapping[str, str] = EXPECTED_FROZEN_HASHES,
    expected_manifest_sha256: str = EXPECTED_FROZEN_MANIFEST_SHA256,
    expected_script_sha256: str = EXPECTED_RUN_SCRIPT_SHA256,
) -> dict[str, Any]:
    if not root.is_dir():
        raise ValidationError(f"archive is not a directory: {root}")
    frozen = parse_frozen_build(
        root,
        expected_hashes=expected_hashes,
        expected_manifest_sha256=expected_manifest_sha256,
        expected_script_sha256=expected_script_sha256,
    )
    run_manifest, retained, run_manifest_sha256 = validate_run_manifest(root)

    observations: list[Observation] = []
    common_cost_fields: dict[str, str] | None = None
    common_objective: Decimal | None = None
    common_paths: dict[str, str] | None = None
    previous_end: datetime | None = None
    for spec in RUN_SPECS:
        err_raw = retained[f"{spec.run_id}.err"]
        if err_raw:
            raise ValidationError(f"{spec.run_id}.err is nonempty")
        meta = parse_meta(retained[f"{spec.run_id}.meta"], spec, frozen)
        timing = parse_timing(
            retained[f"{spec.run_id}.time"],
            f"{spec.run_id}.time",
        )
        validate_command(spec, timing, meta)
        metadata_window = Decimal(str((meta.ended_utc - meta.started_utc).total_seconds()))
        # The metadata window includes dependency hashing and one second of
        # vmstat setup, while its timestamps have one-second resolution.
        if timing.wall_seconds > metadata_window + Decimal("1.01"):
            raise ValidationError(
                f"{spec.run_id}: GNU wall time exceeds the recorded metadata window"
            )
        if metadata_window > timing.wall_seconds + Decimal("30"):
            raise ValidationError(
                f"{spec.run_id}: metadata window is implausibly longer than GNU wall time"
            )
        if previous_end is not None and meta.started_utc < previous_end:
            raise ValidationError(
                f"{spec.run_id} started before the preceding run ended; B,D,D,B order is not verified"
            )
        previous_end = meta.ended_utc
        shared_paths = {role: meta.paths[role] for role in ("binder", "input", "codon")}
        if common_paths is None:
            common_paths = shared_paths
        elif shared_paths != common_paths:
            raise ValidationError(f"{spec.run_id}.meta dependency paths differ across runs")

        output_raw = retained[f"{spec.run_id}.out"]
        cost: CostResult | None = None
        cli: CliResult | None = None
        if spec.full_cli:
            if common_objective is None:
                raise ValidationError("full CLI run appeared before validated cost runs")
            cli = parse_cli_output(output_raw, f"{spec.run_id}.out", common_objective)
        else:
            cost = parse_cost_output(output_raw, f"{spec.run_id}.out", timing)
            comparable = {name: value for name, value in cost.fields.items() if name != "solve_seconds"}
            objective = decimal_value(cost.fields["objective_units"], "common objective")
            if common_cost_fields is None:
                common_cost_fields = comparable
                common_objective = objective
            elif comparable != common_cost_fields or objective != common_objective:
                raise ValidationError(
                    f"{spec.run_id}.out cost metadata/objective differs from the other cost runs"
                )
        observations.append(
            Observation(
                spec,
                meta,
                timing,
                len(output_raw),
                run_manifest[f"{spec.run_id}.out"],
                cost,
                cli,
            )
        )

    if common_cost_fields is None or common_objective is None or common_paths is None:
        raise ValidationError("no complete cost comparison was found")
    cost_observations = [item for item in observations if item.cost is not None]
    if [item.spec.layout for item in cost_observations] != [
        "B-square",
        "D-packed",
        "D-packed",
        "B-square",
    ]:
        raise ValidationError("cost-run layout order is not B,D,D,B")

    layout_summary: dict[str, dict[str, Any]] = {}
    for layout in ("B-square", "D-packed"):
        selected = [item for item in cost_observations if item.spec.layout == layout]
        if len(selected) != 2:
            raise ValidationError(f"layout {layout} does not have exactly two cost runs")
        solve_values = [float(item.cost.solve_seconds) for item in selected if item.cost]
        wall_values = [float(item.timing.wall_seconds) for item in selected]
        rss_values = [item.timing.max_rss_kib for item in selected]
        layout_summary[layout] = {
            "runs": len(selected),
            "run_ids": [item.spec.run_id for item in selected],
            "median_solve_seconds": statistics.median(solve_values),
            "median_wall_seconds": statistics.median(wall_values),
            "median_max_rss_kib": statistics.median(rss_values),
            "median_max_rss_gib": statistics.median(rss_values) / (1024 * 1024),
        }
    square = layout_summary["B-square"]
    packed = layout_summary["D-packed"]
    ratios = {
        "packed_over_square_solve_time": packed["median_solve_seconds"] / square["median_solve_seconds"],
        "packed_over_square_wall_time": packed["median_wall_seconds"] / square["median_wall_seconds"],
        "packed_over_square_max_rss": packed["median_max_rss_kib"] / square["median_max_rss_kib"],
        "square_over_packed_solve_speedup": square["median_solve_seconds"] / packed["median_solve_seconds"],
        "square_over_packed_wall_speedup": square["median_wall_seconds"] / packed["median_wall_seconds"],
        "median_rss_saving_kib": square["median_max_rss_kib"] - packed["median_max_rss_kib"],
        "median_rss_saving_gib": (
            square["median_max_rss_kib"] - packed["median_max_rss_kib"]
        )
        / (1024 * 1024),
        "median_rss_reduction_fraction": 1 - packed["median_max_rss_kib"] / square["median_max_rss_kib"],
    }
    adjacent_pairs: list[dict[str, Any]] = []
    # Chronological adjacent comparisons: B1/D2 and D3/B4.  The second pair is
    # intentionally reversed in execution order, but ratios always use D/B.
    for pair_number, (square_index, packed_index) in enumerate(((0, 1), (3, 2)), 1):
        square_run = cost_observations[square_index]
        packed_run = cost_observations[packed_index]
        if square_run.cost is None or packed_run.cost is None:
            raise ValidationError("internal pairing selected a non-cost observation")
        square_solve = float(square_run.cost.solve_seconds)
        packed_solve = float(packed_run.cost.solve_seconds)
        square_wall = float(square_run.timing.wall_seconds)
        packed_wall = float(packed_run.timing.wall_seconds)
        square_rss = square_run.timing.max_rss_kib
        packed_rss = packed_run.timing.max_rss_kib
        adjacent_pairs.append(
            {
                "pair": pair_number,
                "square_run": square_run.spec.run_id,
                "packed_run": packed_run.spec.run_id,
                "packed_over_square_solve_time": packed_solve / square_solve,
                "packed_over_square_wall_time": packed_wall / square_wall,
                "packed_over_square_max_rss": packed_rss / square_rss,
                "square_over_packed_solve_speedup": square_solve / packed_solve,
                "rss_saving_kib": square_rss - packed_rss,
            }
        )

    run_rows: list[dict[str, Any]] = []
    for item in observations:
        row: dict[str, Any] = {
            "run_id": item.spec.run_id,
            "layout": item.spec.layout,
            "kind": "full_cli" if item.spec.full_cli else "cost_only",
            "binary_name": item.spec.binary_name,
            "binary_sha256": item.meta.hashes["binary"],
            "started_utc": item.meta.started_utc.isoformat(),
            "ended_utc": item.meta.ended_utc.isoformat(),
            "load_before": list(item.meta.loads_before),
            "load_after": list(item.meta.loads_after),
            "wall_seconds": float(item.timing.wall_seconds),
            "elapsed_text": item.timing.elapsed_text,
            "user_seconds": float(item.timing.user_seconds),
            "system_seconds": float(item.timing.system_seconds),
            "cpu_percent": float(item.timing.cpu_percent),
            "max_rss_kib": item.timing.max_rss_kib,
            "max_rss_gib": item.timing.max_rss_kib / (1024 * 1024),
            "exit_status": item.timing.exit_status,
            "command": item.timing.command,
            "output_bytes": item.output_bytes,
            "output_sha256": item.output_sha256,
        }
        if item.cost is not None:
            row["solve_seconds"] = float(item.cost.solve_seconds)
        if item.cli is not None:
            row["cli"] = {
                "header": item.cli.header,
                "rna_nt": item.cli.rna_nt,
                "base_pairs": item.cli.base_pairs,
                "mfe_kcal_printed": format(item.cli.mfe_kcal, ".2f"),
                "cai_printed": format(item.cli.cai, ".3f"),
                "rna_sha256": item.cli.sequence_sha256,
                "structure_sha256": item.cli.structure_sha256,
            }
        run_rows.append(row)

    full_row = run_rows[-1]
    if full_row["kind"] != "full_cli" or "cli" not in full_row:
        raise ValidationError("the final run is not a complete D full-CLI observation")
    return {
        "schema_version": "packed-layout-summary-v1",
        "status": "pass",
        "archive": str(root.resolve()),
        "protocol": {
            "cost_order": ["B-square", "D-packed", "D-packed", "B-square"],
            "full_run": "D-packed",
            "cost_threads_requested": 16,
            "cost_threads_effective": 16,
            "full_threads_requested": 16,
            "full_threads_effective": None,
            "lambda": "0",
        },
        "frozen": {
            "run_script_sha256": expected_script_sha256,
            "manifest_sha256": expected_manifest_sha256,
            "files": dict(expected_hashes),
            "binaries_used": {
                name: frozen.binary_hashes[name] for name in ("cost-B", "cost-D", "lineardesign-D")
            },
            "binder_sha256": EXPECTED_BINDER_SHA256,
            "input_sha256": frozen.dmd_input_sha256,
            "codon_table_sha256": EXPECTED_CODON_SHA256,
        },
        "run_manifest_sha256": run_manifest_sha256,
        "cost_metadata": common_cost_fields,
        "layouts": layout_summary,
        "ratios": ratios,
        "adjacent_pairs": adjacent_pairs,
        "runs": run_rows,
        "full_cli_run": full_row,
        "validation_scope": {
            "verified": [
                "artifact hashes and completeness",
                "B,D,D,B cost-run chronology and command protocol",
                "empty stderr, zero exit status, GNU-time wall/RSS fields",
                "identical cost metadata and objective",
                "full CLI transcript shape, lengths, printed-MFE/objective consistency, and hashes",
            ],
            "not_independently_verified": [
                "RNA translation to the source protein",
                "folding energy or secondary-structure optimality",
                "CAI recomputation",
                "effective thread count for the full CLI run",
            ],
            "comparison_boundary": (
                "B-square versus D-position-major-packed is a combined layout change; "
                "it does not isolate the packed arena from position-major indexing."
            ),
        },
    }


def markdown_report(summary: Mapping[str, Any]) -> str:
    layouts = summary["layouts"]
    ratios = summary["ratios"]
    metadata = summary["cost_metadata"]
    full = summary["full_cli_run"]
    cli = full["cli"]
    lines = [
        "# Packed-layout dystrophin j16 summary",
        "",
        "Validation status: **PASS**",
        "",
        f"Archive: `{summary['archive']}`",
        "",
        "The four cost runs were verified in the counterbalanced order B, D, D, B; "
        "the full production-interface D run followed them.",
        "",
        "## Cost-only layout comparison",
        "",
        "| Layout | Runs | Median solve (s) | Median wall (s) | Median max RSS (KiB) | Median max RSS (GiB) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for layout in ("B-square", "D-packed"):
        values = layouts[layout]
        lines.append(
            f"| {layout} | {values['runs']} | {values['median_solve_seconds']:.6f} | "
            f"{values['median_wall_seconds']:.2f} | {values['median_max_rss_kib']:.1f} | "
            f"{values['median_max_rss_gib']:.3f} |"
        )
    lines.extend(
        [
            "",
            "Ratios use the medians above; no confidence interval is inferred from two runs per layout.",
            "",
            "| Quantity | Value | Interpretation |",
            "|---|---:|---|",
            f"| Packed / square solve time | {ratios['packed_over_square_solve_time']:.6f} | Below 1 favors packed |",
            f"| Packed / square wall time | {ratios['packed_over_square_wall_time']:.6f} | Below 1 favors packed |",
            f"| Packed / square max RSS | {ratios['packed_over_square_max_rss']:.6f} | Below 1 favors packed |",
            f"| Square / packed solve speedup | {ratios['square_over_packed_solve_speedup']:.6f}x | Above 1 favors packed |",
            f"| Square / packed wall speedup | {ratios['square_over_packed_wall_speedup']:.6f}x | Above 1 favors packed |",
            f"| Median RSS saving | {ratios['median_rss_saving_kib']:.1f} KiB | Square minus packed |",
            f"| Median RSS saving | {ratios['median_rss_saving_gib']:.3f} GiB | Square minus packed |",
            f"| Median RSS reduction | {100 * ratios['median_rss_reduction_fraction']:.3f}% | Relative to square |",
            "",
            "Adjacent counterbalanced pairs (ratios always use packed D / square B):",
            "",
            "| Pair | Square run | Packed run | D/B solve | D/B wall | D/B RSS | RSS saving (KiB) |",
            "|---:|---|---|---:|---:|---:|---:|",
        ]
    )
    for pair in summary["adjacent_pairs"]:
        lines.append(
            f"| {pair['pair']} | `{pair['square_run']}` | `{pair['packed_run']}` | "
            f"{pair['packed_over_square_solve_time']:.6f} | "
            f"{pair['packed_over_square_wall_time']:.6f} | "
            f"{pair['packed_over_square_max_rss']:.6f} | {pair['rss_saving_kib']} |"
        )
    lines.extend(
        [
            "",
            "Common cost metadata:",
            "",
        ]
    )
    for key in COST_FIELDS[:-1]:
        lines.append(f"- `{key}={metadata[key]}`")
    lines.extend(
        [
            "",
            "## Individual runs",
            "",
            "| Order | Run | Layout | Kind | Solve (s) | Wall (s) | Max RSS (KiB) | Load1 before/after |",
            "|---:|---|---|---|---:|---:|---:|---:|",
        ]
    )
    for index, run in enumerate(summary["runs"], 1):
        solve = "—" if "solve_seconds" not in run else f"{run['solve_seconds']:.6f}"
        lines.append(
            f"| {index} | `{run['run_id']}` | {run['layout']} | {run['kind']} | {solve} | "
            f"{run['wall_seconds']:.2f} | {run['max_rss_kib']} | "
            f"{run['load_before'][0]:.2f} / {run['load_after'][0]:.2f} |"
        )
    lines.extend(
        [
            "",
            "## Full production-interface D run",
            "",
            f"- Run: `{full['run_id']}`",
            f"- Binary SHA-256: `{full['binary_sha256']}`",
            f"- Wall time: {full['wall_seconds']:.2f} s",
            f"- Maximum RSS: {full['max_rss_kib']} KiB ({full['max_rss_gib']:.3f} GiB)",
            f"- Output: {full['output_bytes']} bytes, SHA-256 `{full['output_sha256']}`",
            f"- Header: `{cli['header']}`",
            f"- RNA length: {cli['rna_nt']} nt; syntactically balanced base pairs: {cli['base_pairs']}",
            f"- Printed MFE: {cli['mfe_kcal_printed']} kcal/mol; printed CAI: {cli['cai_printed']}",
            f"- RNA SHA-256: `{cli['rna_sha256']}`",
            f"- Structure SHA-256: `{cli['structure_sha256']}`",
            "",
            "## Validation scope",
            "",
            "This report verifies artifact integrity, chronology, commands, resource transcripts, "
            "cost-objective identity, and the full CLI transcript's format and hashes. It does **not** "
            "independently verify RNA translation, folding energy or structure optimality, or CAI; "
            "those require a separate biological/output validator. The full CLI transcript requests "
            "16 threads, but it does not independently report its effective team size.",
            "",
            "The B-to-D comparison is a combined square-to-position-major-packed layout change; it "
            "does not isolate the packed arena from position-major indexing. Peak RSS is whole-process "
            "memory, not a direct measurement of the fixed-cell arena.",
            "",
        ]
    )
    return "\n".join(lines)


def write_text_output(destination: str | None, text: str, label: str) -> None:
    if destination in (None, "-"):
        sys.stdout.write(text)
        return
    path = Path(destination)
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise ValidationError(f"cannot write {label} {path}: {error}") from error


def validate_output_destination(destination: str | None, archive: Path, label: str) -> None:
    if destination in (None, "-"):
        return
    path = Path(destination)
    if path.is_symlink():
        raise ValidationError(f"refusing to write {label} through a symbolic link: {path}")
    resolved = path.resolve()
    archive_root = archive.resolve()
    try:
        resolved.relative_to(archive_root)
    except ValueError:
        return
    raise ValidationError(
        f"refusing to write {label} inside the read-only input archive: {path}"
    )


def synthetic_time(command: Sequence[str], wall: int, rss: int) -> str:
    rendered = shlex.join(list(command))
    user = Decimal(wall) * Decimal("1.5")
    cpu = int(user * 100 / Decimal(wall))
    return (
        f'\tCommand being timed: "{rendered}"\n'
        f"\tUser time (seconds): {user}\n"
        "\tSystem time (seconds): 0\n"
        f"\tPercent of CPU this job got: {cpu}%\n"
        f"\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:{wall:02d}.00\n"
        "\tAverage shared text size (kbytes): 0\n"
        "\tAverage unshared data size (kbytes): 0\n"
        "\tAverage stack size (kbytes): 0\n"
        "\tAverage total size (kbytes): 0\n"
        f"\tMaximum resident set size (kbytes): {rss}\n"
        "\tAverage resident set size (kbytes): 0\n"
        "\tMajor (requiring I/O) page faults: 0\n"
        "\tMinor (reclaiming a frame) page faults: 1\n"
        "\tVoluntary context switches: 1\n"
        "\tInvoluntary context switches: 0\n"
        "\tSwaps: 0\n"
        "\tFile system inputs: 0\n"
        "\tFile system outputs: 0\n"
        "\tSocket messages sent: 0\n"
        "\tSocket messages received: 0\n"
        "\tSignals delivered: 0\n"
        "\tPage size (bytes): 4096\n"
        "\tExit status: 0\n"
    )


def synthetic_meta(
    spec: RunSpec,
    started: datetime,
    ended: datetime,
    binaries: Mapping[str, str],
) -> str:
    binary = f"/tmp/frozen/{spec.binary_name}"
    binder = f"/repo/results/{EXPECTED_BINDER_NAME}"
    protein = f"/repo/data/{EXPECTED_INPUT_NAME}"
    codon = f"/repo/data/{EXPECTED_CODON_NAME}"
    return (
        f"{started.isoformat()}\n"
        " 00:00:00 up 1 day, 1 user, load average: 1.00, 2.00, 3.00\n"
        f"{binaries[spec.binary_name]}  {binary}\n"
        f"{EXPECTED_BINDER_SHA256}  {binder}\n"
        f"{EXPECTED_INPUT_SHA256}  {protein}\n"
        f"{EXPECTED_CODON_SHA256}  {codon}\n"
        "pid 123's current affinity list: 0-15\n"
        "procs -----------memory---------- ---swap-- -----io---- -system-- ------cpu-----\n"
        " r  b   swpd   free   buff  cache   si   so    bi    bo   in   cs us sy id wa st\n"
        " 1  0      0   1000      0   1000    0    0     0     0    1    1 10  0 90  0  0\n"
        " 1  0      0   1000      0   1000    0    0     0     0    1    1 10  0 90  0  0\n"
        f"{ended.isoformat()}\n"
        " 00:00:00 up 1 day, 1 user, load average: 1.50, 2.50, 3.50\n"
    )


def run_self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="packed-layout-self-test-") as temporary:
        root = Path(temporary)
        script = b"#!/bin/sh\n# synthetic packed-layout protocol\n"
        (root / "run-dystrophin-j16.sh").write_bytes(script)

        source_hashes = {name: hashlib.sha256(name.encode()).hexdigest() for name in EXPECTED_SOURCE_NAMES}
        source_hashes["lineardesign-clean/data/codon_usage_freq_table_human.csv"] = EXPECTED_CODON_SHA256
        binary_hashes = {name: hashlib.sha256(("binary:" + name).encode()).hexdigest() for name in EXPECTED_BINARY_NAMES}
        source_text = "".join(f"{digest}  {name}\n" for name, digest in sorted(source_hashes.items()))
        binary_text = "".join(f"{digest}  {name}\n" for name, digest in sorted(binary_hashes.items()))
        build_text = (
            "compiler: synthetic\n"
            "common: -std=c++11 -O3 -flto -Wall -fopenmp\n"
            "A: -DLDCLEAN_SQUARE_MEMOS -DLDCLEAN_NO_SAME_LAYER_PRESEED\n"
            "B: -DLDCLEAN_SQUARE_MEMOS\n"
            "C: -DLDCLEAN_PACKED_RANK_MAP\n"
            "D: (no layout macro; position-major packed default)\n"
            "CLI sources: synthetic\n"
            "cost source: synthetic\n"
        )
        validation_text = (
            "Normal default-D suite: PASS\n"
            "ASan default-D suite: PASS\n"
            "UBSan default-D suite: PASS\n"
            f"A/B/C/D CLI output SHA-256: {'1' * 64}\n"
            f"A/B/C/D cost canonical-metadata SHA-256: {'2' * 64}\n"
            f"DMD input SHA-256: {EXPECTED_INPUT_SHA256}\n"
        )
        frozen_contents = {
            "SOURCE_SHA256SUMS": source_text.encode(),
            "BUILD_FLAGS": build_text.encode(),
            "BIN_SHA256SUMS": binary_text.encode(),
            "VALIDATION": validation_text.encode(),
        }
        frozen_hashes = {name: sha256_bytes(raw) for name, raw in frozen_contents.items()}
        for name, raw in frozen_contents.items():
            (root / f"frozen-{name}").write_bytes(raw)
        frozen_manifest = "".join(f"{frozen_hashes[name]}  {name}\n" for name in frozen_hashes)
        (root / "frozen-MANIFEST_SHA256SUMS").write_text(frozen_manifest, encoding="utf-8")

        walls = [13, 11, 12, 15, 20]
        rss_values = [1000, 600, 620, 1100, 650]
        solves = [Decimal("12"), Decimal("10"), Decimal("11"), Decimal("14")]
        started = datetime(2026, 7, 14, tzinfo=timezone.utc)
        for index, spec in enumerate(RUN_SPECS):
            binary = f"/tmp/frozen/{spec.binary_name}"
            binder = f"/repo/results/{EXPECTED_BINDER_NAME}"
            codon = f"/repo/data/{EXPECTED_CODON_NAME}"
            protein = f"/repo/data/{EXPECTED_INPUT_NAME}"
            command = ["env", *OMP_ASSIGNMENTS, binder, binary]
            if spec.full_cli:
                command.extend(["-l", "0", "-j", "16", "-c", codon])
                sequence = "AUG" * EXPECTED_PROTEIN_AA
                output = (
                    f"{EXPECTED_INPUT_HEADER}\n"
                    f"mRNA sequence:  {sequence}\n"
                    f"mRNA structure: {'.' * EXPECTED_RNA_NT}\n"
                    "mRNA folding free energy: -7161.40 kcal/mol; mRNA CAI: 0.500\n"
                )
            else:
                command.extend(
                    ["--mode", "sparse", "-l", "0", "-j", "16", "-c", codon, protein]
                )
                output = (
                    "mode=sparse lambda=0 threads_requested=16 threads_effective=16 "
                    f"protein_aa={EXPECTED_PROTEIN_AA} rna_nt={EXPECTED_RNA_NT} "
                    f"lattice_nodes={EXPECTED_LATTICE_NODES} "
                    f"objective_units={canonical_decimal(EXPECTED_OBJECTIVE_UNITS)} "
                    f"solve_seconds={solves[index]}\n"
                )
            (root / f"{spec.run_id}.out").write_text(output, encoding="utf-8")
            (root / f"{spec.run_id}.err").write_bytes(b"")
            (root / f"{spec.run_id}.time").write_text(
                synthetic_time(command, walls[index], rss_values[index]), encoding="utf-8"
            )
            ended = started + timedelta(seconds=walls[index] + 2)
            (root / f"{spec.run_id}.meta").write_text(
                synthetic_meta(spec, started, ended, binary_hashes), encoding="utf-8"
            )
            started = ended + timedelta(seconds=1)
        run_manifest = "".join(
            f"{sha256_bytes((root / name).read_bytes())}  {root / name}\n"
            for name in sorted(EXPECTED_RUN_FILES)
        )
        (root / "SHA256SUMS").write_text(run_manifest, encoding="utf-8")

        summary = analyze_archive(
            root,
            expected_hashes=frozen_hashes,
            expected_manifest_sha256=sha256_bytes(frozen_manifest.encode()),
            expected_script_sha256=sha256_bytes(script),
        )
        if summary["status"] != "pass":
            raise AssertionError("synthetic valid archive did not pass")
        if summary["layouts"]["B-square"]["median_solve_seconds"] != 13.0:
            raise AssertionError("square median is wrong")
        if summary["layouts"]["D-packed"]["median_solve_seconds"] != 10.5:
            raise AssertionError("packed median is wrong")
        if "does **not** independently verify" not in markdown_report(summary):
            raise AssertionError("Markdown scope disclaimer is missing")

        # A retained transcript change must fail even when its syntax still looks valid.
        (root / "dmd-02-D-packed.err").write_text("unexpected\n", encoding="utf-8")
        try:
            analyze_archive(
                root,
                expected_hashes=frozen_hashes,
                expected_manifest_sha256=sha256_bytes(frozen_manifest.encode()),
                expected_script_sha256=sha256_bytes(script),
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("tampered synthetic archive was accepted")
    print("packed-layout analyzer self-test: PASS")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "archive",
        nargs="?",
        type=Path,
        default=DEFAULT_ARCHIVE,
        help=f"result directory (default: {DEFAULT_ARCHIVE})",
    )
    parser.add_argument(
        "-o", "--output", metavar="MARKDOWN", help="write Markdown here; use '-' for stdout"
    )
    parser.add_argument(
        "--summary-json", metavar="JSON", help="also write summary JSON here; use '-' for stdout"
    )
    parser.add_argument(
        "--self-test", action="store_true", help="run an internal synthetic valid/tampered-fixture test"
    )
    args = parser.parse_args(argv)
    if args.output == "-" and args.summary_json == "-":
        parser.error("Markdown and JSON cannot both be written to standard output")
    if (
        args.output not in (None, "-")
        and args.summary_json not in (None, "-")
        and Path(args.output).resolve() == Path(args.summary_json).resolve()
    ):
        parser.error("Markdown and JSON output paths must differ")
    if args.self_test and (args.output is not None or args.summary_json is not None):
        parser.error("--self-test cannot be combined with output options")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.self_test:
            run_self_test()
            return 0
        validate_output_destination(args.output, args.archive, "Markdown report")
        validate_output_destination(args.summary_json, args.archive, "JSON summary")
        summary = analyze_archive(args.archive)
        markdown = markdown_report(summary)
        json_text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
        if args.output is None and args.summary_json is None:
            sys.stdout.write(markdown)
        else:
            if args.output is not None:
                write_text_output(args.output, markdown, "Markdown report")
            if args.summary_json is not None:
                write_text_output(args.summary_json, json_text, "JSON summary")
        return 0
    except ValidationError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
