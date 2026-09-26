#!/usr/bin/env python3
r"""Validate and summarize the archived dystrophin benchmark logs.

The analyzer is read-only with respect to the archive.  It pairs every selected
``.out`` file with its GNU ``time -v`` ``.time`` file, verifies the raw SHA-256
digests against ``SHA256SUMS`` when that manifest is present, validates exactly
one complete design record per solver output, and checks the resource fields in
the time log.  All selected designs must report the same MFE and RNA length.

By default all four filename groups are analyzed:

* ``clean_j*``: primary sparse-solver thread ladder (group ``clean``)
* ``ld_j*``: primary local dense-fork thread ladder (group ``ld``)
* ``rep*``: clean-solver j16 repeats (group ``rep``)
* ``nc_j*_n*_r*``: clean-solver node comparison (group ``nc``)

Each selected group must match its complete archived grid (5 clean, 5 dense,
5 repeat, or 12 node-comparison pairs); missing and unexpected run stems are
fatal even when the checksum manifest was removed or edited.

Repeat ``--group`` to analyze a subset.  The CSV is written to standard output
unless ``--output`` is supplied; ``--markdown`` optionally writes human-readable
tables.  For example:

  python3 research/analyze_dystrophin_bench.py \
      --output /tmp/dystrophin.csv --markdown /tmp/dystrophin.md

  python3 research/analyze_dystrophin_bench.py \
      --group clean --group ld --output /tmp/primary.csv

At lambda zero, CAI is not part of the optimized objective.  The analyzer
therefore records CAI but deliberately does not require the clean and dense
designs to have the same CAI.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import io
from pathlib import Path
import re
import shlex
import statistics
import sys
from typing import Iterable, Sequence, TextIO


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = REPO_ROOT / "research" / "results" / "dystrophin-2026-07-13"
GROUPS = ("clean", "ld", "rep", "nc")
EXPECTED_STEMS = {
    "clean": frozenset(f"clean_j{threads}" for threads in (1, 2, 4, 8, 16)),
    "ld": frozenset(f"ld_j{threads}" for threads in (1, 2, 4, 8, 16)),
    "rep": frozenset(f"rep{repeat}" for repeat in range(1, 6)),
    "nc": frozenset(
        f"nc_j{threads}_n{node}_r{repeat}"
        for threads in (4, 8)
        for node in (0, 1)
        for repeat in range(1, 4)
    ),
}
CAI_NOTE = (
    "At lambda 0, CAI is outside the objective; clean and dense CAI may differ "
    "without an objective disagreement."
)
HISTORICAL_CODON_TABLE = "codon_usage_freq_table_human.csv"

STEM_PATTERNS = (
    ("clean", re.compile(r"^clean_j(?P<threads>[1-9][0-9]*)$")),
    ("ld", re.compile(r"^ld_j(?P<threads>[1-9][0-9]*)$")),
    ("rep", re.compile(r"^rep(?P<repeat>[1-9][0-9]*)$")),
    (
        "nc",
        re.compile(
            r"^nc_j(?P<threads>[1-9][0-9]*)_n(?P<node>[0-9]+)_r"
            r"(?P<repeat>[1-9][0-9]*)$"
        ),
    ),
)
PAIR_FILE = re.compile(r"^(?P<stem>.+)\.(?P<extension>out|time)$")
MANIFEST_LINE = re.compile(r"^(?P<digest>[0-9a-fA-F]{64}) [ *](?P<name>.+)$")
NUMBER = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
SEQUENCE_LINE = re.compile(r"^mRNA sequence:[ \t]*(?P<value>[ACGU]+)[ \t]*$", re.MULTILINE)
STRUCTURE_LINE = re.compile(
    r"^mRNA structure:[ \t]*(?P<value>[().]+)[ \t]*$", re.MULTILINE
)
SUMMARY_LINE = re.compile(
    rf"^mRNA folding free energy:[ \t]*(?P<mfe>{NUMBER})[ \t]+kcal/mol;"
    rf"[ \t]*mRNA CAI:[ \t]*(?P<cai>{NUMBER})[ \t]*$",
    re.MULTILINE,
)
HEADER_LINE = re.compile(r"^>(?P<value>[^\n]+)$", re.MULTILINE)

CSV_FIELDS = (
    "group",
    "run_id",
    "solver",
    "threads",
    "numa_node",
    "repeat",
    "lambda",
    "mfe_kcal",
    "cai",
    "protein_aa",
    "rna_nt",
    "base_pairs",
    "sequence_sha256",
    "structure_sha256",
    "elapsed_text",
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "cpu_percent",
    "max_rss_kib",
    "max_rss_gib",
    "exit_status",
    "command",
    "output_file",
    "output_bytes",
    "output_sha256",
    "time_file",
    "time_bytes",
    "time_sha256",
    "sha256sums_verified",
    "cai_interpretation",
)


class ValidationError(ValueError):
    """An archive or output invariant was violated."""


@dataclass(frozen=True)
class RunSpec:
    group: str
    stem: str
    threads: int
    node: int | None
    repeat: int | None

    @property
    def solver(self) -> str:
        return "local_dense_fork" if self.group == "ld" else "clean_sparse"

    @property
    def sort_key(self) -> tuple[int, int, int, int, str]:
        group_order = GROUPS.index(self.group)
        return (
            group_order,
            self.threads,
            -1 if self.node is None else self.node,
            -1 if self.repeat is None else self.repeat,
            self.stem,
        )


@dataclass(frozen=True)
class Design:
    header: str
    sequence: str
    structure: str
    mfe: Decimal
    mfe_text: str
    cai: Decimal
    cai_text: str
    base_pairs: int


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
class Observation:
    spec: RunSpec
    node: int
    objective_lambda: Decimal
    design: Design
    timing: Timing
    output_path: Path
    output_size: int
    output_sha256: str
    time_path: Path
    time_size: int
    time_sha256: str
    manifest_verified: bool


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "archive",
        nargs="?",
        type=Path,
        default=DEFAULT_ARCHIVE,
        help=f"archive directory (default: {DEFAULT_ARCHIVE})",
    )
    parser.add_argument(
        "--group",
        action="append",
        choices=GROUPS,
        help="filename group to analyze; repeat as needed (default: all groups)",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="CSV",
        help="write machine-readable CSV here instead of stdout; use '-' for stdout",
    )
    parser.add_argument(
        "--markdown",
        metavar="FILE",
        help="also write validated Markdown tables; use '-' for stdout",
    )
    args = parser.parse_args(argv)
    if args.group is None:
        args.groups = GROUPS
    else:
        if len(set(args.group)) != len(args.group):
            parser.error("--group values must not be repeated")
        args.groups = tuple(args.group)
    if args.output in (None, "-") and args.markdown == "-":
        parser.error("CSV and Markdown cannot both be written to standard output")
    return args


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def decimal_value(text: str, context: str, *, nonnegative: bool = False) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation as error:
        raise ValidationError(f"{context} is not a number: {text!r}") from error
    if not value.is_finite() or (nonnegative and value < 0):
        qualifier = "finite and nonnegative" if nonnegative else "finite"
        raise ValidationError(f"{context} must be {qualifier}: {text!r}")
    return value


def canonical_decimal(value: Decimal) -> str:
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"


def parse_manifest(archive: Path) -> dict[str, str] | None:
    path = archive / "SHA256SUMS"
    if not path.exists():
        return None
    if not path.is_file():
        raise ValidationError(f"checksum manifest is not a regular file: {path}")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ValidationError(f"cannot read checksum manifest {path}: {error}") from error
    entries: dict[str, str] = {}
    for line_number, line in enumerate(lines, 1):
        if not line:
            continue
        match = MANIFEST_LINE.fullmatch(line)
        if match is None:
            raise ValidationError(
                f"{path.name} line {line_number}: malformed SHA-256 entry"
            )
        name = match.group("name")
        if name in entries:
            raise ValidationError(
                f"{path.name} line {line_number}: duplicate entry for {name!r}"
            )
        entries[name] = match.group("digest").lower()
    if not entries:
        raise ValidationError(f"checksum manifest is empty: {path}")
    return entries


def parse_stem(stem: str) -> RunSpec | None:
    for group, pattern in STEM_PATTERNS:
        match = pattern.fullmatch(stem)
        if match is None:
            continue
        values = match.groupdict()
        if group == "rep":
            return RunSpec(group, stem, 16, None, int(values["repeat"]))
        return RunSpec(
            group,
            stem,
            int(values["threads"]),
            int(values["node"]) if values.get("node") is not None else None,
            int(values["repeat"]) if values.get("repeat") is not None else None,
        )
    return None


def group_hint(filename: str) -> str | None:
    if filename.startswith("clean_"):
        return "clean"
    if filename.startswith("ld_"):
        return "ld"
    if filename.startswith("rep"):
        return "rep"
    if filename.startswith("nc_"):
        return "nc"
    return None


def discover_pairs(
    archive: Path, groups: Iterable[str], manifest: dict[str, str] | None
) -> list[tuple[RunSpec, Path, Path, bool]]:
    selected = set(groups)
    disk_names: set[str] = set()
    try:
        entries = list(archive.iterdir())
    except OSError as error:
        raise ValidationError(f"cannot list archive {archive}: {error}") from error
    for path in entries:
        if path.is_file() and path.suffix in (".out", ".time"):
            disk_names.add(path.name)

    candidate_names = set(disk_names)
    if manifest is not None:
        candidate_names.update(
            name for name in manifest if Path(name).suffix in (".out", ".time")
        )

    extensions: dict[str, set[str]] = {}
    specs: dict[str, RunSpec] = {}
    for filename in sorted(candidate_names):
        pair_match = PAIR_FILE.fullmatch(filename)
        if pair_match is None:
            continue
        stem = pair_match.group("stem")
        spec = parse_stem(stem)
        if spec is None:
            hinted = group_hint(filename)
            if hinted in selected:
                raise ValidationError(
                    f"selected group {hinted!r} contains unsupported log filename: {filename}"
                )
            continue
        if spec.group not in selected:
            continue
        specs[stem] = spec
        extensions.setdefault(stem, set()).add(pair_match.group("extension"))

    pairs: list[tuple[RunSpec, Path, Path, bool]] = []
    for group in groups:
        group_specs = [spec for spec in specs.values() if spec.group == group]
        observed_stems = {spec.stem for spec in group_specs}
        expected_stems = EXPECTED_STEMS[group]
        if observed_stems != expected_stems:
            missing = sorted(expected_stems - observed_stems)
            unexpected = sorted(observed_stems - expected_stems)
            details: list[str] = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if unexpected:
                details.append("unexpected " + ", ".join(unexpected))
            raise ValidationError(
                f"selected group {group!r} does not match the archived run grid: "
                + "; ".join(details)
            )
        for spec in group_specs:
            present = extensions[spec.stem]
            if present != {"out", "time"}:
                missing = sorted({"out", "time"} - present)
                raise ValidationError(
                    f"{spec.stem}: missing paired {', '.join('.' + item for item in missing)} file"
                )
            output_path = archive / f"{spec.stem}.out"
            time_path = archive / f"{spec.stem}.time"
            for path in (output_path, time_path):
                if not path.is_file():
                    raise ValidationError(f"paired log is not a regular file: {path}")
                if manifest is not None and path.name not in manifest:
                    raise ValidationError(
                        f"selected raw log has no SHA256SUMS entry: {path.name}"
                    )
            pairs.append((spec, output_path, time_path, manifest is not None))
    pairs.sort(key=lambda item: item[0].sort_key)
    return pairs


def read_and_hash(
    path: Path, manifest: dict[str, str] | None
) -> tuple[bytes, str, bool]:
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ValidationError(f"cannot read {path}: {error}") from error
    digest = sha256_bytes(raw)
    verified = manifest is not None
    if manifest is not None:
        expected = manifest.get(path.name)
        if expected is None:
            raise ValidationError(f"{path.name} has no SHA256SUMS entry")
        if digest != expected:
            raise ValidationError(
                f"{path.name}: SHA-256 mismatch (manifest {expected}, observed {digest})"
            )
    return raw, digest, verified


def one_match(pattern: re.Pattern[str], text: str, label: str, filename: str) -> re.Match[str]:
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValidationError(
            f"{filename}: expected exactly one {label}, observed {len(matches)}"
        )
    return matches[0]


def validate_dot_bracket(structure: str, filename: str) -> int:
    depth = 0
    pairs = 0
    for offset, character in enumerate(structure):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            pairs += 1
            if depth < 0:
                raise ValidationError(
                    f"{filename}: structure closes an unopened pair at offset {offset}"
                )
    if depth:
        raise ValidationError(f"{filename}: structure has {depth} unclosed base pairs")
    return pairs


def parse_design(raw: bytes, filename: str) -> Design:
    try:
        text = raw.decode("utf-8")
    except UnicodeError as error:
        raise ValidationError(f"{filename}: solver output is not valid UTF-8") from error
    if "\x00" in text:
        raise ValidationError(f"{filename}: solver output contains a NUL byte")
    # The dense comparator prints progress using carriage returns.  Normalize a
    # parsing copy only; the checksum above always covers the untouched bytes.
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    header_match = one_match(HEADER_LINE, normalized, "FASTA header", filename)
    sequence_match = one_match(SEQUENCE_LINE, normalized, "mRNA sequence line", filename)
    structure_match = one_match(
        STRUCTURE_LINE, normalized, "mRNA structure line", filename
    )
    summary_match = one_match(SUMMARY_LINE, normalized, "design summary line", filename)
    if not (
        header_match.start()
        < sequence_match.start()
        < structure_match.start()
        < summary_match.start()
    ):
        raise ValidationError(f"{filename}: design record fields are out of order")

    header = header_match.group("value")
    if not header.strip() or header != header.strip():
        raise ValidationError(f"{filename}: FASTA header is empty or padded")
    sequence = sequence_match.group("value")
    structure = structure_match.group("value")
    if not sequence:
        raise ValidationError(f"{filename}: mRNA sequence is empty")
    if len(sequence) != len(structure):
        raise ValidationError(
            f"{filename}: sequence length {len(sequence)} != structure length {len(structure)}"
        )
    if len(sequence) % 3:
        raise ValidationError(
            f"{filename}: coding RNA length {len(sequence)} is not divisible by three"
        )
    base_pairs = validate_dot_bracket(structure, filename)
    mfe_text = summary_match.group("mfe")
    cai_text = summary_match.group("cai")
    if re.fullmatch(r"[+-]?[0-9]+\.[0-9]{2}", mfe_text) is None:
        raise ValidationError(f"{filename}: MFE is not printed to two decimal places")
    if re.fullmatch(r"[+-]?[0-9]+\.[0-9]{3}", cai_text) is None:
        raise ValidationError(f"{filename}: CAI is not printed to three decimal places")
    mfe = decimal_value(mfe_text, f"{filename} MFE")
    cai = decimal_value(cai_text, f"{filename} CAI", nonnegative=True)
    if cai > 1:
        raise ValidationError(f"{filename}: CAI must be in [0, 1], observed {cai_text}")
    return Design(
        header,
        sequence,
        structure,
        mfe,
        mfe_text,
        cai,
        cai_text,
        base_pairs,
    )


def elapsed_seconds(text: str, filename: str) -> Decimal:
    parts = text.split(":")
    if len(parts) not in (2, 3) or any(not part for part in parts):
        raise ValidationError(f"{filename}: malformed elapsed time {text!r}")
    if not all(re.fullmatch(r"[0-9]+", part) for part in parts[:-1]):
        raise ValidationError(f"{filename}: malformed elapsed time {text!r}")
    seconds = decimal_value(parts[-1], f"{filename} elapsed seconds", nonnegative=True)
    if seconds >= 60:
        raise ValidationError(f"{filename}: elapsed seconds field must be below 60")
    if len(parts) == 2:
        minutes = Decimal(parts[0])
        hours = Decimal(0)
    else:
        hours = Decimal(parts[0])
        minutes = Decimal(parts[1])
        if minutes >= 60:
            raise ValidationError(f"{filename}: elapsed minutes field must be below 60")
    return hours * 3600 + minutes * 60 + seconds


def unique_time_value(
    text: str, pattern: re.Pattern[str], label: str, filename: str
) -> str:
    matches = [match.group("value").strip() for match in pattern.finditer(text)]
    if len(matches) != 1:
        raise ValidationError(
            f"{filename}: expected exactly one GNU time {label}, observed {len(matches)}"
        )
    if not matches[0]:
        raise ValidationError(f"{filename}: GNU time {label} is empty")
    return matches[0]


def parse_timing(raw: bytes, filename: str) -> Timing:
    try:
        text = raw.decode("utf-8")
    except UnicodeError as error:
        raise ValidationError(f"{filename}: GNU time log is not valid UTF-8") from error
    patterns = {
        "command": re.compile(
            r"^[ \t]*Command being timed:[ \t]*(?P<value>.+)$", re.MULTILINE
        ),
        "user": re.compile(
            r"^[ \t]*User time \(seconds\):[ \t]*(?P<value>\S+)[ \t]*$",
            re.MULTILINE,
        ),
        "system": re.compile(
            r"^[ \t]*System time \(seconds\):[ \t]*(?P<value>\S+)[ \t]*$",
            re.MULTILINE,
        ),
        "cpu": re.compile(
            r"^[ \t]*Percent of CPU this job got:[ \t]*(?P<value>\S+)[ \t]*$",
            re.MULTILINE,
        ),
        "elapsed": re.compile(
            r"^[ \t]*Elapsed \(wall clock\) time \([^\n]*\):"
            r"[ \t]*(?P<value>\S+)[ \t]*$",
            re.MULTILINE,
        ),
        "rss": re.compile(
            r"^[ \t]*Maximum resident set size \(kbytes\):"
            r"[ \t]*(?P<value>\S+)[ \t]*$",
            re.MULTILINE,
        ),
        "exit": re.compile(
            r"^[ \t]*Exit status:[ \t]*(?P<value>\S+)[ \t]*$", re.MULTILINE
        ),
    }
    values = {
        name: unique_time_value(text, pattern, name, filename)
        for name, pattern in patterns.items()
    }
    cpu_text = values["cpu"]
    if re.fullmatch(r"[0-9]+%", cpu_text) is None:
        raise ValidationError(f"{filename}: malformed CPU percentage {cpu_text!r}")
    user = decimal_value(values["user"], f"{filename} user time", nonnegative=True)
    system = decimal_value(values["system"], f"{filename} system time", nonnegative=True)
    cpu = decimal_value(cpu_text[:-1], f"{filename} CPU percentage", nonnegative=True)
    wall = elapsed_seconds(values["elapsed"], filename)
    if wall == 0:
        raise ValidationError(f"{filename}: elapsed wall time must be positive")
    calculated_cpu = (user + system) * 100 / wall
    if int(calculated_cpu) != int(cpu):
        raise ValidationError(
            f"{filename}: CPU percentage {cpu_text} is inconsistent with user+system "
            f"time ({canonical_decimal(calculated_cpu)}%)"
        )
    for label in ("rss", "exit"):
        if re.fullmatch(r"[0-9]+", values[label]) is None:
            raise ValidationError(f"{filename}: GNU time {label} is not an integer")
    rss = int(values["rss"])
    exit_status = int(values["exit"])
    if exit_status != 0:
        raise ValidationError(f"{filename}: solver exit status is {exit_status}, expected 0")
    return Timing(
        values["command"], user, system, cpu, values["elapsed"], wall, rss, exit_status
    )


def command_arguments(command: str, filename: str) -> list[str]:
    value = command.strip()
    # GNU time quotes the rendered full command in the archive.  Strip that
    # outer display pair before tokenizing the arguments it encloses.
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    try:
        arguments = shlex.split(value)
    except ValueError as error:
        raise ValidationError(f"{filename}: malformed command rendering") from error
    if not arguments:
        raise ValidationError(f"{filename}: empty command rendering")
    return arguments


def executable_index(arguments: Sequence[str], basename: str, filename: str) -> int:
    indices = [index for index, value in enumerate(arguments) if Path(value).name == basename]
    if len(indices) != 1:
        raise ValidationError(
            f"{filename}: expected exactly one {basename} executable, observed {len(indices)}"
        )
    return indices[0]


def command_node(arguments: Sequence[str], binder: str, filename: str) -> int:
    indices = [index for index, value in enumerate(arguments) if Path(value).name == binder]
    if len(indices) != 1:
        raise ValidationError(
            f"{filename}: expected exactly one {binder} binding helper, observed {len(indices)}"
        )
    if binder == "numabind0":
        return 0
    index = indices[0]
    if index + 1 >= len(arguments) or re.fullmatch(r"[0-9]+", arguments[index + 1]) is None:
        raise ValidationError(f"{filename}: numabind command has no numeric node")
    return int(arguments[index + 1])


def validate_binding_prefix(
    arguments: Sequence[str], executable_position: int, binder: str, filename: str
) -> int:
    prefix = arguments[:executable_position]
    node = command_node(prefix, binder, filename)
    expected_length = 1 if binder == "numabind0" else 2
    if len(prefix) != expected_length:
        raise ValidationError(
            f"{filename}: expected only the {binder} binding prefix before the solver"
        )
    if Path(prefix[0]).name != binder:
        raise ValidationError(f"{filename}: malformed {binder} binding prefix")
    if binder == "numabind" and prefix[1] != str(node):
        raise ValidationError(f"{filename}: malformed numabind node argument")
    return node


def validate_command(spec: RunSpec, timing: Timing, filename: str) -> tuple[int, Decimal]:
    arguments = command_arguments(timing.command, filename)
    if spec.group == "ld":
        index = executable_index(arguments, "LinearDesign_2D", filename)
        solver_arguments = arguments[index + 1 :]
        if len(solver_arguments) != 4:
            raise ValidationError(
                f"{filename}: dense command must have four solver arguments"
            )
        if solver_arguments[0] != "0":
            raise ValidationError(
                f"{filename}: dense command first positional argument must be 0"
            )
        if Path(solver_arguments[2]).name != HISTORICAL_CODON_TABLE:
            raise ValidationError(
                f"{filename}: dense command uses an unexpected codon-table path"
            )
        thread_text = solver_arguments[3]
        lambda_text = solver_arguments[1]
        node = validate_binding_prefix(arguments, index, "numabind0", filename)
    else:
        index = executable_index(arguments, "lineardesign-clean", filename)
        solver_arguments = arguments[index + 1 :]
        if len(solver_arguments) != 6 or (
            solver_arguments[0] != "-l"
            or solver_arguments[2] != "-c"
            or solver_arguments[4] != "-j"
        ):
            raise ValidationError(
                f"{filename}: clean command must use exactly '-l VALUE -c TABLE -j THREADS'"
            )
        if Path(solver_arguments[3]).name != HISTORICAL_CODON_TABLE:
            raise ValidationError(
                f"{filename}: clean command uses an unexpected codon-table path"
            )
        lambda_text = solver_arguments[1]
        thread_text = solver_arguments[5]
        binder = "numabind0" if spec.group == "clean" else "numabind"
        node = validate_binding_prefix(arguments, index, binder, filename)
    if re.fullmatch(r"[1-9][0-9]*", thread_text) is None:
        raise ValidationError(f"{filename}: malformed command thread count {thread_text!r}")
    command_threads = int(thread_text)
    if command_threads != spec.threads:
        raise ValidationError(
            f"{filename}: filename requests j{spec.threads}, command requests j{command_threads}"
        )
    if spec.node is not None and node != spec.node:
        raise ValidationError(
            f"{filename}: filename requests node {spec.node}, command requests node {node}"
        )
    objective_lambda = decimal_value(lambda_text, f"{filename} lambda", nonnegative=True)
    if objective_lambda != 0:
        raise ValidationError(
            f"{filename}: archived dystrophin comparison must use lambda 0, observed {lambda_text}"
        )
    return node, objective_lambda


def collect_observations(
    archive: Path, groups: Sequence[str]
) -> tuple[list[Observation], bool]:
    if not archive.is_dir():
        raise ValidationError(f"archive is not a directory: {archive}")
    manifest = parse_manifest(archive)
    pairs = discover_pairs(archive, groups, manifest)
    observations: list[Observation] = []
    for spec, output_path, time_path, _ in pairs:
        output_raw, output_hash, output_verified = read_and_hash(output_path, manifest)
        time_raw, time_hash, time_verified = read_and_hash(time_path, manifest)
        design = parse_design(output_raw, output_path.name)
        timing = parse_timing(time_raw, time_path.name)
        node, objective_lambda = validate_command(spec, timing, time_path.name)
        observations.append(
            Observation(
                spec,
                node,
                objective_lambda,
                design,
                timing,
                output_path,
                len(output_raw),
                output_hash,
                time_path,
                len(time_raw),
                time_hash,
                output_verified and time_verified,
            )
        )
    validate_cross_run_consistency(observations)
    return observations, manifest is not None


def validate_cross_run_consistency(observations: Sequence[Observation]) -> None:
    if not observations:
        raise ValidationError("no observations were selected")
    mfes = {item.design.mfe for item in observations}
    if len(mfes) != 1:
        detail = ", ".join(
            f"{item.spec.stem}={item.design.mfe_text}" for item in observations
        )
        raise ValidationError(f"selected runs report inconsistent MFE values: {detail}")
    lengths = {len(item.design.sequence) for item in observations}
    if len(lengths) != 1:
        detail = ", ".join(
            f"{item.spec.stem}={len(item.design.sequence)}" for item in observations
        )
        raise ValidationError(f"selected runs report inconsistent RNA lengths: {detail}")
    lambdas = {item.objective_lambda for item in observations}
    if lambdas != {Decimal(0)}:
        raise ValidationError("selected runs do not share the lambda-zero objective")


def observation_row(item: Observation) -> dict[str, object]:
    design = item.design
    timing = item.timing
    rna_nt = len(design.sequence)
    return {
        "group": item.spec.group,
        "run_id": item.spec.stem,
        "solver": item.spec.solver,
        "threads": item.spec.threads,
        "numa_node": item.node,
        "repeat": "" if item.spec.repeat is None else item.spec.repeat,
        "lambda": canonical_decimal(item.objective_lambda),
        "mfe_kcal": design.mfe_text,
        "cai": design.cai_text,
        "protein_aa": rna_nt // 3,
        "rna_nt": rna_nt,
        "base_pairs": design.base_pairs,
        "sequence_sha256": sha256_bytes(design.sequence.encode("ascii")),
        "structure_sha256": sha256_bytes(design.structure.encode("ascii")),
        "elapsed_text": timing.elapsed_text,
        "wall_seconds": canonical_decimal(timing.wall_seconds),
        "user_seconds": canonical_decimal(timing.user_seconds),
        "system_seconds": canonical_decimal(timing.system_seconds),
        "cpu_percent": canonical_decimal(timing.cpu_percent),
        "max_rss_kib": timing.max_rss_kib,
        "max_rss_gib": f"{Decimal(timing.max_rss_kib) / Decimal(1024 * 1024):.6f}",
        "exit_status": timing.exit_status,
        "command": timing.command,
        "output_file": item.output_path.name,
        "output_bytes": item.output_size,
        "output_sha256": item.output_sha256,
        "time_file": item.time_path.name,
        "time_bytes": item.time_size,
        "time_sha256": item.time_sha256,
        "sha256sums_verified": str(item.manifest_verified).lower(),
        "cai_interpretation": CAI_NOTE,
    }


def render_csv(observations: Sequence[Observation]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    for item in observations:
        writer.writerow(observation_row(item))
    return output.getvalue()


def markdown_decimal(value: Decimal, places: int = 2) -> str:
    return f"{value:.{places}f}"


def median_decimal(values: Sequence[Decimal]) -> Decimal:
    if not values:
        raise ValueError("median requires at least one value")
    return Decimal(str(statistics.median(values)))


def render_markdown(
    observations: Sequence[Observation], archive: Path, manifest_present: bool
) -> str:
    common = observations[0]
    lines = [
        "# Dystrophin benchmark validation",
        "",
        f"- Archive: `{archive}`",
        f"- Validated run pairs: {len(observations)}",
        f"- Protein/RNA length: {len(common.design.sequence) // 3:,} aa / "
        f"{len(common.design.sequence):,} nt",
        f"- Common MFE: {common.design.mfe_text} kcal/mol",
        "- Raw-log SHA-256 manifest: "
        + (
            "present and verified"
            if manifest_present
            else "not present; computed hashes are in CSV"
        ),
        "",
        f"_CAI interpretation: {CAI_NOTE}_",
        "",
    ]

    primary = [item for item in observations if item.spec.group in ("clean", "ld")]
    if primary:
        lines.extend(
            [
                "## Primary thread ladders",
                "",
                "| Solver | Threads | Wall | Wall seconds | CPU | Max RSS KiB | MFE | CAI |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for item in primary:
            solver = "clean sparse" if item.spec.group == "clean" else "local dense fork"
            lines.append(
                f"| {solver} | {item.spec.threads} | {item.timing.elapsed_text} | "
                f"{canonical_decimal(item.timing.wall_seconds)} | "
                f"{canonical_decimal(item.timing.cpu_percent)}% | "
                f"{item.timing.max_rss_kib:,} | {item.design.mfe_text} | "
                f"{item.design.cai_text} |"
            )
        lines.append("")

    repeats = [item for item in observations if item.spec.group == "rep"]
    if repeats:
        best = min(item.timing.wall_seconds for item in repeats)
        lines.extend(
            [
                "## Clean j16 repeats",
                "",
                "| Run | Node | Wall | Wall seconds | CPU | Max RSS KiB | Best |",
                "| --- | ---: | ---: | ---: | ---: | ---: | :---: |",
            ]
        )
        for item in repeats:
            lines.append(
                f"| `{item.spec.stem}` | {item.node} | {item.timing.elapsed_text} | "
                f"{canonical_decimal(item.timing.wall_seconds)} | "
                f"{canonical_decimal(item.timing.cpu_percent)}% | "
                f"{item.timing.max_rss_kib:,} | "
                f"{'yes' if item.timing.wall_seconds == best else ''} |"
            )
        lines.append("")

    node_comparison = [item for item in observations if item.spec.group == "nc"]
    if node_comparison:
        cells: dict[tuple[int, int], list[Observation]] = {}
        for item in node_comparison:
            cells.setdefault((item.spec.threads, item.node), []).append(item)
        lines.extend(
            [
                "## Node comparison",
                "",
                "| Threads | Node | Runs | Median wall seconds | Min wall seconds | "
                "Max wall seconds | Median RSS KiB |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for (threads, node), items in sorted(cells.items()):
            walls = [item.timing.wall_seconds for item in items]
            rss = [item.timing.max_rss_kib for item in items]
            lines.append(
                f"| {threads} | {node} | {len(items)} | "
                f"{markdown_decimal(median_decimal(walls))} | "
                f"{markdown_decimal(min(walls))} | {markdown_decimal(max(walls))} | "
                f"{statistics.median(rss):,.0f} |"
            )
        lines.extend(
            [
                "",
                "### Node-comparison observations",
                "",
                "| Threads | Node | Repeat | Wall | Wall seconds | CPU | Max RSS KiB |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for item in node_comparison:
            lines.append(
                f"| {item.spec.threads} | {item.node} | {item.spec.repeat} | "
                f"{item.timing.elapsed_text} | {canonical_decimal(item.timing.wall_seconds)} | "
                f"{canonical_decimal(item.timing.cpu_percent)}% | "
                f"{item.timing.max_rss_kib:,} |"
            )
        lines.append("")
    return "\n".join(lines)


def reject_output_collisions(
    archive: Path, output: str | None, markdown: str | None
) -> None:
    destinations: list[tuple[str, Path]] = []
    for option, value in (("--output", output), ("--markdown", markdown)):
        if value is None or value == "-":
            continue
        path = Path(value).resolve()
        try:
            path.relative_to(archive.resolve())
        except ValueError:
            pass
        else:
            raise ValidationError(f"{option} must not write inside the read-only archive")
        destinations.append((option, path))
    if len(destinations) == 2 and destinations[0][1] == destinations[1][1]:
        raise ValidationError("--output and --markdown must not name the same file")


def write_text(destination: str | None, content: str, fallback: TextIO) -> None:
    if destination is None or destination == "-":
        fallback.write(content)
        return
    try:
        Path(destination).write_text(content, encoding="utf-8")
    except OSError as error:
        raise ValidationError(f"cannot write {destination}: {error}") from error


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    reject_output_collisions(args.archive, args.output, args.markdown)
    observations, manifest_present = collect_observations(args.archive, args.groups)
    csv_text = render_csv(observations)
    markdown_text = (
        render_markdown(observations, args.archive, manifest_present)
        if args.markdown is not None
        else None
    )
    write_text(args.output, csv_text, sys.stdout)
    if markdown_text is not None:
        write_text(args.markdown, markdown_text, sys.stdout)
    return 0


def main() -> int:
    try:
        return run()
    except (OSError, UnicodeError, csv.Error, ValidationError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
