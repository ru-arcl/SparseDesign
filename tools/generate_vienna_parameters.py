#!/usr/bin/env python3
"""Regenerate src/vienna/*.h from ViennaRNA's own Turner-2004 parameter set.

WHY THIS EXISTS. The parameter tables under src/vienna/ are numeric constants
published with the Turner 2004 nearest-neighbour model and distributed by
ViennaRNA. They were previously carried as a header transcribed by a third party,
which left the tree asserting it contained no code from a restricted distribution
while shipping bytes identical to that distribution's copy. This tool removes the
question by deriving the same numbers directly from ViennaRNA, so the provenance
is a hash-pinned upstream parameter set plus the transformation below.

SOURCE. `RNA.parameter_set_rna_turner2004` from the ViennaRNA Python package,
which is the standard `RNAfold parameter file v2.0` text. Its complete SHA-256
digest must equal TURNER2004_PARAMETER_SHA256 in this self-contained script.

TRANSFORMATION. Two mechanical steps, both required to match the C layout:

  1. Index offset. ViennaRNA's file lists 7 base-pair types (CG GC GU UG AU UA
     NN). The C arrays are dimensioned [NBPAIRS+1] with index 0 reserved for "no
     pair", filled with VIE_INF. Nucleotide axes are already 5 wide in the source
     with index 0 = N, except in int22.

  2. int22 unknown-base fill. ViennaRNA stores int22 only for the 6 canonical
     pairs and 4 canonical bases (6*6*4^4 = 9216 values), and fills the N and NN
     entries with the MAXIMUM over the corresponding known index -- a
     conservative choice, since an unknown base must not look favourable. This
     tool reproduces that cascade, innermost axis first.

VERIFICATION. `--check` regenerates in memory and compares complete committed
headers, ignoring only the ViennaRNA version in the provenance comment. A changed
upstream parameter digest, numeric value, declaration, or transformation fails.
The comparison does not modify any files and needs no C++ compiler.

Usage:
    python3 tools/generate_vienna_parameters.py --vienna-python PATH [--write | --check]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import textwrap

HERE = Path(__file__).resolve().parent
VIENNA_DIR = HERE.parent / "src" / "vienna"

INF_NAME = "VIE_INF"
INF = 10000000
NBPAIRS = 7
TURNER2004_PARAMETER_SHA256 = "2a43345a495850cfd2e0a78c57c6e02085e6df3c53496fe3289dc294d21732ad"

PROBE = r"""
import json, RNA
print(json.dumps({"version": RNA.__version__,
                  "text": RNA.parameter_set_rna_turner2004}))
"""

BANNER = """\
/*
 * ViennaRNA Turner 2004 nearest-neighbour parameters.
 *
 * GENERATED FILE -- do not edit by hand.
 *   generator : tools/generate_vienna_parameters.py
 *   source    : ViennaRNA {version}, parameter_set_rna_turner2004
 *   digest    : sha256 {digest}
 *
 * These are published thermodynamic constants transcribed from the upstream
 * ViennaRNA parameter set. Index 0 of each base-pair axis is the "no pair"
 * slot and holds {inf_name}; int22's unknown-base entries carry ViennaRNA's
 * maximum-fill. Regenerate with --write and verify with --check.
 */
"""


def section_map(text: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    current: str | None = None
    buffer: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            if current is not None:
                out[current] = buffer
            current = line.strip("# ").strip()
            buffer = []
        elif current is not None:
            buffer.append(line)
    if current is not None:
        out[current] = buffer
    return out


def values(sections: dict[str, list[str]], name: str) -> list[int]:
    out: list[int] = []
    for line in sections.get(name, []):
        line = re.sub(r"/\*.*?\*/", "", line)
        for token in line.split():
            if re.fullmatch(r"-?\d+", token):
                out.append(int(token))
            elif token.upper() == "INF":
                out.append(INF)
    return out


def loop_table(sections: dict[str, list[str]], name: str) -> tuple[list[str], list[int]]:
    sequences: list[str] = []
    energies: list[int] = []
    for line in sections.get(name, []):
        parts = line.split()
        if len(parts) >= 2 and re.fullmatch(r"[ACGU]+", parts[0]):
            sequences.append(parts[0])
            energies.append(int(parts[1]))
    return sequences, energies


def pair_indexed(sections: dict[str, list[str]], name: str, stride: int) -> list[int]:
    """Prepend the no-pair row: ViennaRNA rows are pairs 1..7."""
    return [INF] * stride + values(sections, name)


def pair_pair_indexed(
    sections: dict[str, list[str]], name: str, stride: int
) -> list[int]:
    source = values(sections, name)
    out: list[int] = []
    for first in range(NBPAIRS + 1):
        for second in range(NBPAIRS + 1):
            if first == 0 or second == 0:
                out += [INF] * stride
            else:
                start = ((first - 1) * NBPAIRS + (second - 1)) * stride
                out += source[start : start + stride]
    return out


def int22(sections: dict[str, list[str]]) -> list[int]:
    source = values(sections, "int22")
    table = [
        [[[[[INF] * 5 for _ in range(5)] for _ in range(5)] for _ in range(5)]
         for _ in range(NBPAIRS + 1)]
        for _ in range(NBPAIRS + 1)
    ]
    cursor = 0
    canonical = range(1, 7)
    for first in canonical:
        for second in canonical:
            for a in range(1, 5):
                for b in range(1, 5):
                    for c in range(1, 5):
                        for d in range(1, 5):
                            table[first][second][a][b][c][d] = source[cursor]
                            cursor += 1
    if cursor != len(source):
        raise ValueError(f"int22 consumed {cursor} of {len(source)} values")

    # Unknown-base and unknown-pair fill, innermost axis first so each cascade
    # sees the zeros written by the previous one.
    for first in canonical:
        for second in canonical:
            block = table[first][second]
            for a in range(1, 5):
                for b in range(1, 5):
                    for c in range(1, 5):
                        block[a][b][c][0] = max(block[a][b][c][d] for d in range(1, 5))
            for a in range(1, 5):
                for b in range(1, 5):
                    for d in range(5):
                        block[a][b][0][d] = max(block[a][b][c][d] for c in range(1, 5))
            for a in range(1, 5):
                for c in range(5):
                    for d in range(5):
                        block[a][0][c][d] = max(block[a][b][c][d] for b in range(1, 5))
            for b in range(5):
                for c in range(5):
                    for d in range(5):
                        block[0][b][c][d] = max(block[a][b][c][d] for a in range(1, 5))
    for first in canonical:
        for a in range(5):
            for b in range(5):
                for c in range(5):
                    for d in range(5):
                        table[first][7][a][b][c][d] = max(
                            table[first][s][a][b][c][d] for s in canonical)
    for second in range(1, 8):
        for a in range(5):
            for b in range(5):
                for c in range(5):
                    for d in range(5):
                        table[7][second][a][b][c][d] = max(
                            table[f][second][a][b][c][d] for f in canonical)
    for a in range(5):
        for b in range(5):
            for c in range(5):
                for d in range(5):
                    table[7][7][a][b][c][d] = max(
                        table[7][s][a][b][c][d] for s in canonical)

    flat: list[int] = []
    for first in range(NBPAIRS + 1):
        for second in range(NBPAIRS + 1):
            for a in range(5):
                for b in range(5):
                    for c in range(5):
                        for d in range(5):
                            flat.append(table[first][second][a][b][c][d])
    return flat


def build_all(text: str) -> dict[str, object]:
    s = section_map(text)
    arrays: dict[str, object] = {}
    stack = values(s, "stack")
    flat = [INF] * (NBPAIRS + 1)
    for pair in range(NBPAIRS):
        flat += [INF] + stack[pair * NBPAIRS : (pair + 1) * NBPAIRS]
    arrays["stack37"] = flat
    for target, source in (
        ("mismatchH37", "mismatch_hairpin"),
        ("mismatchI37", "mismatch_internal"),
        ("mismatchM37", "mismatch_multi"),
        ("mismatch1nI37", "mismatch_internal_1n"),
        ("mismatch23I37", "mismatch_internal_23"),
        ("mismatchExt37", "mismatch_exterior"),
    ):
        arrays[target] = pair_indexed(s, source, 25)
    arrays["dangle5_37"] = pair_indexed(s, "dangle5", 5)
    arrays["dangle3_37"] = pair_indexed(s, "dangle3", 5)
    arrays["hairpin37"] = values(s, "hairpin")
    arrays["bulge37"] = values(s, "bulge")
    arrays["internal_loop37"] = values(s, "internal")
    arrays["int11_37"] = pair_pair_indexed(s, "int11", 25)
    arrays["int21_37"] = pair_pair_indexed(s, "int21", 125)
    arrays["int22_37"] = int22(s)
    ml = values(s, "ML_params"); ninio = values(s, "NINIO"); misc = values(s, "Misc")
    arrays["ML_BASE37"] = ml[0]
    arrays["ML_closing37"] = ml[2]
    arrays["ML_intern37"] = ml[4]
    arrays["ninio37"] = ninio[0]
    arrays["MAX_NINIO"] = ninio[2]
    arrays["TerminalAU37"] = misc[2]
    for target, source in (
        ("Triloop", "Triloops"),
        ("Tetraloop", "Tetraloops"),
        ("Hexaloop", "Hexaloops"),
    ):
        sequences, energies = loop_table(s, source)
        arrays[target + "s"] = sequences
        arrays[target + "37"] = energies
    return arrays




# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------

def fmt(value: int) -> str:
    return INF_NAME if value == INF else str(value)


def wrap_flat(name: str, decl: str, data: list[int], per_line: int = 12) -> str:
    body = []
    for start in range(0, len(data), per_line):
        body.append("    " + ", ".join(fmt(v) for v in data[start : start + per_line]))
    return f"{decl} = {{\n" + ",\n".join(body) + "\n};\n"


def nest(data: list[int], dims: list[int], depth: int = 0) -> str:
    """Render a flat list as a braced C initialiser of the given shape."""
    pad = "    " * (depth + 1)
    if len(dims) == 1:
        return "{" + ", ".join(fmt(v) for v in data) + "}"
    stride = 1
    for d in dims[1:]:
        stride *= d
    parts = [nest(data[i * stride : (i + 1) * stride], dims[1:], depth + 1)
             for i in range(dims[0])]
    joiner = ",\n" + pad
    return "{\n" + pad + joiner.join(parts) + "\n" + "    " * depth + "}"


def c_string_array(name: str, size: int, sequences: list[str]) -> str:
    joined = " ".join(sequences) + " "
    return f'char {name}[{size}] =\n    "{joined}";\n'


def emit(arrays: dict, version: str, digest: str) -> dict[str, str]:
    banner = BANNER.format(version=version, digest=digest, inf_name=INF_NAME)
    P1 = NBPAIRS + 1

    head = [banner]
    head.append(f"#ifndef {INF_NAME}\n#define {INF_NAME} {INF}\n#endif\n")
    head.append(f"#ifndef NBPAIRS\n#define NBPAIRS {NBPAIRS}  // NP CG GC GU UG AU UA NN\n#endif\n")
    head.append("#ifndef SPECIAL_HP\n#define SPECIAL_HP\n#endif\n")
    head.append(
        "// Scalar model parameters. All but lxc37 come from the upstream\n"
        "// ML_params, NINIO and Misc sections; lxc37 is ViennaRNA's built-in\n"
        "// logarithmic loop-extrapolation constant and is stated, not derived.\n"
        "double lxc37 = 107.856;\n"
        f"int ML_intern37 = {arrays['ML_intern37']};\n"
        f"int ML_closing37 = {arrays['ML_closing37']};\n"
        f"int ML_BASE37 = {arrays['ML_BASE37']};\n"
        f"int MAX_NINIO = {arrays['MAX_NINIO']};\n"
        f"int ninio37 = {arrays['ninio37']};\n"
        f"int TerminalAU37 = {arrays['TerminalAU37']};\n"
    )
    for stem, size in (("Triloop", 241), ("Tetraloop", 281), ("Hexaloop", 361)):
        head.append(c_string_array(stem + "s", size, arrays[stem + "s"]))
        vals = arrays[stem + "37"]
        head.append(f"int {stem}37[{len(vals)}] = {{" + ", ".join(str(v) for v in vals) + "};\n")
    head.append(f"int stack37[NBPAIRS+1][NBPAIRS+1] =\n    {nest(arrays['stack37'], [P1, P1])};\n")
    for name in ("hairpin37", "bulge37", "internal_loop37"):
        head.append(f"int {name}[31] = {{" + ", ".join(fmt(v) for v in arrays[name]) + "};\n")
    for name in ("mismatchI37", "mismatchH37", "mismatchM37",
                 "mismatch1nI37", "mismatch23I37", "mismatchExt37"):
        head.append(f"int {name}[NBPAIRS+1][5][5] =\n    {nest(arrays[name], [P1, 5, 5])};\n")
    for name in ("dangle5_37", "dangle3_37"):
        head.append(f"int {name}[NBPAIRS+1][5] =\n    {nest(arrays[name], [P1, 5])};\n")

    guard = (f"#ifndef {INF_NAME}\n#define {INF_NAME} {INF}\n#endif\n"
             f"#ifndef NBPAIRS\n#define NBPAIRS {NBPAIRS}\n#endif\n")
    return {
        "energy_parameter.h": "\n".join(head),
        "intl11.h": banner + "\n" + guard + "\n"
            + f"int int11_37[NBPAIRS+1][NBPAIRS+1][5][5] =\n    {nest(arrays['int11_37'], [P1, P1, 5, 5])};\n",
        "intl21.h": banner + "\n" + guard + "\n"
            + f"int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5] =\n    {nest(arrays['int21_37'], [P1, P1, 5, 5, 5])};\n",
        "intl22.h": banner + "\n" + guard + "\n"
            + f"int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5] =\n    {nest(arrays['int22_37'], [P1, P1, 5, 5, 5, 5])};\n",
    }


def probe(interpreter: str) -> tuple[str, str, str]:
    import hashlib
    import json

    completed = subprocess.run(
        [interpreter, "-c", PROBE], capture_output=True, timeout=300
    )
    if completed.returncode != 0:
        raise SystemExit(
            "error: could not read parameters from "
            f"{interpreter}: {completed.stderr.decode(errors='replace')[:300]}"
        )
    payload = json.loads(completed.stdout.decode())
    digest = hashlib.sha256(payload["text"].encode("utf-8")).hexdigest()
    return payload["version"], payload["text"], digest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vienna-python", default=os.environ.get("SPARSEDESIGN_VIENNA_PYTHON"))
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=VIENNA_DIR,
                        help="header directory to write/check (default: src/vienna)")
    arguments = parser.parse_args(argv)
    if not arguments.vienna_python:
        raise SystemExit(
            "error: pass --vienna-python PATH or set SPARSEDESIGN_VIENNA_PYTHON"
        )
    version, text, digest = probe(arguments.vienna_python)
    if digest != TURNER2004_PARAMETER_SHA256:
        raise SystemExit(f"error: Turner-2004 source digest mismatch: expected "
                         f"{TURNER2004_PARAMETER_SHA256}, observed {digest}")
    arrays = build_all(text)
    print(f"ViennaRNA {version}, parameter digest {digest[:16]}...")
    files = emit(arrays, version, digest)
    if arguments.check:
        mismatches = []
        for name, content in files.items():
            path = arguments.output_dir / name
            if not path.is_file():
                mismatches.append(name + " (missing)")
                continue
            # Source release versions may differ while the hash-pinned text is identical.
            normalize = lambda value: re.sub(r"(?m)^ \*   source    : ViennaRNA [^,]+,",
                                            " *   source    : ViennaRNA VERSION,", value)
            if normalize(path.read_text(encoding="ascii")) != normalize(content):
                mismatches.append(name)
            else:
                print(f"  verified {name}")
        if mismatches:
            raise SystemExit("error: generated parameter headers differ: " + ", ".join(mismatches))
        print("parameter verification: PASS (pinned source and all four headers)")
    elif arguments.write:
        arguments.output_dir.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            (arguments.output_dir / name).write_text(content, encoding="ascii")
            print(f"  wrote src/vienna/{name} ({len(content)} bytes)")
    else:
        for name, content in files.items():
            print(f"  would write src/vienna/{name} ({len(content)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
