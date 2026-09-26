#!/usr/bin/env python3
"""Reproducible release validation, external refolding, constraints and sensitivity.

Python 3.10+, ViennaRNA==2.7.2; compiler needed for --build. The solver traversal
and dense modes share the energy implementation. Only ViennaRNA computations
below are an independent folding implementation. Timings are instrumented
validation costs, not controlled performance comparisons.
"""
import argparse
import csv
import fcntl
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import time
import zipfile


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def load_table(path):
    rows = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) != 3 or row[0].startswith("#") or row[0] == "codon":
                continue
            rows.append((row[0].replace("T", "U"), row[1], float(row[2])))
    maxima = {aa: max(f for _, a, f in rows if a == aa) for _, aa, _ in rows}
    weights = {c: f / maxima[a] for c, a, f in rows}
    translate = {c: a for c, a, _ in rows}
    options = {aa: sorted(c for c, a, _ in rows if a == aa) for aa in maxima}
    return weights, translate, options


def prepare_inputs(root, out, ward):
    out.mkdir(parents=True, exist_ok=True)
    (out / "inputs").mkdir(exist_ok=True)
    target = out / "inputs" / "human-rounded-original.csv"
    shutil.copyfile(root / "lineardesign-clean/data/codon_usage_freq_table_human.csv", target)
    inputs = {"schema_version": 1, "selection": "Prespecified deterministic manageable validation panel, not a representative performance sample.", "cases": [], "tables": {"human": "inputs/human-rounded-original.csv"}}
    case = inputs["cases"]
    # Natural inputs are the lexicographically first accession in each selected
    # length band of the archived, length-balanced 400-protein panel.
    manifest_path = root / "research/data/uniprot-sprot-2026_02-sparse-balanced-400.manifest.json"
    archive_path = root / "research/data/uniprot-sprot-2026_02-sparse-balanced-400.zip"
    manifest = json.loads(manifest_path.read_text())
    inputs["natural_source"] = {"manifest_sha256": digest(manifest_path), "archive_sha256": digest(archive_path)}
    natural = []
    with zipfile.ZipFile(archive_path) as archive:
        for lo, hi in ((10, 24), (25, 49), (50, 74), (75, 99), (100, 120)):
            selected = min((r for r in manifest["selected"] if lo <= r["aa_length"] <= hi), key=lambda r: r["accession"])
            raw = archive.read(selected["fasta_path"])
            assert hashlib.sha256(raw).hexdigest() == selected["fasta_sha256"]
            protein = "".join(x for x in raw.decode().splitlines() if not x.startswith(">"))
            natural.append((selected["accession"], protein))
    for accession, protein in natural:
        for lam in (0, 0.25, 0.5, 1, 2, 4, 8):
            case.append(dict(id=f"natural-{accession}-l{lam}", group="natural_lambda", protein=protein, lambda_user=lam, table="human", motifs=[], ensemble=True))
    # Complete published failing-case panel, retaining its documented lambda.
    aa_map = dict(zip("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val End".split(), "ARNDCQEGHILKMFPSTWYV*"))
    if ward is not None:
        (out / "inputs/ward").mkdir(exist_ok=True)
        for name in ("lineardesign_beaking_cases.csv", "derna_breaking_cases.csv", "homosapiens.txt"):
            shutil.copyfile(ward / "data" / name, out / "inputs/ward" / name)
        shutil.copyfile(ward / "LICENSE", out / "inputs/ward/LICENSE")
        inputs["ward_source"] = {"url": "https://github.com/maxhwardg/mrna_folding_comparison", "commit": subprocess.check_output(["git", "-C", str(ward), "rev-parse", "HEAD"], text=True).strip(), "license": "MIT; inputs/ward/LICENSE", "paper": "https://doi.org/10.1093/bib/bbaf386"}
        rows = []
        for line in (ward / "data/homosapiens.txt").read_text().splitlines():
            bits = line.split()
            if len(bits) >= 3:
                rows.append((bits[1].replace("T", "U"), aa_map[bits[0]], float(bits[2])))
        totals = {a: sum(f for _, aa, f in rows if aa == a) for _, a, _ in rows}
        for precision in ("full", "rounded2"):
            name = f"ward-human-{precision}"
            rel = f"inputs/{name}.csv"
            with (out / rel).open("w") as handle:
                handle.write("# Derived from Ward homosapiens.txt counts; normalized within amino acid.\n")
                for c, a, f in rows:
                    value = f / totals[a]
                    if precision == "rounded2":
                        value = round(value, 2)
                    handle.write(f"{c},{a},{value:.17g}\n")
            inputs["tables"][name] = rel
        for filename in ("lineardesign_beaking_cases.csv", "derna_breaking_cases.csv"):
            with (ward / "data" / filename).open() as handle:
                for i, row in enumerate(csv.DictReader(handle, skipinitialspace=True)):
                    case.append(dict(id=f"ward-{filename.split('_')[0]}-{i+1}", group="ward", protein=row["protein"], lambda_user=float(row["lambda"]), table="ward-human-full", motifs=[], ensemble=len(row["protein"]) <= 120, historical_failure=row["notes"]))
        for accession, protein in natural:
            for lam in (0, 1, 4):
                for precision in ("full", "rounded2"):
                    case.append(dict(id=f"precision-{accession}-l{lam}-{precision}", group="codon_precision", protein=protein, lambda_user=lam, table=f"ward-human-{precision}", motifs=[], ensemble=True))
    # Longer previously high-density families, including the two lambda endpoints.
    for family in ("S", "YI"):
        for length in (32, 64, 128, 256):
            protein = (family * length)[:length]
            for lam in (0, 1, 4):
                case.append(dict(id=f"family-{family}-{length}-l{lam}", group="synthetic_family", protein=protein, lambda_user=lam, table="human", motifs=[], ensemble=length <= 64))
    # Pattern inventories are nested; overlapping UUUUU/UUUUUU and cross-codon
    # restriction-site spellings exercise DFA state history, not codon filtering.
    motif_sets = [[], ["UUUUU"], ["UUUUU", "UUUUUU", "GGCGCC", "GCGGCCGC"], ["UUUUU", "UUUUUU", "GGCGCC", "GCGGCCGC", "GAAUUC", "GGUUCC", "CUCGAG", "AAGCUU"], ["UUUUU", "UUUUUU", "GGCGCC", "GCGGCCGC", "GAAUUC", "GGUUCC", "CUCGAG", "AAGCUU", "GGAUCC", "CCCGGG", "GUCGAC", "ACCGGU", "AGAUCU", "UCUAGA", "GGGCCC", "CCUAGG"]]
    for accession, protein in (natural[2], natural[4]):
        for motifs in motif_sets:
            for lam in (0, 1, 4):
                case.append(dict(id=f"motif-{accession}-m{len(motifs)}-l{lam}", group="constraints", protein=protein, lambda_user=lam, table="human", motifs=motifs, ensemble=False))
    for protein, motifs in (("M" * 50, ["AUG"]), ("K" * 60, ["AAA", "AAG"]), ("F" * 80, ["UUUU"])):
        case.append(dict(id=f"constraint-{protein[0]}-infeasible" if protein[0] != "F" else "constraint-F-overlap-feasible", group="constraint_edge", protein=protein, lambda_user=1, table="human", motifs=motifs, expected_infeasible=protein[0] != "F", ensemble=False))
    inputs["table_sha256"] = {key: digest(out / path) for key, path in inputs["tables"].items()}
    dump(out / "inputs.json", inputs)
    return inputs


def build(source, out):
    sources = [source / "src" / x for x in ("fold_turner.cc", "energy.cc", "dfa.cc", "codon_table.cc")]
    driver = Path(__file__).with_suffix(".cc").resolve()
    command = ["g++", "-std=c++11", "-O3", "-DNDEBUG", "-DLDCLEAN_SQUARE_MEMOS", "-fopenmp", "-I" + str(source / "src"), str(driver), *map(str, sources), "-o", str(out / "validation-driver")]
    result = subprocess.run(command, text=True, capture_output=True)
    (out / "build.log").write_text(result.stdout + result.stderr)
    result.check_returncode()
    import RNA
    dump(out / "build.json", dict(command=command, compiler=subprocess.check_output(["g++", "--version"], text=True), python=sys.version, viennarna=RNA.__version__, platform=platform.platform(), source_sha256={str(p.relative_to(source)): digest(p) for p in sorted((source / "src").rglob("*")) if p.is_file()}, driver_sha256=digest(driver), script_sha256=digest(__file__), binary_sha256=digest(out / "validation-driver"), parameter_string_sha256=hashlib.sha256(RNA.parameter_set_rna_turner2004.encode()).hexdigest(), model={"temperature_celsius": 37, "dangles": 0, "special_hp": True, "noLP": False, "min_loop_size": 3, "internal_bulge_unpaired_cap": int(RNA.MAXLOOP), "salt_molar": 1.021, "energy_parameter_set": "Turner2004"}))


def run_driver(binary, mode, table=None, sequence=None, lam=0, motifs=(), timeout=300):
    command = [str(binary), mode]
    if mode == "eval":
        command += [sequence, motifs]
    elif mode != "special":
        command += [str(table), sequence, str(lam), *motifs]
    start = time.monotonic()
    try:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
        value = json.loads(completed.stdout)
        if isinstance(value, dict):
            value.update(process_seconds=time.monotonic() - start, returncode=completed.returncode, stderr=completed.stderr)
        return value
    except subprocess.TimeoutExpired:
        return dict(status="timeout", process_seconds=time.monotonic() - start, timeout_seconds=timeout)


def fold_compound(rna, dangles=0):
    import RNA
    md = RNA.md()
    md.temperature = 37
    md.dangles = dangles
    md.special_hp = 1
    md.noLP = 0
    md.noGU = 0
    md.noGUclosure = 0
    md.min_loop_size = 3
    md.gquad = 0
    return RNA.fold_compound(rna, md)


def structure_metrics(rna, structure):
    stack, pairs = [], {}
    if len(rna) != len(structure):
        raise ValueError("length mismatch")
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if not stack:
                raise ValueError("unmatched close")
            j = stack.pop()
            if rna[j] + rna[i] not in {"AU", "UA", "GC", "CG", "GU", "UG"}:
                raise ValueError("noncanonical base pair")
            pairs[j] = i
        elif c != ".":
            raise ValueError("bad dot bracket character")
    if stack:
        raise ValueError("unmatched open")
    max_two_loop, smallest_hairpin = 0, None
    for i, j in pairs.items():
        children = []
        k = i + 1
        while k < j:
            if k in pairs:
                children.append((k, pairs[k]))
                k = pairs[k] + 1
            else:
                k += 1
        if len(children) == 1:
            p, q = children[0]
            max_two_loop = max(max_two_loop, p - i - 1 + j - q - 1)
        elif not children:
            smallest_hairpin = min(smallest_hairpin or (j - i - 1), j - i - 1)
    return dict(base_pairs=len(pairs), max_two_loop_unpaired=max_two_loop, min_hairpin_unpaired=smallest_hairpin)


def verify(case, value, table, ensemble=True):
    if value["status"] != "ok":
        value["checks_pass"] = bool(case.get("expected_infeasible") and value["status"] == "error" and "eliminate every synonymous" in value.get("message", ""))
        return value
    weights, translation, _ = load_table(table)
    rna, structure = value["rna"], value["structure"]
    metrics = structure_metrics(rna, structure)
    protein = "".join(translation[rna[i:i + 3]] for i in range(0, len(rna), 3))
    penalty = -sum(math.log(weights[rna[i:i + 3]]) for i in range(0, len(rna), 3))
    cai = math.exp(-penalty / (len(rna) // 3))
    recomputed = value["energy_centikcal"] + 100 * case["lambda_user"] * penalty
    tolerance = 1e-6 + 1e-10 * abs(value["dp_centikcal"])
    checks = dict(translation=protein == case["protein"], motifs=not any(m.replace("T", "U") in rna for m in case["motifs"]), cai=abs(cai - value["cai"]) <= 1e-12, codon_penalty=abs(penalty - value["codon_penalty"]) <= 1e-10 * (1 + penalty), dp_objective=abs(recomputed - value["dp_centikcal"]) <= tolerance, admissible=metrics["max_two_loop_unpaired"] <= 30 and (metrics["min_hairpin_unpaired"] is None or metrics["min_hairpin_unpaired"] >= 3))
    for key in ("dense_reference_centikcal", "dense_wavefront_centikcal", "dense_right_normal_centikcal"):
        if key in value:
            checks[key] = abs(value[key] - value["dp_centikcal"]) <= tolerance
    fc = fold_compound(rna)
    external_eval = float(fc.eval_structure(structure))
    external_structure, external_mfe = fc.mfe()
    checks["external_structure_energy"] = abs(external_eval - value["energy_centikcal"] / 100) <= 0.001
    checks["external_fixed_sequence_mfe"] = abs(external_mfe - value["energy_centikcal"] / 100) <= 0.001
    value.update(structure_metrics=metrics, checks=checks, checks_pass=all(checks.values()), objective_error_centikcal=recomputed - value["dp_centikcal"], external_d0=dict(eval_structure_kcal=external_eval, mfe_kcal=float(external_mfe), mfe_structure=external_structure))
    if ensemble:
        value["ensemble"] = {}
        for dangles in (0, 2):
            efc = fc if dangles == 0 else fold_compound(rna, dangles)
            _, mfe = efc.mfe()
            efc.exp_params_rescale(mfe)
            _, free_energy = efc.pf()
            probabilities = efc.bpp()
            pair_sum = sum(probabilities[i][j] for i in range(1, len(rna) + 1) for j in range(i + 1, len(rna) + 1))
            aup = 1 - 2 * pair_sum / len(rna)
            value["ensemble"][f"d{dangles}"] = dict(mfe_kcal=float(mfe), free_energy_kcal=float(free_energy), average_unpaired_probability=aup, emitted_structure_kcal=float(efc.eval_structure(structure)))
            if not (-1e-9 <= aup <= 1 + 1e-9) or free_energy > mfe + 1e-3:
                value["checks_pass"] = False
    return value


def exhaustive(binary, table, out):
    """Independent synonymous enumeration + ViennaRNA folding, including breakpoints."""
    weights, _, options = load_table(table)
    records, comparisons = [], []
    for protein, motifs in (("MFL", []), ("RSL", []), ("SSSS", []), ("MWFMWF", []), ("MSWFWFM", []), ("SSSS", ["TAG", "UCG", "CAG"]), ("SLS", ["UCUC", "CUCU", "UCAG"])):
        possibilities = []
        for codons in itertools.product(*(options[aa] for aa in protein)):
            rna = "".join(codons)
            if any(m.replace("T", "U") in rna for m in motifs):
                continue
            _, mfe = fold_compound(rna).mfe()
            penalty = -sum(math.log(weights[c]) for c in codons)
            # Recover the exact 0.01-kcal integer from ViennaRNA's float return.
            e = round(mfe * 100) / 100
            possibilities.append((e, penalty))
            records.append(dict(protein=protein, motifs=motifs, rna=rna, mfe_kcal=e, codon_penalty=penalty))
        unique = sorted(set(possibilities))
        # Dominated penalty/energy points cannot support an optimum. Keep only
        # the best energy at each penalty before enumerating envelope crossings.
        by_penalty = {}
        for e, p in unique:
            by_penalty[p] = min(e, by_penalty.get(p, math.inf))
        points = [(e, p) for p, e in by_penalty.items()]
        support = set()
        for (e1, p1), (e2, p2) in itertools.combinations(points, 2):
            if abs(p1 - p2) < 1e-12:
                continue
            lam = (e2 - e1) / (p1 - p2)
            if 1e-8 < lam <= 100 and e1 + lam * p1 <= min(e + lam * p for e, p in points) + 1e-8:
                support.add(round(lam, 12))
        lambdas = {0., 1., 4.}
        for lam in sorted(support):
            lambdas.update((lam - 1e-8, lam, lam + 1e-8))
        for lam in sorted(lambdas):
            value = run_driver(binary, "design", table, protein, lam, motifs)
            case = dict(protein=protein, motifs=motifs, lambda_user=lam)
            value = verify(case, value, table, ensemble=False)
            optimum = min(e + lam * p for e, p in possibilities) * 100
            value.update(case=case, feasible_cds=len(possibilities), independent_brute_centikcal=optimum, independent_brute_error_centikcal=value["dp_centikcal"] - optimum, breakpoint=any(abs(lam - b) < 1e-7 for b in support))
            value["checks_pass"] &= abs(value["dp_centikcal"] - optimum) < 1e-5
            comparisons.append(value)
    with (out / "exhaustive-sequences.jsonl").open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    dump(out / "exhaustive.json", comparisons)
    return comparisons


def loop_validation(binary, table, out):
    import RNA
    cases = []
    for special in run_driver(binary, "special"):
        rna = special["sequence"]
        cases.append(dict(name=f"special-{rna}", rna=rna, structure="(" + "." * (len(rna) - 2) + ")", kind="special_hairpin"))
    for left, right in ((0, 29), (0, 30), (0, 31), (15, 14), (15, 15), (15, 16)):
        rna = "G" * 6 + "A" * left + "C" * 6 + "AAA" + "G" * 6 + "A" * right + "C" * 6
        structure = "(" * 6 + "." * left + "(" * 6 + "..." + ")" * 6 + "." * right + ")" * 6
        cases.append(dict(name=f"twoloop-{left}-{right}", rna=rna, structure=structure, kind="cap_boundary", total_unpaired=left + right))
    results = []
    for case in cases:
        value = run_driver(binary, "eval", sequence=case["rna"], motifs=case["structure"])
        fc = fold_compound(case["rna"])
        value.update(case=case, external_eval_kcal=float(fc.eval_structure(case["structure"])))
        value["checks_pass"] = abs(value["energy_centikcal"] / 100 - value["external_eval_kcal"]) <= .001
        fixed = run_driver(binary, "fixed", table, case["rna"])
        _, external_mfe = fc.mfe()
        value.update(fixed_sparse_centikcal=fixed["dp_centikcal"], external_unconstrained_mfe_kcal=float(external_mfe))
        value["checks_pass"] &= abs(fixed["dp_centikcal"] / 100 - external_mfe) <= .001
        if case["kind"] == "cap_boundary":
            # Forced 31-unpaired loops can be evaluated, but are outside the
            # folding recurrence. This explicitly tests the external boundary;
            # our cap is checked by emitted-structure traversal and matched MFE.
            fc.hc_add_from_db(case["structure"].replace(".", "x"), RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
            _, forced = fc.mfe()
            value["external_forced_mfe_kcal"] = float(forced)
            value["external_forced_feasible"] = forced < 99999
            value["checks_pass"] &= (forced < 99999) == (case["total_unpaired"] <= 30)
        results.append(value)
    dump(out / "loops.json", results)
    return results


def summary(out):
    records = [json.loads(x) for x in (out / "designs.jsonl").read_text().splitlines()]
    exhaustive_records = json.loads((out / "exhaustive.json").read_text())
    loops = json.loads((out / "loops.json").read_text())
    groups = {}
    for group in sorted({r["case"]["group"] for r in records}):
        selected = [r for r in records if r["case"]["group"] == group]
        groups[group] = dict(cases=len(selected), passed=sum(r.get("checks_pass", False) for r in selected), statuses={s: sum(r["status"] == s for r in selected) for s in sorted({r["status"] for r in selected})})
    ok = [r for r in records if r["status"] == "ok"]
    precision = [r for r in ok if r["case"]["group"] == "codon_precision"]
    precision_pairs = []
    for r in precision:
        if r["case"]["table"].endswith("full"):
            other = next(s for s in precision if s["case"]["protein"] == r["case"]["protein"] and s["case"]["lambda_user"] == r["case"]["lambda_user"] and s["case"]["table"].endswith("rounded2"))
            precision_pairs.append(dict(protein=r["case"]["id"].split("-")[1], lambda_user=r["case"]["lambda_user"], same_sequence=r["rna"] == other["rna"], mfe_difference_kcal=other["mfe_kcal"] - r["mfe_kcal"], candidate_ratio=other["candidates"] / r["candidates"] if r["candidates"] else None))
    natural = [r for r in ok if r["case"]["group"] == "natural_lambda"]
    inversion = {"d0_vs_d2_mfe": 0, "d0_mfe_vs_d0_ensemble": 0, "d0_mfe_vs_d0_aup": 0, "strict_mfe_pairs": 0}
    for a, b in itertools.combinations(natural, 2):
        if a["case"]["protein"] != b["case"]["protein"] or a["rna"] == b["rna"]:
            continue
        diff = a["mfe_kcal"] - b["mfe_kcal"]
        if abs(diff) <= .001:
            continue
        inversion["strict_mfe_pairs"] += 1
        for key, dangle, measure in (("d0_vs_d2_mfe", "d2", "mfe_kcal"), ("d0_mfe_vs_d0_ensemble", "d0", "free_energy_kcal"), ("d0_mfe_vs_d0_aup", "d0", "average_unpaired_probability")):
            if diff * (a["ensemble"][dangle][measure] - b["ensemble"][dangle][measure]) < -1e-8:
                inversion[key] += 1
    result = dict(schema_version=1, groups=groups, total_design_cases=len(records), passed_design_cases=sum(r.get("checks_pass", False) for r in records), failures=[r["case"]["id"] for r in records if not r.get("checks_pass")], exhaustive=dict(comparisons=len(exhaustive_records), passed=sum(r["checks_pass"] for r in exhaustive_records), feasible_sequence_records=sum(1 for _ in (out / "exhaustive-sequences.jsonl").open()), breakpoint_comparisons=sum(r["breakpoint"] for r in exhaustive_records)), loops=dict(cases=len(loops), passed=sum(r["checks_pass"] for r in loops)), max_objective_error_centikcal=max(abs(r["objective_error_centikcal"]) for r in ok), max_external_energy_error_kcal=max(abs(r["mfe_kcal"] - r["external_d0"]["eval_structure_kcal"]) for r in ok), max_external_refold_error_kcal=max(abs(r["mfe_kcal"] - r["external_d0"]["mfe_kcal"]) for r in ok), precision_pairs=precision_pairs, ensemble_ranking_inversions=inversion, timing_boundary="One instrumented validation run; elapsed includes lattice, solve, traceback and candidate statistics; reported RSS sampled before optional dense checks. Never a controlled performance ratio.")
    dump(out / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--source", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--ward", type=Path, help="Optional pinned mrna_folding_comparison checkout; required only when generating fresh inputs")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--skip-designs", action="store_true")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--lock", type=Path, help="Shared flock held during all computations")
    args = parser.parse_args()
    root = args.root.resolve()
    out = (args.out or root / "research/results/publication-2026-09-06/validation").resolve()
    source = (args.source or root / "lineardesign-clean").resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.prepare:
        prepare_inputs(root, out, args.ward)
    if args.build:
        build(source, out)
    import RNA
    if RNA.__version__ != "2.7.2":
        raise RuntimeError("Pinned validation requires ViennaRNA==2.7.2")
    RNA.params_load_RNA_Turner2004()
    inputs = json.loads((out / "inputs.json").read_text())
    for key, rel in inputs["tables"].items():
        if digest(out / rel) != inputs["table_sha256"][key]:
            raise RuntimeError("table hash drift: " + key)
    lock = args.lock or out.parent / "compute.lock"
    lock.parent.mkdir(parents=True, exist_ok=True)
    with lock.open("a") as handle:
        print("Waiting for validation compute lock", flush=True)
        fcntl.flock(handle, fcntl.LOCK_EX)
        print("Acquired validation compute lock", flush=True)
        binary = out / "validation-driver"
        output = out / "designs.jsonl"
        previous = [json.loads(line) for line in output.read_text().splitlines()] if output.exists() else []
        completed = {r["case"]["id"] for r in previous}
        if not args.skip_designs:
            with output.open("a") as stream:
                for case in inputs["cases"]:
                    if case["id"] in completed:
                        continue
                    value = run_driver(binary, "design", out / inputs["tables"][case["table"]], case["protein"], case["lambda_user"], case["motifs"], args.timeout)
                    value = verify(case, value, out / inputs["tables"][case["table"]], ensemble=case["ensemble"])
                    value["case"] = case
                    value["binary_sha256"] = digest(binary)
                    stream.write(json.dumps(value, sort_keys=True) + "\n")
                    stream.flush()
                    print(case["id"], value["status"], "PASS" if value["checks_pass"] else "FAIL", flush=True)
        if not (out / "exhaustive.json").exists():
            exhaustive(binary, out / inputs["tables"]["human"], out)
        if not (out / "loops.json").exists():
            loop_validation(binary, out / inputs["tables"]["human"], out)
        result = summary(out)
        print(json.dumps(result, indent=2), flush=True)
        return int(bool(result["failures"]) or result["loops"]["passed"] != result["loops"]["cases"] or result["exhaustive"]["passed"] != result["exhaustive"]["comparisons"])


if __name__ == "__main__":
    raise SystemExit(main())
