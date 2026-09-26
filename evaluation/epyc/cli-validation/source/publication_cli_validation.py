#!/usr/bin/env python3
"""Translate and independently score the 70 planned controlled CLI outcomes.

ViennaRNA 2.7.2 evaluation-only mode; no design, MFE or partition-function call.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
import analyze_publication_campaign as strict
import analyze_publication_baselines as baseline
import publication_baselines as native


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def load_weights(path):
    frequencies, translation = {}, {}
    with path.open() as handle:
        for row in csv.reader(handle):
            if len(row) == 3 and row[0] in native.STANDARD_CODE:
                require(row[0] not in frequencies, 'Duplicate codon')
                frequencies[row[0]], translation[row[0]] = float(row[2]), row[1]
    require(translation == native.STANDARD_CODE, 'Expected complete standard genetic code')
    require(all(math.isfinite(v) and v > 0 for v in frequencies.values()), 'Invalid codon frequency')
    maxima = {aa: max(frequencies[c] for c in frequencies if translation[c] == aa)
              for aa in set(translation.values())}
    return {c: v / maxima[translation[c]] for c, v in frequencies.items()}


def parse_cli(text):
    values = [re.findall(pattern, text, re.M) for pattern in (
        r'^mRNA sequence:\s+([ACGU]+)$', r'^mRNA structure:\s+([.()]+)$',
        r'^mRNA folding free energy: (-?\d+\.\d{2}) kcal/mol; mRNA CAI: (\d+\.\d{3})$')]
    require(all(len(v) == 1 for v in values), 'Expected exactly one complete CLI result')
    rna, structure, energy = (v[0] for v in values)
    return dict(rna=rna, structure=structure,
                printed_energy_centikcal=int(Decimal(energy[0]) * 100))


def configure_rna(RNA):
    require(RNA.__version__ == '2.7.2', 'Pinned ViennaRNA 2.7.2 required')
    require(int(RNA.MAXLOOP) == 30, 'Unexpected ViennaRNA two-loop cap')
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    settings = dict(temperature=37., dangles=0, special_hp=1, noLP=0, noGU=0,
                    noGUclosure=0, min_loop_size=3, gquad=0, logML=0, circ=0, salt=1.021)
    for key, value in settings.items():
        setattr(md, key, value)
    model = dict(viennarna=RNA.__version__, parameter_set='Turner2004',
                 parameter_string_sha256=hashlib.sha256(RNA.parameter_set_rna_turner2004.encode()).hexdigest(),
                 settings=settings, internal_loop_unpaired_cap=30,
                 constructor_option='OPTION_EVAL_ONLY', scorer='eval_structure_pt',
                 energy_unit='integer centikcal/mol (10 cal/mol)',
                 mfe_calls=0, partition_function_calls=0)
    return md, model


def score_pair(RNA, md, rna, structure):
    compound = RNA.fold_compound(rna, md, RNA.OPTION_EVAL_ONLY)
    require(compound.matrices is None and compound.exp_matrices is None,
            'Evaluation-only compound unexpectedly allocated DP matrices')
    energy = compound.eval_structure_pt(RNA.ptable(structure))
    require(isinstance(energy, int), 'Expected integer structure energy')
    return dict(rna=rna, structure=structure, energy_centikcal=energy,
                mfe_matrices_allocated=False, partition_matrices_allocated=False)


def check_result(parsed, protein, weights, metrics, lam, score):
    validation = native.validate_output(parsed, protein, weights)
    require(all(validation['validation_checks'].values()), 'Invalid translated RNA/structure')
    loops = baseline.structure_metrics(parsed['structure'])
    energy = score['energy_centikcal']
    penalty = validation['common_codon_penalty']
    objective = energy + 100 * lam * penalty
    checks = dict(validation['validation_checks'])
    checks.update(two_loop_cap30=loops['max_two_loop_unpaired'] <= 30,
                  minimum_hairpin=loops['min_hairpin_unpaired'] is None or loops['min_hairpin_unpaired'] >= 3,
                  printed_energy_exact=parsed['printed_energy_centikcal'] == energy,
                  weighted_dp_objective=abs(objective - metrics['dp_cost']) <= 1e-4,
                  weighted_realized_objective=abs(objective - metrics['realized_cost']) <= 1e-4,
                  codon_penalty=abs(penalty - metrics['codon_penalty']) <= 1e-8,
                  cai=abs(validation['common_cai'] - metrics['cai']) <= 1e-12)
    if lam == 0:
        checks['lambda0_dp_exact_integer'] = metrics['dp_cost'] == energy
        checks['lambda0_realized_exact_source_roundtrip'] = metrics['realized_cost'] == (energy / 100.0) * 100.0
    return dict(checks=checks, all_checks_pass=all(checks.values()), structure_metrics=loops,
                energy_centikcal=energy, independent_codon_penalty=penalty,
                independent_cai=validation['common_cai'], independent_objective_centikcal=objective,
                dp_difference_centikcal=objective - metrics['dp_cost'],
                realized_difference_centikcal=objective - metrics['realized_cost'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--lock', type=Path, required=True)
    parser.add_argument('--allow-partial', action='store_true', help='Explicitly provisional output only')
    args = parser.parse_args()
    require(not args.output.exists(), 'Use a new validation output directory')
    protocol, build, records, missing = strict.verified_records(args.input, not args.allow_partial)
    planned = [t for t in protocol['tasks'] if t['phase'] in ('scaling', 'workstation')]
    require(Counter(t['phase'] for t in planned) == {'scaling': 60, 'workstation': 10}
            and all(t['mode'] == 'cli' for t in planned), 'Unexpected planned CLI grid')
    by_id = {r['task']['id']: r for r in records}
    table = args.input / 'source/data/codon_usage_freq_table_human.csv'
    require(sha(table) == build['source_sha256']['data/codon_usage_freq_table_human.csv'], 'Codon table changed')
    weights = load_weights(table)
    cache, rows = {}, []
    with strict.analysis_lock(args.lock, True):
        import RNA
        md, model = configure_rna(RNA)
        module_files = [Path(RNA.__file__)]
        extension = getattr(RNA, '_RNA', None)
        if extension is not None and getattr(extension, '__file__', None):
            module_files.append(Path(extension.__file__))
        model['installed_module_sha256'] = {p.name: sha(p) for p in module_files}
        for task in planned:
            record = by_id.get(task['id'])
            row = dict(task=task, status=record['status'] if record else 'missing',
                       record_sha256=record['_record_sha256'] if record else None)
            if record and record['status'] == 'ok':
                path = args.input / 'timings/runs' / task['id']
                parsed = parse_cli((path / 'stdout.txt').read_text())
                fasta = args.input / task['fasta']
                require(sha(fasta) == task['fasta_sha256'], 'Protein input changed')
                protein = ''.join(line.strip() for line in fasta.read_text().splitlines() if not line.startswith('>'))
                key = hashlib.sha256((parsed['rna'] + '\n' + parsed['structure']).encode()).hexdigest()
                if key not in cache:
                    cache[key] = score_pair(RNA, md, parsed['rna'], parsed['structure'])
                result = check_result(parsed, protein, weights, record['metrics'], task['lambda'], cache[key])
                row.update(score_cache_key=key, raw_file_sha256=record['output_sha256'], **result)
            rows.append(row)
    successful = [r for r in rows if r['status'] == 'ok']
    summary = dict(schema='sparsedesign-controlled-cli-validation-v1', planned=70,
                   recorded=sum(r['status'] != 'missing' for r in rows),
                   original_study_complete=not missing, statuses=dict(Counter(r['status'] for r in rows)),
                   successful_outputs_checked=len(successful), distinct_scored_rna_structure_pairs=len(cache),
                   all_successful_outputs_pass=all(r['all_checks_pass'] for r in successful) if successful else None,
                   failed_check_task_ids=[r['task']['id'] for r in successful if not r['all_checks_pass']],
                   maximum_absolute_dp_difference_centikcal=max((abs(r['dp_difference_centikcal']) for r in successful), default=None),
                   maximum_absolute_realized_difference_centikcal=max((abs(r['realized_difference_centikcal']) for r in successful), default=None),
                   protocol_sha256=sha(args.input / 'protocol.json'), build_sha256=sha(args.input / 'build.json'),
                   codon_table_sha256=sha(table), model=model, outcomes=rows,
                   source_sha256={Path(p).name: sha(p) for p in (__file__, strict.__file__, baseline.__file__, native.__file__)},
                   interpretation='Independent translation and returned-structure energy/admissibility checks only. No fixed-sequence refolding, new sequence design, ensemble claim or independent optimality proof.')
    args.output.mkdir(parents=True)
    save(args.output / 'summary.json', summary)
    save(args.output / 'score-cache.json', dict(model=model, values=cache))
    save(args.output / 'input-records.json', dict(protocol_sha256=summary['protocol_sha256'], build_sha256=summary['build_sha256'],
                                               record_sha256={r['task']['id']: r['record_sha256'] for r in rows if r['record_sha256']}))
    report = ['# Independent checks of controlled CLI outputs', '', summary['interpretation'], '',
              f"Original study complete: {summary['original_study_complete']}. CLI outcomes: {summary['recorded']}/70; statuses: {summary['statuses']}.",
              f"Checked {len(successful)} successful executions and {len(cache)} distinct RNA/structure pairs; every successful output passes: {summary['all_successful_outputs_pass']}.", '',
              'The frozen raw-record verifier runs first. Translation uses the standard genetic code; codon penalty and CAI are recomputed from the original hashed table. Structures must have canonical pairs, hairpins of at least three unpaired bases, and at most 30 unpaired bases in each bulge/internal loop.', '',
              'ViennaRNA 2.7.2 loads Turner 2004 at 37 C, dangles 0, special hairpins enabled and salt 1.021 M. OPTION_EVAL_ONLY must leave both DP matrix pointers empty. eval_structure_pt returns integer centikcal/mol. Printed energies must match exactly; lambda-0 DP energies must equal that integer and realized diagnostics must equal the exact source divide/multiply round trip. Weighted objectives allow 1e-4 centikcal/mol numerical tolerance; codon penalty and CAI tolerances are 1e-8 and 1e-12. All individual checks and differences remain in JSON.', '',
              '| Phase | Accession | Lambda | Threads | Layout | Repeat | Status | Checks pass |',
              '|---|---|---:|---:|---|---:|---|---|']
    for row in rows:
        t = row['task']
        report.append(f"| {t['phase']} | {t['accession']} | {t['lambda']} | {t['threads']} | {t['layout']} | {t['repeat']} | {row['status']} | {row.get('all_checks_pass', 'not applicable')} |")
    (args.output / 'REPORT.md').write_text('\n'.join(report) + '\n')
    (args.output / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(args.output.iterdir()) if p.is_file()))
    print(json.dumps({k: summary[k] for k in ('recorded', 'statuses', 'successful_outputs_checked', 'distinct_scored_rna_structure_pairs', 'all_successful_outputs_pass')}))
    return 0 if not summary['failed_check_task_ids'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
