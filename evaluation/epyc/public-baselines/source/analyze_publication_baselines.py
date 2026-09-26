#!/usr/bin/env python3
"""Independently translate and rescore retained native baseline outputs.

Run with ViennaRNA 2.7.2. This program never redesigns proteins. Optional fixed
sequence refolding and all structure rescoring occur outside native timing,
under the same advisory compute lock. Native outcomes remain immutable.
"""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

from publication_baselines import STANDARD_CODE, parse_output, sha, validate_output, write_json


def task_key(record):
    return (record['tool'],record['accession'],record['lambda_value'],record['repetition'])


def centikcal(kcal):
    # The audited model sums integer centikcal terms; ViennaRNA's Python API
    # returns a floating kcal value. Recover the exact model-energy quantum.
    return int(round(100*kcal))


def sparse_lambda0_objective_match(row):
    energy = centikcal(row['structure_energy_kcal'])
    # The native DP sums integer centikcal terms at lambda 0. Its independent
    # diagnostic explicitly divides the integer rescore by 100.0, then multiplies
    # by 100.0 (fold_turner.cc). Require that exact double round trip, which can
    # differ from the integer by one floating-point ULP; no energy tolerance.
    return (row.get('diagnostic_dp_cost') == energy
            and row.get('diagnostic_realized_cost') == (energy / 100.0) * 100.0)


def manifest_completion(planned, records):
    expected = Counter(task_key(r) for r in planned)
    observed = Counter(task_key(r) for r in records)
    if any(n != 1 for n in expected.values()):
        raise ValueError('Planned task keys are not unique')
    if any(n != 1 for n in observed.values()) or observed-expected:
        raise ValueError('Observed tasks include duplicates or tasks outside the manifest')
    return dict(planned=len(planned),recorded=len(records),complete=observed==expected,
                missing_task_keys=[list(k) for k in sorted(expected-observed)])


def compare_phase_metadata(first, second):
    for key in ('panel_sha256','table_sha256','derna_table_sha256','wrapper_sha256','cpu',
                'memory_limit_kind','memory_limit_bytes','timeout_seconds','shared_lock','lambdas','tools','seed'):
        if first[key] != second[key]:
            raise ValueError(f'Feasibility and repeated protocol differs: {key}')
    if set(first['pins']) != set(second['pins']):
        raise ValueError('Feasibility and repeated tool identities differ')
    changes = {}
    for tool in first['pins']:
        a,b = first['pins'][tool],second['pins'][tool]
        if tool != 'sparsedesign':
            if a != b:
                raise ValueError(f'Feasibility and repeated upstream provenance differs: {tool}')
            continue
        for key in ('url','commit','executable_sha256'):
            if a[key] != b[key]:
                raise ValueError(f'Feasibility and repeated SparseDesign identity differs: {key}')
        aa,bb = a['tracked_file_sha256'],b['tracked_file_sha256']
        if {k:v for k,v in aa.items() if k.startswith('src/')} != {k:v for k,v in bb.items() if k.startswith('src/')}:
            raise ValueError('Feasibility and repeated SparseDesign compiled src/ inputs differ')
        delta = dict(added=sorted(bb.keys()-aa.keys()),removed=sorted(aa.keys()-bb.keys()),
                     changed=sorted(k for k in aa.keys()&bb.keys() if aa[k]!=bb[k]))
        if any(delta.values()):
            changes[tool] = delta
    return changes


def structure_metrics(structure):
    stack, children, pairs, parent = [], defaultdict(list), {}, {}
    for i, ch in enumerate(structure):
        if ch == '(':
            par = stack[-1] if stack else -1
            children[par].append(i)
            parent[i] = par
            stack.append(i)
        elif ch == ')':
            pairs[stack.pop()] = i
    hairpins, two_loops, multi_flanks, multi_gaps = [], [], [], []
    for i, j in pairs.items():
        cc = children[i]
        unpaired = j-i-1-sum(pairs[k]-k+1 for k in cc)
        if not cc:
            hairpins.append(unpaired)
        elif len(cc) == 1:
            two_loops.append(unpaired)
        else:
            multi_flanks.append(cc[0]-i-1+j-pairs[cc[-1]]-1)
            multi_gaps.extend(cc[k+1]-pairs[cc[k]]-1 for k in range(len(cc)-1))
    return dict(min_hairpin_unpaired=min(hairpins, default=None),
                max_two_loop_unpaired=max(two_loops, default=0),
                max_multiloop_outer_flanks=max(multi_flanks, default=0),
                max_multiloop_interbranch_gap=max(multi_gaps, default=0))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runs', type=Path, nargs='+', required=True)
    ap.add_argument('--panel', type=Path, required=True)
    ap.add_argument('--table', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--lock', type=Path, required=True)
    ap.add_argument('--refold-max-nt', type=int, default=3000)
    ap.add_argument('--model-audit', type=Path, help='Required before generating conditional DERNA lambda0 timing ratios')
    args = ap.parse_args()
    import RNA
    if RNA.__version__ != '2.7.2':
        raise ValueError('The frozen common-scoring protocol requires ViennaRNA 2.7.2')
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature, md.dangles, md.special_hp, md.noLP, md.min_loop_size = 37., 0, 1, 0, 3
    md.salt = 1.021
    model = dict(viennarna=RNA.__version__, parameters='Turner2004', temperature=37,
                 dangles=0, special_hp=True, noLP=False, min_loop_size=3,
                 internal_loop_unpaired_cap=int(RNA.MAXLOOP), salt=md.salt,
                 parameter_string_sha256=hashlib.sha256(RNA.parameter_set_rna_turner2004.encode()).hexdigest())
    entries = {e['accession']: e for e in json.loads(args.panel.read_text())['panel']}
    proteins = {}
    for accession, entry in entries.items():
        fasta = args.panel.parent.parent/entry['fasta_path']
        if sha(fasta) != entry['fasta_sha256']:
            raise ValueError(f'FASTA digest mismatch: {fasta}')
        proteins[accession] = ''.join(x.strip() for x in fasta.read_text().splitlines() if not x.startswith('>'))
    freqs, translation = {}, {}
    for row in csv.reader(args.table.open()):
        if len(row) == 3 and row[0] in STANDARD_CODE:
            freqs[row[0]], translation[row[0]] = float(row[2]), row[1]
    if translation != STANDARD_CODE:
        raise ValueError('Table does not specify the standard genetic code')
    maxima = {aa: max(freqs[c] for c in freqs if translation[c] == aa) for aa in set(translation.values())}
    weights = {c: freqs[c]/maxima[translation[c]] for c in freqs}
    args.output.mkdir(parents=True, exist_ok=True)
    cache_path = args.output/'rescore-cache.json'
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else dict(model=model, values={})
    if cache['model'] != model:
        raise ValueError('Rescore model changed; use a new output directory')
    rows, groups, completion, final_metadata = [], defaultdict(list), {}, {}
    for folder in args.runs:
        metadata = json.loads((folder/'metadata.json').read_text())
        if metadata['panel_sha256'] != sha(args.panel) or metadata['table_sha256'] != sha(args.table):
            raise ValueError(f'Input provenance mismatch: {folder}')
        if folder.name in ('feasibility','repeated'):
            final_metadata[folder.name] = metadata
            if metadata['schema'] != 'sparsedesign-public-baselines-v2':
                raise ValueError('Final analysis requires the verified version-2 runner metadata')
        planned = json.loads((folder/'planned-order.json').read_text())
        native_records = [(p,json.loads(p.read_text())) for p in sorted(folder.glob('*/record.json'))]
        completion[folder.name] = manifest_completion(planned,[r for _,r in native_records])
        completion[folder.name].update(planned_manifest_sha256=sha(folder/'planned-order.json'),metadata_sha256=sha(folder/'metadata.json'))
        for record_path, native in native_records:
            if metadata['schema'] == 'sparsedesign-public-baselines-v2':
                for key in ('wrapper_sha256','table_sha256'):
                    if native[key] != metadata[key]:
                        raise ValueError(f'Native record provenance differs: {record_path}/{key}')
                if native['status'] != 'not_attempted_after_cap':
                    pin = metadata['pins'][native['tool']]
                    if native['executable_sha256'] != pin['executable_sha256'] or native['source_commit'] != pin['commit']:
                        raise ValueError(f'Native executable/source identity differs: {record_path}')
                    if native['input_fasta_sha256'] != entries[native['accession']]['fasta_sha256']:
                        raise ValueError(f'Native input identity differs: {record_path}')
            for name, expected in native.get('raw_file_sha256', {}).items():
                if sha(record_path.parent/name) != expected:
                    raise ValueError(f'Raw output digest mismatch: {record_path.parent/name}')
            accession, lam, tool = native['accession'], native['lambda_value'], native['tool']
            if accession not in entries:
                raise ValueError('A record is outside the supplied approved panel')
            expected_label = f"{accession}.lambda{lam:g}.{tool}.rep{native['repetition']}"
            if record_path.parent.name != native['label'] or native['label'] != expected_label or native['aa_length'] != entries[accession]['aa_length']:
                raise ValueError('Record identity or length differs from its approved task')
            row = {k: native.get(k) for k in ('label', 'tool', 'accession', 'aa_length', 'lambda_value', 'repetition', 'status', 'wall_seconds', 'user_seconds', 'system_seconds', 'max_rss_kib', 'native_temporary_bin_bytes')}
            row.update(phase=folder.name, record_sha256=sha(record_path), native_record=str(record_path))
            if native['status'] in ('ok', 'invalid_or_unparsed_output', 'invalid_output'):
                raw = (record_path.parent/'stdout.txt').read_text(errors='replace')
                result = record_path.parent/'native-result.txt'
                if result.exists():
                    raw += '\n'+result.read_text(errors='replace')
                parsed = parse_output(tool, raw)
                verification = validate_output(parsed, proteins[accession], weights)
                row.update(verification)
                row['native_mfe_kcal'] = parsed.get('native_mfe')
                row['validation_pass'] = all(verification['validation_checks'].values())
                if row['validation_pass']:
                    rna, structure = parsed['rna'], parsed['structure']
                    row.update(rna_sha256=hashlib.sha256(rna.encode()).hexdigest(), structure_sha256=hashlib.sha256(structure.encode()).hexdigest())
                    row.update(structure_metrics(structure))
                    key = hashlib.sha256((rna+'\n'+structure).encode()).hexdigest()
                    needs_refold = len(rna) <= args.refold_max_nt and 'fixed_sequence_mfe_kcal' not in cache['values'].get(key,{})
                    if key not in cache['values'] or needs_refold:
                        with args.lock.open('a') as lock:
                            fcntl.flock(lock, fcntl.LOCK_EX)
                            compound = RNA.fold_compound(rna, md)
                            score = dict(structure_energy_kcal=float(compound.eval_structure(structure)))
                            if len(rna) <= args.refold_max_nt:
                                _, mfe = compound.mfe()
                                score['fixed_sequence_mfe_kcal'] = float(mfe)
                            cache['values'][key] = score
                            fcntl.flock(lock, fcntl.LOCK_UN)
                        write_json(cache_path, cache)
                    row.update(cache['values'][key])
                    row['structure_energy_centikcal'] = centikcal(row['structure_energy_kcal'])
                    row['common_objective_kcal'] = row['structure_energy_centikcal']/100 + lam*row['common_codon_penalty']
                    row['admissible_cap30'] = row['max_two_loop_unpaired'] <= 30
                    if row['native_mfe_kcal'] is not None:
                        row['native_energy_minus_common_structure_kcal'] = row['native_mfe_kcal']-row['structure_energy_kcal']
                        row['native_energy_matches_common_centikcal'] = centikcal(row['native_mfe_kcal'])==centikcal(row['structure_energy_kcal'])
                    if 'fixed_sequence_mfe_kcal' in row:
                        row['traceback_minus_refold_kcal'] = row['structure_energy_kcal']-row['fixed_sequence_mfe_kcal']
                    if tool == 'sparsedesign':
                        stderr = (record_path.parent/'stderr.txt').read_text()
                        import re
                        for k,v in re.findall(r'\b(dp_cost|realized_cost|objective_delta|codon_penalty|cai)=([-+.eE0-9]+)', stderr):
                            row['diagnostic_'+k] = float(v)
                        if 'diagnostic_realized_cost' in row:
                            row['common_minus_realized_centikcal'] = 100*row['common_objective_kcal']-row['diagnostic_realized_cost']
                        if lam == 0:
                            row['sparsedesign_dp_matches_returned_centikcal'] = sparse_lambda0_objective_match(row)
            rows.append(row)
            groups[folder.name, tool, accession, lam].append(row)
    non_compilation_source_changes = {}
    if set(final_metadata) == {'feasibility','repeated'}:
        non_compilation_source_changes = compare_phase_metadata(final_metadata['feasibility'],final_metadata['repeated'])
        if final_metadata['repeated']['repetitions'] != 5:
            raise ValueError('The final repeated protocol requires five repetitions')
    summaries = []
    for (phase,tool,accession,lam), rr in sorted(groups.items()):
        ok = [r for r in rr if r['status'] == 'ok' and r.get('validation_pass')]
        summary = dict(phase=phase,tool=tool,accession=accession,lambda_value=lam,aa_length=entries[accession]['aa_length'],
                       outcomes=dict(Counter(r['status'] for r in rr)), successful=len(ok))
        for metric in ('wall_seconds','max_rss_kib','common_objective_kcal','structure_energy_kcal','common_cai'):
            values = [r[metric] for r in ok if r.get(metric) is not None]
            if values:
                summary[metric] = dict(median=statistics.median(values), minimum=min(values), maximum=max(values), n=len(values))
        summary['distinct_rna_outputs'] = len({r['rna_sha256'] for r in ok})
        summaries.append(summary)
    write_json(args.output/'results.json', rows)
    ratio_rows = []
    if args.model_audit:
        audit = json.loads(args.model_audit.read_text())
        if not audit['all_canonical_constants_equal']:
            raise ValueError('Parameter audit does not support matched thermodynamic ratios')
        if 'repeated' not in final_metadata:
            raise ValueError('Model-audited ratios require the repeated-run provenance')
        pins = final_metadata['repeated']['pins']
        if audit['pins']['derna'][1] != pins['derna']['commit']:
            raise ValueError('Audited DERNA revision differs from the measured revision')
        bound = {'derna':set(),'sparsedesign':set()}
        for path, digest in audit['source_sha256'].items():
            if '/derna/src/' in path:
                tool, relative = 'derna','src/'+path.split('/derna/src/',1)[1]
            elif '/src/vienna/' in path:
                tool, relative = 'sparsedesign','src/vienna/'+path.split('/src/vienna/',1)[1]
            else:
                raise ValueError(f'Unrecognized audited source identity: {path}')
            if pins[tool]['tracked_file_sha256'].get(relative) != digest:
                raise ValueError(f'Audited parameter source differs from measured source: {path}')
            bound[tool].add(relative)
        if bound['derna'] != {'src/default.cpp','src/params/intl11.h','src/params/intl21.h','src/params/intl22.h'} or bound['sparsedesign'] != {'src/vienna/energy_parameter.h','src/vienna/intl11.h','src/vienna/intl21.h','src/vienna/intl22.h'}:
            raise ValueError('Parameter audit does not cover the required thermodynamic sources')
        for accession in entries:
            rr = {tool:[r for r in rows if r['phase']=='repeated' and r['tool']==tool and r['accession']==accession and r['lambda_value']==0] for tool in ('derna','sparsedesign')}
            values = rr['derna']+rr['sparsedesign']
            eligible = all(len(v)==5 and {r['repetition'] for r in v}=={1,2,3,4,5} for v in rr.values()) and all(r['status']=='ok' and r.get('validation_pass') and r.get('admissible_cap30') and r.get('wall_seconds',0)>0 and r.get('native_energy_matches_common_centikcal') and (r['tool']!='sparsedesign' or r.get('sparsedesign_dp_matches_returned_centikcal')) for r in values)
            if eligible:
                energies = [r['structure_energy_kcal'] for r in values]
                eligible = len({centikcal(e) for e in energies})==1
            if eligible:
                ratio_rows.append(dict(accession=accession,aa_length=entries[accession]['aa_length'],
                     derna_median_seconds=statistics.median(r['wall_seconds'] for r in rr['derna']),
                     sparsedesign_median_seconds=statistics.median(r['wall_seconds'] for r in rr['sparsedesign']),
                     derna_over_sparsedesign=statistics.median(r['wall_seconds'] for r in rr['derna'])/statistics.median(r['wall_seconds'] for r in rr['sparsedesign']),
                     common_energy_kcal=energies[0],repetitions_each=5))
    write_json(args.output/'matched-lambda0-ratios.json',dict(model_audit_sha256=sha(args.model_audit) if args.model_audit else None,
               interpretation='Ratio of medians for common successful lambda0 cases with matching thermodynamic parameters and independently matching returned energies; complete native CLI boundary, conditional on completing both tools.',rows=ratio_rows))
    write_json(args.output/'summary.json', dict(model=model, panel_sha256=sha(args.panel), table_sha256=sha(args.table),
               script_sha256=sha(__file__), python=sys.version, groups=summaries, completion=completion, refold_max_nt=args.refold_max_nt,
               non_compilation_source_changes=non_compilation_source_changes,
               notes=['LD/LCD speed ratios are excluded because native grammars differ; DERNA lambda0 ratios require an explicit model audit and matching verified outputs.',
                      'Feasibility and pilot phases are excluded from repeated timing estimates.',
                      'GNU time reports centisecond wall resolution; zero observations are not valid ratio denominators.',
                      'DERNA lambda>0 reports native fixed-sequence refolding energy, which need not equal its traceback structure energy.']))
    flat = [{k:v for k,v in r.items() if not isinstance(v,dict)} for r in rows]
    keys = sorted(set().union(*(r.keys() for r in flat)))
    with (args.output/'results.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(flat)
    table_rows = []
    for s in summaries:
        if s['phase'] != 'repeated':
            continue
        table_rows.append(dict(tool=s['tool'],accession=s['accession'],aa_length=s['aa_length'],lambda_value=s['lambda_value'],
                              successful=s['successful'],wall_median=s.get('wall_seconds',{}).get('median'),
                              wall_min=s.get('wall_seconds',{}).get('minimum'),wall_max=s.get('wall_seconds',{}).get('maximum'),
                              rss_median_kib=s.get('max_rss_kib',{}).get('median'),objective_median_kcal=s.get('common_objective_kcal',{}).get('median')))
    if table_rows:
        with (args.output/'repeated-summary.csv').open('w') as handle:
            writer = csv.DictWriter(handle,fieldnames=list(table_rows[0]));writer.writeheader();writer.writerows(table_rows)
    print(json.dumps(dict(records=len(rows), validated=sum(r.get('validation_pass',False) for r in rows),
                         statuses=dict(Counter(r['status'] for r in rows)))))


if __name__ == '__main__':
    main()
