#!/usr/bin/env python3
"""Verify and summarize the separate EPYC Q8VIM6 follow-up; no new solves."""
import hashlib
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent / 'followup-Q8VIM6'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    build = json.loads((ROOT / 'build.json').read_text())
    for name, digest in build['source_sha256'].items():
        assert sha(ROOT / 'source' / name) == digest, name
    tasks = protocol['tasks']
    assert len(tasks) == 15
    assert len({t['id'] for t in tasks}) == 15
    blocks, records = {}, {}
    for task in tasks:
        assert (task['accession'], task['lambda'], task['threads']) == ('Q8VIM6', 4, 1)
        directory = ROOT / 'timings/runs' / task['id']
        record = json.loads((directory / 'record.json').read_text())
        assert record['task'] == task and record['status'] == 'ok'
        assert record['protocol_sha256'] == sha(ROOT / 'protocol.json')
        assert record['binary_sha256'] == build['binary_sha256']['publication-bench']
        assert record['collector_sha256'] == build['source_sha256']['publication_campaign_final.py']
        assert record['codon_table_sha256'] == build['source_sha256']['data/codon_usage_freq_table_human.csv']
        assert sha(ROOT / task['fasta']) == task['fasta_sha256']
        for name, digest in record['output_sha256'].items():
            assert sha(directory / name) == digest, (task['id'], name)
        blocks.setdefault(task['repeat'], {})[task['mode']] = record['metrics']
        records[task['id']] = sha(directory / 'record.json')
    assert set(blocks) == set(range(1, 6))
    assert all(set(b) == {'dense_scalar', 'dense_typed', 'sparse'} for b in blocks.values())
    objectives = [m['objective_units'] for b in blocks.values() for m in b.values()]
    assert max(objectives) - min(objectives) <= 1e-4
    ratios = {}
    for a, b in [('dense_scalar', 'sparse'), ('dense_typed', 'sparse'), ('dense_typed', 'dense_scalar')]:
        values = [blocks[r][a]['solve_seconds'] / blocks[r][b]['solve_seconds'] for r in sorted(blocks)]
        ratios[a + '_over_' + b] = dict(paired_values=values, median=statistics.median(values),
                                        min=min(values), max=max(values))
    summary = dict(platform='AMD EPYC 7313', successful=15, complete_triplets=5,
                   protocol_sha256=sha(ROOT / 'protocol.json'), build_sha256=sha(ROOT / 'build.json'),
                   generator_sha256=sha(Path(__file__)), record_sha256=records, ratios=ratios,
                   interpretation='Separate descriptive follow-up; never pooled with primary ablation. '
                   'The original i9 quiet-host gate was not applied to these EPYC measurements.')
    (ROOT / 'analysis').mkdir(exist_ok=True)
    (ROOT / 'analysis/summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps(ratios, indent=2))


if __name__ == '__main__':
    main()
