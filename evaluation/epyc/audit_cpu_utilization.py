#!/usr/bin/env python3
"""Audit whole-host mean CPU use from kept EPYC run snapshots, without rerunning solves."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def utilization(record):
    snapshots = [record.get(stage, {}).get('stat') for stage in ('before', 'after')]
    if not all(snapshots):
        return None
    counters = []
    for snapshot in snapshots:
        lines = snapshot.splitlines()
        assert lines[0].split()[0] == 'cpu'
        assert sum(line.split()[0][3:].isdigit() for line in lines if line.startswith('cpu')) == 32
        counters.append([int(value) for value in lines[0].split()[1:9]])
    delta = [after - before for before, after in zip(*counters)]
    assert len(delta) == 8 and all(value >= 0 for value in delta) and sum(delta) > 0
    return 100 * sum(delta[index] for index in (0, 1, 2, 5, 6, 7)) / sum(delta)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    sources = ('timings/timings/runs', 'followup-Q8VIM6/timings/runs',
               'public-baselines/feasibility', 'public-baselines/repeated')
    for source in sources:
        for path in sorted((args.input / source).glob('*/record.json')):
            raw = path.read_bytes()
            record = json.loads(raw)
            task = record.get('task', {})
            rows.append(dict(source=source, id=path.parent.name,
                             phase=task.get('phase', source.split('/')[-1]),
                             threads=task.get('threads', 1), status=record['status'],
                             attempted=record['status'] != 'not_attempted_after_cap',
                             host_cpu_percent=utilization(record),
                             record_sha256=sha(raw)))
    groups = {}
    for source in sources:
        selected = [row for row in rows if row['source'] == source]
        measured = [row for row in selected if row['host_cpu_percent'] is not None]
        non16 = [row for row in measured if row['threads'] != 16]
        groups[source] = dict(records=len(selected), statuses=dict(Counter(row['status'] for row in selected)),
                              attempted=sum(row['attempted'] for row in selected),
                              snapshots=len(measured), non16_snapshots=len(non16),
                              non16_max_cpu_percent=max((row['host_cpu_percent'] for row in non16), default=None),
                              non16_over90_ids=[row['id'] for row in non16 if row['host_cpu_percent'] > 90],
                              attempted_without_snapshots=sum(row['attempted'] and row['host_cpu_percent'] is None
                                                              for row in selected))
    paper_non16 = [row for row in rows if row['host_cpu_percent'] is not None
                   and row['threads'] != 16 and row['phase'] != 'pilot']
    validation_path = args.input / 'validation/designs.jsonl'
    validation_raw = validation_path.read_bytes()
    validation = [json.loads(line) for line in validation_raw.splitlines()]
    result = dict(schema='sparsedesign-epyc-cpu-audit-v1',
                  script_sha256=sha(Path(__file__).read_bytes()), groups=groups,
                  definition='Run-average busy / total delta from the aggregate cpu line in before/after /proc/stat; '
                             '32 CPUs including our processes. Busy=user+nice+system+irq+softirq+steal; '
                             'total=busy+idle+iowait. Guest counters are not added a second time.',
                  paper_timing_non16=dict(records=len(paper_non16),
                                          maximum_cpu_percent=max(row['host_cpu_percent'] for row in paper_non16),
                                          over90_ids=[row['id'] for row in paper_non16 if row['host_cpu_percent'] > 90]),
                  validation=dict(designs=len(validation), designs_sha256=sha(validation_raw),
                                  snapshots=sum(utilization(row) is not None for row in validation)),
                  limitations=['Native records have starting load averages but no per-run CPU snapshots; '
                               'load averages cannot establish utilization.',
                               'Validation design records have no before/after CPU snapshots.',
                               'Snapshot deltas establish run averages, not instantaneous peak CPU use.',
                               'Do not extend the timing-suite CPU claim to native comparisons or validation.'],
                  records=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
