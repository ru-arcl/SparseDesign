#!/usr/bin/env python3
"""Reproducible serialized publication timings; see --help and protocol.json.

Uses a separate process for each solve, GNU time for peak RSS, CPU affinity,
fixed OpenMP settings, a common advisory compute lock, and hash-checked resume.
All failures remain explicit records. Profiling is never mixed with timings.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import resource
import signal
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT if (ROOT / 'src/fold_turner.cc').exists() else ROOT / 'lineardesign-clean'
DEFAULT = ROOT / 'research/results/publication-2026-09-08'
MODES = ('dense_typed', 'dense_scalar', 'sparse')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def parse_kv(text):
    result = {}
    for item in text.strip().split():
        if '=' in item:
            k, v = item.split('=', 1)
            try:
                result[k] = float(v) if any(c in v for c in '.eE') else int(v)
            except ValueError:
                result[k] = v
    return result


def proc_stats():
    return {'load': list(os.getloadavg()),
            'meminfo': Path('/proc/meminfo').read_text(),
            'stat': Path('/proc/stat').read_text(),
            'vmstat': Path('/proc/vmstat').read_text()}


def make_protocol(output):
    panel = json.loads((output / 'inputs/panel.json').read_text())['panel']
    rng = random.Random(20260906)
    tasks = []
    def add(phase, item, lam, threads, mode, repeat, layout='square'):
        name = f"{phase}-{item['accession']}-l{lam}-j{threads}-{mode}-{layout}-r{repeat}"
        tasks.append({'id': name, 'phase': phase, 'accession': item['accession'],
                      'aa_length': item['aa_length'], 'fasta': item['fasta_path'],
                      'fasta_sha256': item['fasta_sha256'], 'lambda': lam,
                      'threads': threads, 'mode': mode, 'repeat': repeat, 'layout': layout})
    for item in [panel[i] for i in (0, 2, 6, 10)]:
        for threads in (1, 16):
            for mode in MODES:
                add('pilot', item, 0, threads, mode, 0)
    conditions = list(itertools.product(panel, (0, 4), (1, 16)))
    rng.shuffle(conditions)
    permutations = list(itertools.permutations(MODES))
    for item, lam, threads in conditions:
        order = rng.randrange(6)
        for repeat in range(1, 6):
            for mode in permutations[(order + repeat - 1) % 6]:
                add('ablation', item, lam, threads, mode, repeat)
    for item in panel:
        for lam in (0, 4):
            for threads in (1, 16):
                for mode in MODES:
                    add('profile', item, lam, threads, mode, 1)
    # Complete-CLI repeated strong scaling on one medium and one long natural protein.
    for item in (panel[6], panel[10]):
        conditions = list(itertools.product((0, 4), (1, 2, 4, 8, 16), range(1, 4)))
        rng.shuffle(conditions)
        for lam, threads, repeat in conditions:
            add('scaling', item, lam, threads, 'cli', repeat)
    dmd = {'accession': 'NP_000100.3', 'aa_length': 3677,
           'fasta_path': 'inputs/NP_000100.3.fasta',
           'fasta_sha256': sha(output / 'inputs/NP_000100.3.fasta')}
    # Five paired full-length layout repetitions, not best-of-five endpoints.
    for repeat in range(1, 6):
        for layout in (('square', 'packed') if repeat % 2 else ('packed', 'square')):
            add('workstation', dmd, 0, 16, 'cli', repeat, layout)
    return {'schema': 'sparsedesign-timing-protocol-v1', 'created_utc': utc(),
            'panel_sha256': sha(output / 'inputs/panel.json'), 'seed': 20260906,
            'address_space_limit_gib': 32, 'timeout_seconds': 600,
            'workstation_timeout_seconds': 1800,
            'openmp': {'OMP_DYNAMIC': 'false', 'OMP_PROC_BIND': 'close',
                       'OMP_PLACES': 'cores', 'OMP_WAIT_POLICY': 'active', 'OMP_THREAD_LIMIT': '16'},
            'affinity': 'nested logical CPUs 0..threads-1; i9 hybrid P/E topology recorded, not homogeneous cores',
            'timing_boundary': 'Cost calls: solver construction, fill, destruction; excludes DFA and traceback. CLI: complete process. Profile separate.',
            'interference_control': 'One cooperating evaluation process at a time via compute.lock; host load and proc snapshots retained.',
            'tasks': tasks}


def run_task(task, protocol, output, binaries):
    location = output / 'timings/runs' / task['id']
    record_path = location / 'record.json'
    binary = binaries['bench'] if task['mode'] != 'cli' else binaries[task['layout']]
    binary_digest = sha(binary)
    if record_path.exists():
        old = json.loads(record_path.read_text())
        if old['task'] != task or old['binary_sha256'] != binary_digest:
            raise RuntimeError(f"resume identity mismatch: {task['id']}")
        for filename, digest in old['output_sha256'].items():
            if sha(location / filename) != digest:
                raise RuntimeError(f"resume output mismatch: {task['id']} {filename}")
        return old
    location.mkdir(parents=True, exist_ok=True)
    fasta = output / task['fasta']
    if sha(fasta) != task['fasta_sha256']:
        raise RuntimeError('input checksum mismatch')
    table = SOLVER / 'data/codon_usage_freq_table_human.csv'
    if task['mode'] == 'cli':
        command = [str(binary), '-l', str(task['lambda']), '-j', str(task['threads']), '-c', str(table), '--diagnostics']
    else:
        command = [str(binary), task['mode'], str(task['lambda']), str(task['threads']), str(table), str(fasta)]
        if task['phase'] == 'profile':
            command.append('profile')
    cpus = list(range(task['threads']))
    invocation = ['/usr/bin/time', '-v', '-o', str(location / 'time.txt'),
                  'taskset', '-c', ','.join(map(str, cpus)), *command]
    env = dict(os.environ, **protocol['openmp'])
    def limits():
        cap = protocol['address_space_limit_gib'] * 2**30
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    with (output / 'compute.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        before, started = proc_stats(), utc()
        clock_start = time.perf_counter()
        with fasta.open('rb') as stdin, (location / 'stdout.txt').open('wb') as stdout, (location / 'stderr.txt').open('wb') as stderr:
            process = subprocess.Popen(invocation, cwd=SOLVER, env=env, stdin=stdin,
                                       stdout=stdout, stderr=stderr, start_new_session=True, preexec_fn=limits)
            timed_out = False
            try:
                cap_seconds = protocol.get('workstation_timeout_seconds', protocol['timeout_seconds']) if task['phase'] == 'workstation' else protocol['timeout_seconds']
                code = process.wait(timeout=cap_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                code = process.wait()
        wall = time.perf_counter() - clock_start
        after = proc_stats()
        fcntl.flock(lock, fcntl.LOCK_UN)
    stdout = (location / 'stdout.txt').read_text()
    stderr = (location / 'stderr.txt').read_text()
    timer = (location / 'time.txt').read_text() if (location / 'time.txt').exists() else ''
    time_values = {}
    for line in timer.splitlines():
        if ': ' in line:
            k, v = line.strip().rsplit(': ', 1)
            time_values[k] = v
    metrics = parse_kv(stdout) if task['mode'] != 'cli' else {}
    if task['mode'] == 'cli':
        for line in stderr.splitlines():
            if 'dp_cost=' in line:
                metrics.update(parse_kv(line))
        if 'dp_cost' in metrics:
            metrics['objective_units'] = metrics['dp_cost']
    status = 'timeout' if timed_out else ('ok' if code == 0 else 'failed')
    if status == 'ok' and task['mode'] != 'cli':
        if not math.isfinite(metrics.get('objective_units', float('nan'))) or metrics.get('threads_effective') != task['threads']:
            status = 'invalid_output'
    if status == 'ok' and task['mode'] == 'cli':
        if not math.isfinite(metrics.get('objective_units', float('nan'))) or abs(metrics.get('objective_delta', float('inf'))) > 1e-4:
            status = 'invalid_output'
    record = {'task': task, 'status': status, 'exit_code': code, 'started_utc': started,
              'ended_utc': utc(), 'wall_seconds': wall, 'binary_sha256': binary_digest,
              'codon_table_sha256': sha(table), 'command': invocation,
              'metrics': metrics, 'gnu_time': time_values, 'before': before, 'after': after,
              'output_sha256': {name: sha(location / name) for name in ('stdout.txt', 'stderr.txt', 'time.txt') if (location / name).exists()}}
    save(record_path, record)
    print(f"{task['id']} {status} wall={wall:.3f}s", flush=True)
    return record


def summarize(output, protocol):
    groups, records = {}, []
    for task in protocol['tasks']:
        path = output / 'timings/runs' / task['id'] / 'record.json'
        if not path.exists():
            continue
        r = json.loads(path.read_text()); records.append(r)
        groups.setdefault(task['phase'], []).append(r)
    pairs = {}
    for r in groups.get('ablation', []):
        t = r['task']; key = (t['accession'], t['lambda'], t['threads'], t['repeat'])
        pairs.setdefault(key, {})[t['mode']] = r
    cells = {}
    for key, modes in pairs.items():
        if len(modes) != 3 or any(r['status'] != 'ok' for r in modes.values()):
            continue
        values = [r['metrics']['objective_units'] for r in modes.values()]
        if max(values) - min(values) > 1e-4:
            raise RuntimeError(f'objective disagreement: {key} {values}')
        dense, scalar, sparse = [modes[m]['metrics']['solve_seconds'] for m in MODES]
        cells.setdefault(key[:3], []).append({'repeat': key[3], 'typed_over_sparse': dense/sparse,
                                            'scalar_over_sparse': scalar/sparse,
                                            'typed_over_scalar': dense/scalar})
    summaries = []
    for (accession, lam, threads), items in sorted(cells.items()):
        summaries.append({'accession': accession, 'lambda': lam, 'threads': threads,
                          'pairs': len(items), 'observations': items,
                          **{name: statistics.median(x[name] for x in items)
                             for name in ('typed_over_sparse','scalar_over_sparse','typed_over_scalar')}})
    result = {'created_utc': utc(), 'protocol_sha256': sha(output/'protocol.json'),
              'phase_completion': {phase: {'expected': sum(t['phase'] == phase for t in protocol['tasks']),
                                          'recorded': len(rs), 'ok': sum(r['status'] == 'ok' for r in rs)} for phase, rs in groups.items()},
              'ablation': summaries, 'failures': [{'task': r['task']['id'], 'status': r['status']} for r in records if r['status'] != 'ok']}
    save(output / 'timings/summary.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('plan', 'run', 'summarize'))
    p.add_argument('--output', type=Path, default=DEFAULT)
    p.add_argument('--phase', choices=('pilot','ablation','profile','scaling','workstation'), default='pilot')
    p.add_argument('--limit', type=int)
    args = p.parse_args(); output = args.output.resolve()
    protocol_path = output / 'protocol.json'
    if args.action == 'plan':
        if protocol_path.exists():
            raise SystemExit('Protocol already exists; do not silently replace a running study.')
        save(protocol_path, make_protocol(output))
        return
    protocol = json.loads(protocol_path.read_text())
    if sha(output / 'inputs/panel.json') != protocol['panel_sha256']:
        raise RuntimeError('panel changed since protocol freeze')
    if args.action == 'run':
        binaries = {'bench': output/'bin/publication-bench', 'square': output/'bin/sparsedesign-square', 'packed': output/'bin/sparsedesign-packed'}
        tasks = [t for t in protocol['tasks'] if t['phase'] == args.phase]
        if args.limit is not None:
            tasks = tasks[:args.limit]
        for task in tasks:
            run_task(task, protocol, output, binaries)
            summarize(output, protocol)
    else:
        print(json.dumps(summarize(output, protocol), indent=2))


if __name__ == '__main__':
    main()
