#!/usr/bin/env python3
"""Run timing-study tasks highest-thread-count first; optionally gate and redo runs slowed by contention.

By default each task runs once, in thread order, with no waiting and no retries. The load gate
and contention retry below apply only with --contention-control.

Uses the study's own, unchanged collector (imported, not copied), so every record keeps the
same collector and protocol hashes as a plain `publication_campaign_final.py run`.

Order: phases are merged and stably sorted by descending thread count, so tasks with equal
threads keep their protocol order and each ablation block stays adjacent in its registered
kernel order. The analyzer checks task identity and non-overlapping windows, not order.

Contention retry (user request, 2026-09-24): after each run, compare its wall time with the
same task on the i9. The EPYC is slower even when quiet, by an amount that depends on the
thread count, so a run is "slow" when

    epyc_wall / i9_wall > TOLERANCE * QUIET_RATIO[threads]

A slow run is moved to <study>/rejected-runs/<task>/attempt-N/ (never deleted) and the task is
rerun after a back-off, before moving on, until a run passes. Termination guard: once the last
3 attempts of a task agree within 10% (consistent, not contention-like), the task is treated as
intrinsically slower on this host and its fastest attempt is accepted and flagged. Every
decision is appended to <study>/retry-log.jsonl.

Load gate (user request, 2026-09-24): before every attempt, wait until the host is quiet enough
for this run. Host CPU activity is sampled from /proc/stat while none of our solves is running,
so it is other users' load; the gate passes when two consecutive GATE_WINDOW-second windows have
busy CPUs <= 32 - threads - GATE_HEADROOM (12 for a 16-thread run, 27 for a serial run). The
post-run ratio check above remains the acceptance test. Tasks without a successful i9 reference (the
i9 timeouts) are accepted as they are.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import time

PHASES = ('pilot', 'ablation', 'profile', 'scaling', 'workstation')
REPO = Path(__file__).resolve().parents[2]
I9_RUNS = REPO / 'publication/code/evaluation/timings/timings/runs'
I9_FOLLOWUP_RUNS = REPO / 'publication/code/evaluation/quiet-replication/runs'
# Quiet EPYC-node-1 / i9 wall ratios: j16 = median of the first 20 j16 ablation runs (load 12-26),
# j1 = median of the load-2 serial pilots. j2-j8 are interpolated in log2(threads).
QUIET_J1, QUIET_J16 = 2.17, 1.72
TOLERANCE = 1.3
CONSISTENT_ATTEMPTS, CONSISTENT_SPREAD = 3, 0.10
BACKOFF_SECONDS = (60, 120, 300, 600, 900)
HOST_CPUS, GATE_HEADROOM, GATE_WINDOW, GATE_WINDOWS, GATE_POLL = 32, 4, 15, 2, 30
BUSY_FIELDS = (0, 1, 2, 5, 6, 7)  # user nice system irq softirq steal


def busy_cpus(seconds):
    def snap():
        busy = total = 0
        for line in open('/proc/stat'):
            if line.startswith('cpu') and line[3] != ' ':
                v = [int(x) for x in line.split()[1:]]
                busy += sum(v[i] for i in BUSY_FIELDS); total += sum(v)
        return busy, total
    b0, t0 = snap(); time.sleep(seconds); b1, t1 = snap()
    return HOST_CPUS * (b1 - b0) / max(1, t1 - t0)


def wait_for_load(threads, log, task_id):
    allowed = HOST_CPUS - threads - GATE_HEADROOM
    passed, waited_since, last_note = 0, None, 0.0
    while passed < GATE_WINDOWS:
        busy = busy_cpus(GATE_WINDOW)
        if busy <= allowed:
            passed += 1
            continue
        passed = 0
        if waited_since is None:
            waited_since = time.time()
        if time.time() - last_note >= 600:
            last_note = time.time()
            log(task=task_id, decision='waiting_for_load', busy_cpus=round(busy, 2), allowed=allowed)
            print(f"{utc()} {task_id} waiting: other users busy {busy:.1f} CPUs > {allowed}", flush=True)
        time.sleep(GATE_POLL)
    if waited_since is not None:
        log(task=task_id, decision='load_ok', waited_seconds=round(time.time() - waited_since), allowed=allowed)


def quiet_ratio(threads):
    x = math.log2(threads) / 4
    return QUIET_J1 + (QUIET_J16 - QUIET_J1) * x


def utc():
    return datetime.now(timezone.utc).isoformat()


def i9_wall(task_id):
    for base in (I9_RUNS, I9_FOLLOWUP_RUNS):
        path = base / task_id / 'record.json'
        if path.exists():
            r = json.loads(path.read_text())
            return r['wall_seconds'] if r['status'] == 'ok' else None
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', type=Path, required=True, help='adapted timing (or follow-up) study directory')
    # Default excludes 'pilot': the paper reports only the ablation, profile, scaling and
    # workstation measurements (pilots are mentioned only as excluded).
    ap.add_argument('--phases', nargs='+', choices=PHASES, default=[p for p in PHASES if p != 'pilot'])
    ap.add_argument('--dry-run', action='store_true', help='print the remaining order and exit')
    # Off by default since 2026-09-24 (user: run the whole batch now, no waiting, no retries).
    ap.add_argument('--contention-control', action='store_true',
                    help='enable the pre-run load gate and the slow-run retry described above')
    args = ap.parse_args()
    output = args.output.resolve()
    path = output / 'source/publication_campaign_final.py'
    spec = importlib.util.spec_from_file_location('publication_campaign_final', path)
    collector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(collector)

    protocol = json.loads((output / 'protocol.json').read_text())
    if collector.sha(output / 'inputs/panel.json') != protocol['panel_sha256']:
        raise SystemExit('panel changed since protocol freeze')
    needed = set(range(collector.CPU_BASE, collector.CPU_BASE + 16))
    if not needed <= os.sched_getaffinity(0):
        raise SystemExit(f'CPUs {min(needed)}-{max(needed)} are not all available')
    binaries = {'bench': output / 'bin/publication-bench', 'square': output / 'bin/sparsedesign-square',
                'packed': output / 'bin/sparsedesign-packed'}
    runs, rejected = output / 'timings/runs', output / 'rejected-runs'

    indexed = [(i, t) for i, t in enumerate(protocol['tasks']) if t['phase'] in args.phases]
    ordered = [t for _, t in sorted(indexed, key=lambda it: (-it[1]['threads'], it[0]))]
    if args.dry_run:
        remaining = [t for t in ordered if not (runs / t['id'] / 'record.json').exists()]
        for t in remaining:
            print(t['threads'], f"limit={TOLERANCE * quiet_ratio(t['threads']):.2f}", t['id'])
        print(len(remaining), 'remaining of', len(ordered))
        return

    def log(**entry):
        with (output / 'retry-log.jsonl').open('a') as f:
            f.write(json.dumps(dict(utc=utc(), **entry)) + '\n')

    for task in ordered:
        if not args.contention_control:
            collector.run_task(task, protocol, output, binaries)  # one run per task; verifies/skips kept records
            collector.summarize(output, protocol)
            continue
        reference = i9_wall(task['id'])
        limit = TOLERANCE * quiet_ratio(task['threads'])
        failures = 0
        while True:
            if not (runs / task['id'] / 'record.json').exists():
                wait_for_load(task['threads'], log, task['id'])
            record = collector.run_task(task, protocol, output, binaries)  # verifies/skips a kept record
            ratio = record['wall_seconds'] / reference if reference and record['status'] == 'ok' else None
            if ratio is None or ratio <= limit:
                break
            attempts = sorted(p for p in (rejected / task['id']).glob('attempt-*')) if (rejected / task['id']).exists() else []
            dest = rejected / task['id'] / f'attempt-{len(attempts) + 1}'
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(runs / task['id']), str(dest))
            load = record['before']['load'][0]
            log(task=task['id'], decision='rejected', ratio=round(ratio, 3), limit=round(limit, 3),
                wall=record['wall_seconds'], i9_wall=reference, load1_before=load, moved_to=str(dest.relative_to(output)))
            print(f"{utc()} {task['id']} REJECTED ratio={ratio:.2f} > {limit:.2f} (load1 {load:.1f}); retrying", flush=True)
            # Termination guard: consistent slowness is not contention; accept the fastest attempt.
            recent = sorted(rejected.glob(f"{task['id']}/attempt-*"), key=lambda p: int(p.name.split('-')[1]))[-CONSISTENT_ATTEMPTS:]
            if len(recent) == CONSISTENT_ATTEMPTS:
                walls = {p: json.loads((p / 'record.json').read_text())['wall_seconds'] for p in recent}
                low, high = min(walls.values()), max(walls.values())
                if (high - low) / low <= CONSISTENT_SPREAD:
                    best = min(walls, key=walls.get)
                    shutil.move(str(best), str(runs / task['id']))
                    log(task=task['id'], decision='accepted_consistent', restored_from=str(best.relative_to(output)),
                        ratio=round(walls[best] / reference, 3), limit=round(limit, 3), spread=round((high - low) / low, 3))
                    print(f"{utc()} {task['id']} accepted fastest of {CONSISTENT_ATTEMPTS} consistent attempts", flush=True)
                    break
            time.sleep(BACKOFF_SECONDS[min(failures, len(BACKOFF_SECONDS) - 1)])
            failures += 1
        collector.summarize(output, protocol)


if __name__ == '__main__':
    main()
