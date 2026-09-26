#!/usr/bin/env python3
"""Run the remaining single-thread timing tasks in parallel on NUMA node 1, scaled to host load.

User request (2026-09-25): run the remaining runs in parallel when the CPU load is not high,
using at most 12 cores in total: at most 4 on node 0 and at most 8 on node 1. Lowered the same day
(user) to at most 8 cores: at most 2 on node 0 and at most 6 on node 1. Runs already in progress are
never stopped: a new scheduler adopts running workers, counts them against the caps and starts
nothing new until usage is below the caps.

- Slots: at most MAX_WORKERS single-thread solves at once, each pinned to its own CPU with its
  memory bound to that CPU's NUMA node (never cross-socket). CPUs are preferred in CPU_ORDER:
  one per L3 group on node 1, then one per L3 group on node 0, then second cores on node 1, so
  concurrent runs rarely share an L3. A CPU is used only if it was idle (< IDLE_CPU busy) in
  the last sampling window, i.e. not occupied by another user.
- Load control: every POLL seconds, other users' busy CPUs are measured from /proc/stat as
  host busy minus our running solves. A new run starts only while
  running < min(MAX_WORKERS, HOST_TARGET - others_busy). When others use more than
  HOST_TARGET - MAX_WORKERS CPUs, fewer slots are used; at HOST_TARGET or above, nothing new
  starts. Running solves are never interrupted.
- Each task runs once (no retries) in a separate worker process. The worker imports the study's
  unchanged collector and replaces only (a) CPU_BASE, so a j1 task is pinned to its slot CPU, and
  (b) the collector's module-level fcntl, so its global compute.lock no longer serialises runs.
  Record contents and collector/protocol hashes are otherwise identical to a serial run.
- Only single-thread tasks are run here; multi-thread tasks stay with run_by_threads.py.

  run_parallel.py --study DIR [--study DIR2 ...]   schedule every remaining j1 task, in order
  run_parallel.py --worker --study DIR --cpu N --task ID   (internal)
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import types

# Preference: one core per L3 group (node 1, then node 0), then second cores, and so on; every CPU
# is a candidate so a busy core can be replaced by an idle one on the same node.
CPU_ORDER = tuple(c for k in range(4) for base in (16, 0) for c in (base + k, base + 4 + k, base + 8 + k, base + 12 + k))
MAX_WORKERS = 8  # user cap: at most 8 cores in total (was 12)
NODE_CAPS = {0: 2, 1: 6}  # user cap per NUMA node (CPUs 0-15 = node 0, 16-31 = node 1); was {0: 4, 1: 8}
IDLE_CPU = 0.25
HOST_CPUS, HOST_TARGET, POLL, WINDOW = 32, 28, 5, 5
BUSY_FIELDS = (0, 1, 2, 5, 6, 7)


def utc():
    return datetime.now(timezone.utc).isoformat()


def load_collector(study):
    path = study / 'source/publication_campaign_final.py'
    spec = importlib.util.spec_from_file_location('publication_campaign_final', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def binaries(study):
    return {'bench': study / 'bin/publication-bench', 'square': study / 'bin/sparsedesign-square',
            'packed': study / 'bin/sparsedesign-packed'}


def worker(study, cpu, task_id):
    collector = load_collector(study)
    collector.CPU_BASE = cpu  # j1 task -> cpus = [cpu]
    collector.MEMORY_NODE = cpu // 16  # bind memory to the slot CPU's own node
    collector.fcntl = types.SimpleNamespace(flock=lambda *a: None, LOCK_EX=0, LOCK_UN=0)
    protocol = json.loads((study / 'protocol.json').read_text())
    task = next(t for t in protocol['tasks'] if t['id'] == task_id)
    if task['threads'] != 1:
        raise SystemExit('parallel worker runs single-thread tasks only')
    collector.run_task(task, protocol, study, binaries(study))


def sample(running):
    """Return (other users' busy CPUs host-wide, per-CPU busy fraction) over one window."""
    def snap():
        out = {}
        for line in open('/proc/stat'):
            if line.startswith('cpu') and line[3] != ' ':
                f = line.split(); v = [int(x) for x in f[1:]]
                out[int(f[0][3:])] = (sum(v[i] for i in BUSY_FIELDS), sum(v))
        return out
    a = snap(); time.sleep(WINDOW); b = snap()
    per = {c: (b[c][0] - a[c][0]) / max(1, b[c][1] - a[c][1]) for c in a}
    return max(0.0, sum(per.values()) - running), per


class Adopted:
    """A worker started by an earlier scheduler; finished when its process is gone."""
    def __init__(self, pid):
        self.pid, self.returncode = pid, None

    def poll(self):
        if Path(f'/proc/{self.pid}').exists():
            return None
        self.returncode = 0
        return 0


def adopt_running():
    found = {}
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            argv = (proc / 'cmdline').read_bytes().split(b'\0')
        except OSError:
            continue
        argv = [a.decode() for a in argv if a]
        if '--worker' in argv and any(a.endswith('run_parallel.py') for a in argv):
            opt = {argv[i]: argv[i + 1] for i in range(len(argv) - 1) if argv[i].startswith('--')}
            found[int(opt['--cpu'])] = (Adopted(int(proc.name)), Path(opt['--study']), opt['--task'])
    return found


def schedule(studies):
    running = adopt_running()  # cpu -> (process, study, task_id)
    busy_tasks = {tid for _, _, tid in running.values()}
    if running:
        print(f'{utc()} adopted {len(running)} running workers on CPUs {sorted(running)}', flush=True)
    queue = []
    for study in studies:
        protocol = json.loads((study / 'protocol.json').read_text())
        for t in protocol['tasks']:
            if (t['phase'] != 'pilot' and t['threads'] == 1 and t['id'] not in busy_tasks
                    and not (study / 'timings/runs' / t['id'] / 'record.json').exists()):
                queue.append((study, t['id']))
    print(f'{utc()} {len(queue)} single-thread tasks queued', flush=True)
    last_note = None
    while queue or running:
        for cpu, (proc, study, tid) in list(running.items()):
            if proc.poll() is not None:
                del running[cpu]
                if proc.returncode:
                    print(f'{utc()} worker for {tid} exited {proc.returncode}', flush=True)
        busy, per_cpu = sample(len(running))
        allowed = max(0, min(MAX_WORKERS, int(HOST_TARGET - busy)))
        note = (allowed, len(running))
        if note != last_note:
            print(f'{utc()} others busy {busy:.1f} CPUs -> slots {allowed}; running {len(running)}; queued {len(queue)}', flush=True)
            last_note = note
        while queue and len(running) < allowed:
            per_node = {n: sum(1 for c in running if c // 16 == n) for n in NODE_CAPS}
            cpu = next((c for c in CPU_ORDER if c not in running and per_cpu.get(c, 1) < IDLE_CPU
                        and per_node[c // 16] < NODE_CAPS[c // 16]), None)
            if cpu is None:
                break
            study, tid = queue.pop(0)
            proc = subprocess.Popen([sys.executable, __file__, '--worker', '--study', str(study),
                                     '--cpu', str(cpu), '--task', tid])
            running[cpu] = (proc, study, tid)
        time.sleep(POLL)
    for study in studies:
        collector = load_collector(study)
        collector.summarize(study, json.loads((study / 'protocol.json').read_text()))
    print(f'{utc()} all queued single-thread tasks finished', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--study', type=Path, action='append', required=True)
    ap.add_argument('--worker', action='store_true')
    ap.add_argument('--cpu', type=int)
    ap.add_argument('--task')
    args = ap.parse_args()
    studies = [s.resolve() for s in args.study]
    if args.worker:
        return worker(studies[0], args.cpu, args.task)
    if not set(CPU_ORDER) <= os.sched_getaffinity(0):
        raise SystemExit('slot CPUs not all available (do not start this under numabind)')
    schedule(studies)


if __name__ == '__main__':
    main()
