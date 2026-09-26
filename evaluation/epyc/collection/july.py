#!/usr/bin/env python3
"""July i9 lineardesign-clean benchmarks, re-run on EPYC 7313 NUMA node 1.

  july.py prepare   build the committed versions and write the synthetic fixtures (no timing)
  july.py plan      print the task list with the number of runs
  july.py run       run every task not yet recorded (resumable; one solve at a time)

Differences from the i9 originals (see AGENTS.md):
- The i9 synthetic proteins were never archived. The fixtures here are new: prefixes of one
  seeded uniform 20-amino-acid sequence (Met first). MFEs will not match the i9 values, so
  comparisons are between versions on this host, not same-input comparisons with the i9.
- Only committed states are built. The uncommitted i9 prototypes (ext collapse, A1-only, A1+A2,
  intermediate Round-3 states) are not reproduced.
- Node 1 has 16 cores, so i9 -j24 points are dropped; -j16 is the widest run.
"""
import argparse
import ctypes
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'research/results/epyc-node1-2026-09/july'
CPU_BASE, NODE = 16, 1
VERSIONS = {  # version -> commit whose lineardesign-clean/ tree is built
    'v0032': '787d3d7', 'v0033': '50066da', 'v0035': 'f0e2d09', 'v0037': '391fe6b', 'v0039': '72e2a1a'}
LENGTHS = (300, 1000, 2000, 3677)
SEED = 20260923
AMINO = 'ACDEFGHIKLMNPQRSTVWY'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def tasks():
    out = []
    def add(group, version, aa, lam, threads, reps=1):
        for r in range(1, reps + 1):
            out.append(dict(id=f'{group}-{version}-aa{aa}-l{lam}-j{threads}-r{r}', group=group,
                            version=version, aa=aa, lam=lam, threads=threads, repeat=r))
    # [0033] Round 3: 300 aa, lambda 1, best-of-3, baseline (v0032) vs final (v0033).
    for v in ('v0032', 'v0033'):
        for j in (1, 8, 16):
            add('r3', v, 300, 1, j, reps=3)
    # [0033] 8 concurrent -j1 processes (bandwidth test), v0033, 300 aa.
    out.append(dict(id='r3-concurrent8-v0033-aa300-l1-j1', group='concurrent8', version='v0033',
                    aa=300, lam=1, threads=1, repeat=1))
    # [0033] scale validation, lambda 1.
    for v in ('v0032', 'v0033'):
        for j in (8, 16):
            add('scale', v, 1000, 1, j)
        add('scale', v, 2000, 1, 16)
    add('scale', 'v0033', 1000, 1, 1)
    add('scale', 'v0033', 2000, 1, 1)
    add('scale', 'v0033', 2000, 1, 8)
    # [0039] i9 ladder: 3,677 aa, lambda 0.
    for j in (16, 8, 4, 2, 1):
        add('ladder39', 'v0039', 3677, 0, j)
    # [0035] A-stack before/after at lambda 0 (i9 used -j24 for 1000/2000 aa).
    for v in ('v0033', 'v0035'):
        add('stack35', v, 1000, 0, 16)
        add('stack35', v, 2000, 0, 16)
        add('stack35', v, 3677, 0, 16)
        add('stack35', v, 3677, 0, 8)
    # [0037]/[0038] -j1 cross-host: v0035 vs v0037, lambda 0.
    for v in ('v0035', 'v0037'):
        add('serial37', v, 1000, 0, 1)
    add('serial37', 'v0037', 3677, 0, 1)
    add('serial37', 'v0035', 3677, 0, 1)  # also the [0035] -j1 control
    # [0034] scale runs on v0033, 3,677 aa, lambda 1 (its 2000-aa rows reuse the scale group).
    for j in (16, 8, 1):
        add('scale34', 'v0033', 3677, 1, j)
    return out


def fixture(aa):
    return OUT / 'inputs' / f'synthetic-{aa}.fa'


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    seq = 'M' + ''.join(rng.choice(AMINO) for _ in range(max(LENGTHS) - 1))
    (OUT / 'inputs').mkdir(exist_ok=True)
    for aa in LENGTHS:
        fixture(aa).write_text(f'>synthetic-seed{SEED}-{aa}aa\n{seq[:aa]}\n')
    builds = {}
    for version, commit in VERSIONS.items():
        dest = OUT / 'builds' / version
        if not (dest / 'lineardesign-clean').exists():
            dest.mkdir(parents=True, exist_ok=True)
            archive = subprocess.run(['git', '-C', str(REPO), 'archive', commit, 'lineardesign-clean'],
                                     check=True, capture_output=True).stdout
            subprocess.run(['tar', '-x', '-C', str(dest)], input=archive, check=True)
            with (dest / 'build.log').open('w') as log:
                subprocess.run(['make', '-C', str(dest / 'lineardesign-clean'), 'lineardesign-clean'],
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        builds[version] = dict(commit=subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', commit], text=True).strip(),
                               binary_sha256=sha(dest / 'lineardesign-clean/lineardesign-clean'))
    manifest = dict(created_utc=utc(), seed=SEED, generator='random.Random(seed); "M" + uniform choice over ACDEFGHIKLMNPQRSTVWY; prefixes',
                    fixtures={aa: sha(fixture(aa)) for aa in LENGTHS}, builds=builds,
                    compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                    cpu_base=CPU_BASE, memory_node=NODE, tasks=tasks())
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Prepared', OUT, len(manifest['tasks']), 'tasks')


def bind_node():
    mask = ctypes.c_ulong(1 << NODE)  # MPOL_BIND=2, SYS_set_mempolicy=238 (x86_64)
    if ctypes.CDLL(None, use_errno=True).syscall(238, 2, ctypes.byref(mask), 64) != 0:
        os._exit(126)


def launch(task, workdir, cpu_list, tag):
    build = OUT / 'builds' / task['version'] / 'lineardesign-clean'
    cmd = ['/usr/bin/time', '-v', '-o', str(workdir / f'time{tag}.txt'), 'taskset', '-c', cpu_list,
           './lineardesign-clean', '-l', str(task['lam']), '-j', str(task['threads']), '-v']
    stdin = fixture(task['aa']).open('rb')
    out = (workdir / f'stdout{tag}.txt').open('wb'); err = (workdir / f'stderr{tag}.txt').open('wb')
    return cmd, subprocess.Popen(cmd, cwd=build, stdin=stdin, stdout=out, stderr=err, preexec_fn=bind_node)


def run_one(task, manifest):
    workdir = OUT / 'runs' / task['id']
    record_path = workdir / 'record.json'
    if record_path.exists():
        return
    workdir.mkdir(parents=True, exist_ok=True)
    binary = OUT / 'builds' / task['version'] / 'lineardesign-clean/lineardesign-clean'
    if sha(binary) != manifest['builds'][task['version']]['binary_sha256'] or sha(fixture(task['aa'])) != manifest['fixtures'][str(task['aa'])]:
        raise SystemExit(f'binary or fixture drift: {task["id"]}')
    with (OUT / '../compute.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        load, started = os.getloadavg(), utc()
        t0 = time.perf_counter()
        if task['group'] == 'concurrent8':
            procs = [launch(task, workdir, str(CPU_BASE + i), f'-{i}') for i in range(8)]
            codes = [p.wait() for _, p in procs]; cmds = [c for c, _ in procs]
        else:
            cpus = ','.join(str(CPU_BASE + i) for i in range(task['threads']))
            cmd, p = launch(task, workdir, cpus, '')
            codes, cmds = [p.wait()], [cmd]
        wall = time.perf_counter() - t0
        fcntl.flock(lock, fcntl.LOCK_UN)
    record = dict(task=task, status='ok' if all(c == 0 for c in codes) else 'failed', exit_codes=codes,
                  started_utc=started, ended_utc=utc(), wall_seconds=wall, load_before=load,
                  load_after=os.getloadavg(), commands=cmds,
                  binary_sha256=sha(binary), output_sha256={p.name: sha(p) for p in sorted(workdir.iterdir())})
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    print(f"{utc()} {task['id']} {record['status']} wall={wall:.2f}s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('action', choices=('prepare', 'plan', 'run'))
    args = ap.parse_args()
    if args.action == 'prepare':
        return prepare()
    manifest = json.loads((OUT / 'manifest.json').read_text())
    if manifest['tasks'] != tasks():
        raise SystemExit('Task list differs from the prepared manifest.')
    if args.action == 'plan':
        for t in manifest['tasks']:
            print(t['id'])
        print(len(manifest['tasks']), 'tasks')
        return
    if not set(range(CPU_BASE, CPU_BASE + 16)) <= os.sched_getaffinity(0):
        raise SystemExit('CPUs 16-31 are not all available')
    for task in manifest['tasks']:
        run_one(task, manifest)


if __name__ == '__main__':
    main()
