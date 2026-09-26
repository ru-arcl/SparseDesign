#!/usr/bin/env python3
"""Create the Q8VIM6 follow-up study (15 registered runs) for EPYC node 1.

Reuses the adapted timing study's binaries, source, collector and caps. The
tasks and their order come unchanged from the registered i9 follow-up protocol,
except that phase is relabelled 'ablation' so the frozen collector accepts them
(ids keep their 'quiet-' prefix). The i9 idle-host gate is not applied: it was
calibrated for a dedicated 24-CPU workstation, and this host is shared.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--timing-study', type=Path, required=True)
    ap.add_argument('--registered', type=Path, required=True, help='evaluation/quiet-replication/protocol.json')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    study, out = args.timing_study.resolve(), args.output.resolve()
    if out.exists():
        raise SystemExit('Output exists; not overwriting.')
    registered = json.loads(args.registered.read_text())
    protocol = json.loads((study / 'protocol.json').read_text())
    out.mkdir(parents=True)
    for name in ('source', 'inputs', 'bin'):
        shutil.copytree(study / name, out / name)
    for name in ('build.json', 'preflight-protocol.json', 'epyc-node1-adaptation.json'):
        shutil.copy2(study / name, out / name)
    tasks = []
    for task in registered['tasks']:
        task = dict(task, phase='ablation')
        if sha(out / task['fasta']) != task['fasta_sha256']:
            raise SystemExit(f'input mismatch for {task["id"]}')
        tasks.append(task)
    protocol['tasks'] = tasks
    protocol['followup'] = {
        'registered_protocol_sha256': sha(args.registered),
        'note': 'Registered i9 Q8VIM6 lambda-4 serial follow-up order; phase relabelled ablation for the collector; '
                'i9 quiet gate not applied on the shared EPYC host.',
    }
    (out / 'protocol.json').write_text(json.dumps(protocol, indent=2, sort_keys=True) + '\n')
    print('Prepared', out, len(tasks), 'tasks')


if __name__ == '__main__':
    main()
