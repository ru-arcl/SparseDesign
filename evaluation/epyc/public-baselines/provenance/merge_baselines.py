#!/usr/bin/env python3
"""Merge the 8 parallel native-comparison chains into canonical feasibility/ and repeated/ folders.

2026-09-25 (user): the native comparison ran as 8 independent chains (4 tools x 2 lambdas), each a
complete runner invocation in public-baselines/chains/<tool>-l<lambda>/{feasibility,repeated}.
The analyzer identifies phases by folder name, so each phase is merged into
public-baselines/<phase>/:

- record folders are hard-linked unchanged (labels are unique per tool/lambda/accession/repetition);
- planned-order.json is the concatenation of the chain manifests, in chain order;
- metadata.json is the first chain's metadata with tools, lambdas, pins, cpu and shared_lock replaced
  by their union over chains. Every other field must be identical across chains, and the per-chain
  metadata is kept under 'chains'.
"""
import argparse
import json
import os
from pathlib import Path

TOOLS = ['lineardesign', 'linearcdsfold', 'derna', 'sparsedesign']
LAMBDAS = [0.0, 4.0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True, help='public-baselines directory')
    args = ap.parse_args()
    root = args.root.resolve()
    chains = [root / 'chains' / f'{t}-l{int(l)}' for t in TOOLS for l in LAMBDAS]
    for phase in ('feasibility', 'repeated'):
        out = root / phase
        if out.exists():
            raise SystemExit(f'{out} exists; not overwriting')
        metas, planned = [], []
        for chain in chains:
            src = chain / phase
            meta = json.loads((src / 'metadata.json').read_text())
            plan = json.loads((src / 'planned-order.json').read_text())
            records = sorted(src.glob('*/record.json'))
            if len(records) != len(plan):
                raise SystemExit(f'{src}: {len(records)} records for {len(plan)} planned tasks')
            metas.append((chain.name, meta)); planned += plan
        base = dict(metas[0][1])
        varying = {'tools', 'lambdas', 'pins', 'cpu', 'shared_lock', 'panel_selection'}
        for name, meta in metas[1:]:
            for key in set(base) | set(meta):
                if key not in varying and key != 'feasibility_record_sha256' and base.get(key) != meta.get(key):
                    raise SystemExit(f'{name}: metadata differs in {key}')
        merged = dict(base, tools=TOOLS, lambdas=LAMBDAS,
                      pins={k: v for _, m in metas for k, v in m['pins'].items()},
                      cpu='per chain (see chains)', shared_lock='per chain (see chains)',
                      chains={name: m for name, m in metas})
        if phase == 'repeated':
            merged['feasibility_record_sha256'] = {k: v for _, m in metas for k, v in m.get('feasibility_record_sha256', {}).items()}
        out.mkdir()
        for chain in chains:
            src = chain / phase
            for record in sorted(src.glob('*/record.json')):
                dest = out / record.parent.name
                if dest.exists():
                    raise SystemExit(f'duplicate label {dest.name}')
                dest.mkdir()
                for f in record.parent.iterdir():
                    os.link(f, dest / f.name)
        (out / 'metadata.json').write_text(json.dumps(merged, indent=2, sort_keys=True) + '\n')
        (out / 'planned-order.json').write_text(json.dumps(planned, indent=2) + '\n')
        derna = chains[0] / phase / 'codon-table-derna.csv'
        if derna.exists():
            os.link(derna, out / 'codon-table-derna.csv')
        print(f'merged {phase}: {len(planned)} tasks from {len(chains)} chains')


if __name__ == '__main__':
    main()
