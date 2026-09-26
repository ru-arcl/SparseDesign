#!/usr/bin/env python3
"""Audit pinned upstream model constants without copying upstream source code."""
import argparse
import json
from pathlib import Path
import re

from publication_baselines import PINS, sha, write_json


def source_text(path):
    return re.sub(r'/\*.*?\*/|//[^\n]*', '', path.read_text(), flags=re.S)


def array(path, name):
    match = re.search(r'\b'+name+r'\s*(?:\[[^]]*\])+\s*=\s*(\{.*?\});', source_text(path), re.S)
    if not match:
        raise ValueError(f'Missing initializer: {path}/{name}')
    text = match[1].replace('VIE_INF','10000000').replace('INF','10000000')
    text = re.sub(r',\s*}', '}', text)
    return json.loads(text.replace('{','[').replace('}',']'))


def differences(first, second, indices=()):
    if isinstance(first, list):
        if not isinstance(second, list) or len(first) != len(second):
            raise ValueError(f'Array shape mismatch at {indices}')
        for i, (a,b) in enumerate(zip(first,second)):
            yield from differences(a,b,indices+(i,))
    elif first != second:
        yield dict(indices=indices,sparsedesign=first,derna=second)


def motifs(path, name):
    match = re.search(r'\b'+name+r'\s*(?:\[[^]]*\])?\s*=\s*(.*?);',source_text(path),re.S)
    if not match:
        raise ValueError(f'Missing special-loop motifs: {path}/{name}')
    return ' '.join(re.findall(r'"([^"]*)"',match[1])).split()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sources',type=Path,required=True)
    ap.add_argument('--sparsedesign-source',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    sparse = args.sparsedesign_source/'src/vienna'
    derna = args.sources/'derna/src'
    checks, paths = [], set()
    names = ['stack37','hairpin37','bulge37','internal_loop37','mismatchH37',
             'mismatchI37','mismatch1nI37','mismatch23I37','Triloop37','Tetraloop37','Hexaloop37']
    for name, filename in [(n,'energy_parameter.h') for n in names]+[(n+'_37',f) for n,f in [('int11','intl11.h'),('int21','intl21.h'),('int22','intl22.h')]]:
        first = sparse/filename
        second = derna/('default.cpp' if filename == 'energy_parameter.h' else 'params/'+filename)
        paths.update((first,second))
        mismatch = list(differences(array(first,name),array(second,name)))
        pair_axes = 2 if name.startswith(('stack','int11','int21','int22')) else 1 if name.startswith('mismatch') else 0
        relevant = [m for m in mismatch if all(1<=i<=6 for i in m['indices'][:pair_axes]) and all(1<=i<=4 for i in m['indices'][pair_axes:])] if pair_axes else mismatch
        checks.append(dict(name=name,all_entries_equal=not mismatch,all_mismatch_count=len(mismatch),
                           canonical_entries_equal=not relevant,canonical_mismatches=relevant))
    scalars = []
    for name in ['lxc37','ML_intern37','ML_closing37','ML_BASE37','MAX_NINIO','ninio37','TerminalAU37']:
        values = []
        for path in (sparse/'energy_parameter.h',derna/'default.cpp'):
            match = re.search(r'\b'+name+r'\s*=\s*([-+.0-9]+)\s*;',source_text(path))
            if not match:
                raise ValueError(f'Missing scalar: {path}/{name}')
            values.append(float(match[1]))
        scalars.append(dict(name=name,sparsedesign=values[0],derna=values[1],equal=values[0]==values[1]))
    motif_checks = []
    for prefix in ('Tri','Tetra','Hexa'):
        first = motifs(sparse/'energy_parameter.h',prefix+'loops')
        second = motifs(derna/'default.cpp',prefix+'loopSeq')
        motif_checks.append(dict(kind=prefix.lower(),count=len(first),same_identities_and_order=first==second))
    all_equal = all(c['canonical_entries_equal'] for c in checks) and all(s['equal'] for s in scalars) and all(m['same_identities_and_order'] for m in motif_checks)
    write_json(args.output,dict(schema='sparsedesign-native-model-audit-v1',
               source_sha256={str(p):sha(p) for p in sorted(paths)},script_sha256=sha(__file__),
               pins=PINS,comparison='Frozen SparseDesign versus pinned DERNA; canonical pairs 1..6, nucleotides1..4 where applicable',
               checks=checks,scalars=scalars,special_loop_motifs=motif_checks,
               all_canonical_constants_equal=all_equal,
               interpretation='Parameter equality is necessary but does not prove recurrence or traceback correctness; outputs are independently translated and rescored.'))
    print(json.dumps(dict(arrays=len(checks),scalars=len(scalars),special_motif_classes=len(motif_checks),canonical_equal=all_equal)))


if __name__ == '__main__':
    main()
