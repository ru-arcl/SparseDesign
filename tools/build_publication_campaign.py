#!/usr/bin/env python3
"""Build frozen publication kernels into a fresh, separate reproduction directory.

Never modifies recorded measurements. Requires Linux, GCC, GNU Make's usual
build prerequisites, GNU time, and util-linux taskset/flock for later runs.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template',type=Path,required=True,help='frozen evaluation/epyc/timings artifact')
    parser.add_argument('--output',type=Path,required=True,help='new empty reproduction directory')
    parser.add_argument('--compiler',default='g++')
    parser.add_argument('--lock',type=Path,help='shared lock used by other concurrent evaluations')
    args=parser.parse_args()
    template=args.template.resolve(); output=args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise SystemExit('Output must be empty; recorded results are never overwritten.')
    original=json.loads((template/'build.json').read_text())
    for name,digest in original['source_sha256'].items():
        if sha(template/'source'/name)!=digest:
            raise SystemExit(f'Frozen source checksum mismatch: {name}')
    protocol=json.loads((template/'protocol.json').read_text())
    if sha(template/'inputs/panel.json')!=protocol['panel_sha256']:
        raise SystemExit('Frozen input-panel checksum mismatch.')
    for task in protocol['tasks']:
        if sha(template/task['fasta'])!=task['fasta_sha256']:
            raise SystemExit(f'Frozen FASTA checksum mismatch: {task["id"]}')
    output.mkdir(parents=True,exist_ok=True)
    shutil.copytree(template/'source',output/'source')
    shutil.copytree(template/'inputs',output/'inputs')
    shutil.copyfile(template/'protocol.json',output/'protocol.json')
    if (template/'preflight-protocol.json').exists():
        shutil.copyfile(template/'preflight-protocol.json',output/'preflight-protocol.json')
    (output/'bin').mkdir()
    source=output/'source'; lock_path=args.lock or output/'compute.lock'
    lock_path.parent.mkdir(parents=True,exist_ok=True)
    commands=[]
    with lock_path.open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        for name,main,square in [('publication-bench',source/'publication_bench.cc',True),
                                 ('sparsedesign-square',source/'src/main.cc',True),
                                 ('sparsedesign-packed',source/'src/main.cc',False)]:
            command=[args.compiler,'-std=c++11','-O3','-flto','-Wall','-fopenmp']
            if square:command+=['-DLDCLEAN_SQUARE_MEMOS']
            command+=['-I',str(source/'src'),str(main)]
            command += [str(source/'src'/name) for name in ('fold_turner.cc','energy.cc','dfa.cc','codon_table.cc')]
            command+=['-o',str(output/'bin'/name)]
            commands.append(command)
            with (output/f'build-{name}.log').open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            print('Built',name,flush=True)
    build=dict(commands=commands,
               source_sha256={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()},
               binary_sha256={p.name:sha(p) for p in (output/'bin').iterdir()},
               compiler=subprocess.check_output([args.compiler,'--version'],text=True),
               uname=platform.uname()._asdict(),lscpu=subprocess.check_output(['lscpu'],text=True),
               meminfo=Path('/proc/meminfo').read_text(),
               execution_wrappers={'pilot':'source/publication_campaign_final.py','final_timings':'source/publication_campaign_final.py'},
               reproduction=dict(template_protocol_sha256=sha(template/'protocol.json'),
                                 template_build_sha256=sha(template/'build.json'),
                                 created_utc=datetime.now(timezone.utc).isoformat(),
                                 builder_sha256=sha(Path(__file__))))
    (output/'build.json').write_text(json.dumps(build,indent=2)+'\n')
    print('Fresh reproduction directory:',output)
    print('Run its frozen source/publication_campaign_final.py with --output pointing to this directory.')


if __name__=='__main__':
    main()
