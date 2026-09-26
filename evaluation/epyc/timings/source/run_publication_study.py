#!/usr/bin/env python3
"""Resume the frozen publication study and stop on any analysis/integrity failure.

The optional --await-pilot coordinates with an already running pilot; it does not
launch duplicate workers. Only the owner lock controls this driver. Individual
solves use the common compute.lock shared with public baseline evaluations.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--await-pilot',action='store_true')
    parser.add_argument('--phases',nargs='+',choices=('pilot','ablation','profile','scaling','workstation'),
                        default=['ablation','profile','scaling','workstation'])
    args=parser.parse_args()
    output=args.output.resolve()
    scripts=Path(__file__).resolve().parent
    protocol=json.loads((output/'protocol.json').read_text())
    def status(state,**kwargs):
        data=dict(state=state,updated_utc=datetime.now(timezone.utc).isoformat(),**kwargs)
        (output/'driver-status.json').write_text(json.dumps(data,indent=2)+'\n')
        print(json.dumps(data),flush=True)
    with (output/'driver.lock').open('a+') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('Another study driver already owns this output.')
        if args.await_pilot:
            status('awaiting_existing_pilot')
            pilot=[t for t in protocol['tasks'] if t['phase']=='pilot']
            while any(not (output/'timings/runs'/t['id']/'record.json').exists() for t in pilot):
                time.sleep(10)
        for phase in args.phases:
            status('running',phase=phase)
            cmd=[sys.executable,str(scripts/'publication_campaign.py'),'run','--output',str(output),'--phase',phase]
            with (output/f'{phase}-progress.log').open('a') as log:
                rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
            if rc:
                status('failed',phase=phase,exit_code=rc)
                raise SystemExit(rc)
        status('verifying_analysis')
        with (output/'compute.lock').open('a+') as compute:
            fcntl.flock(compute,fcntl.LOCK_EX)
            cmd=[sys.executable,str(scripts/'analyze_publication_campaign.py'),'--input',str(output),'--require-complete']
            with (output/'analysis-verification.log').open('w') as log:
                rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
        if rc:
            status('analysis_failed',exit_code=rc)
            raise SystemExit(rc)
        status('complete')


if __name__=='__main__':
    main()
