#!/usr/bin/env python3
"""Verify frozen timing records and generate summaries/figures without new solves."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from datetime import datetime
import fcntl
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import re
import statistics
import sys

MODES = ('dense_typed', 'dense_scalar', 'sparse')
RATIOS = (('dense_typed','sparse'), ('dense_scalar','sparse'), ('dense_typed','dense_scalar'))
ELAPSED = 'Elapsed (wall clock) time (h:mm:ss or m:ss)'
RSS = 'Maximum resident set size (kbytes)'
PHASE_SECONDS = ('construction_seconds','closed_seconds','multiloop_seconds','external_seconds')


class ValidationError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def read_json(path):
    def pairs(items):
        result = {}
        for k,v in items:
            require(k not in result, f'duplicate JSON key: {k}')
            result[k] = v
        return result
    def bad(value):
        raise ValidationError(f'nonfinite JSON constant: {value}')
    return json.loads(path.read_text(), object_pairs_hook=pairs, parse_constant=bad)


def contained(root, relative):
    path = Path(relative)
    require(not path.is_absolute() and '..' not in path.parts, f'unsafe path: {relative}')
    result = root/path
    require(result.resolve().is_relative_to(root.resolve()), f'path escapes archive: {relative}')
    return result


def number(value, label, minimum=None):
    require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value), f'{label}: not finite numeric')
    require(minimum is None or value >= minimum, f'{label}: below {minimum}')
    return value


def parse_kv(text, strict=False):
    result = {}
    for token in text.strip().split():
        if '=' not in token:
            require(not strict, 'unexpected harness token')
            continue
        k,v = token.split('=',1)
        require(not strict or k not in result, f'duplicate harness field: {k}')
        try:
            result[k] = float(v) if any(c in v for c in '.eE') else int(v)
        except ValueError:
            result[k] = v
    return result


def elapsed_seconds(text):
    require(isinstance(text,str) and re.fullmatch(r'\d+:\d{2}(?::\d{2})?(?:\.\d+)?',text), f'malformed GNU elapsed time: {text!r}')
    parts = [float(x) for x in text.split(':')]
    require(parts[-1]<60 and (len(parts)==2 or parts[-2]<60), f'GNU elapsed component out of range: {text}')
    return sum(x*60**i for i,x in enumerate(reversed(parts)))


def parse_time(text, strict=True):
    result = {}
    for line in text.splitlines():
        if ': ' in line:
            k,v = line.lstrip().split(': ',1)
            v=v.strip()
            require(k not in result, f'duplicate GNU-time field: {k}')
            result[k] = v
    if strict and ELAPSED in result:
        elapsed_seconds(result[ELAPSED])
    for k in (RSS,'Exit status','Swaps'):
        if strict and k in result:
            require(re.fullmatch(r'\d+',result[k]), f'invalid GNU integer: {k}')
    for k in ('User time (seconds)','System time (seconds)'):
        if strict and k in result:
            number(float(result[k]),k,0)
    return result


def expected_tasks(panel, dmd_hash):
    """Reconstruct the declared v1 scheduling algorithm without executing its runner."""
    require(len(panel)==12, 'version-1 protocol requires the frozen 12-case panel')
    tasks=[]
    def add(phase,item,lam,threads,mode,repeat,layout='square'):
        tasks.append({'id':f"{phase}-{item['accession']}-l{lam}-j{threads}-{mode}-{layout}-r{repeat}",
                      'phase':phase,'accession':item['accession'],'aa_length':item['aa_length'],
                      'fasta':item['fasta_path'],'fasta_sha256':item['fasta_sha256'],
                      'lambda':lam,'threads':threads,'mode':mode,'repeat':repeat,'layout':layout})
    for item,threads,mode in itertools.product([panel[i] for i in (0,2,6,10)],(1,16),MODES):
        add('pilot',item,0,threads,mode,0)
    rng=random.Random(20260906)
    conditions=list(itertools.product(panel,(0,4),(1,16)))
    rng.shuffle(conditions)
    permutations=list(itertools.permutations(MODES))
    for item,lam,threads in conditions:
        order=rng.randrange(6)
        for repeat in range(1,6):
            for mode in permutations[(order+repeat-1)%6]:
                add('ablation',item,lam,threads,mode,repeat)
    for item,lam,threads,mode in itertools.product(panel,(0,4),(1,16),MODES):
        add('profile',item,lam,threads,mode,1)
    for item in (panel[6],panel[10]):
        conditions=list(itertools.product((0,4),(1,2,4,8,16),range(1,4)))
        rng.shuffle(conditions)
        for lam,threads,repeat in conditions:
            add('scaling',item,lam,threads,'cli',repeat)
    dmd=dict(accession='NP_000100.3',aa_length=3677,fasta_path='inputs/NP_000100.3.fasta',fasta_sha256=dmd_hash)
    for repeat in range(1,6):
        for layout in (('square','packed') if repeat%2 else ('packed','square')):
            add('workstation',dmd,0,16,'cli',repeat,layout)
    return tasks


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


GENERATOR_SHA256=sha(Path(__file__))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def quantile(xs, p):
    require(bool(xs) and 0<=p<=1, 'quantile requires observations and probability in [0,1]')
    xs = sorted(xs)
    position = (len(xs) - 1) * p
    lo = int(position)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (position - lo)


def describe(xs):
    return {'n': len(xs), 'median': statistics.median(xs), 'min': min(xs),
            'max': max(xs), 'q25': quantile(xs, .25), 'q75': quantile(xs, .75)} if xs else None


def gnu_seconds(record):
    value = record['gnu_time'].get(ELAPSED)
    try:
        return elapsed_seconds(value) if value is not None else None
    except ValidationError:
        require(record['status']!='ok','malformed successful-run GNU elapsed time')
        return None


def rss_gib(record):
    value = record['gnu_time'].get(RSS)
    if value is None:return None
    if not re.fullmatch(r'\d+',value):
        require(record['status']!='ok','malformed successful-run GNU RSS')
        return None
    return int(value)/2**20


def vm_delta(record, key):
    def lookup(text):
        rows=[line.split() for line in text.splitlines()]
        require(all(len(r)==2 for r in rows) and len(dict(rows))==len(rows),'malformed vmstat snapshot')
        return dict(rows)[key]
    delta=int(lookup(record['after']['vmstat'])) - int(lookup(record['before']['vmstat']))
    require(delta>=0,'host swap counter decreased')
    return delta


def verified_records(root, require_complete=False):
    root=Path(root)
    protocol = read_json(root/'protocol.json')
    build = read_json(root/'build.json')
    panel = read_json(root/'inputs/panel.json')['panel']
    require(protocol['schema']=='sparsedesign-timing-protocol-v1' and protocol['seed']==20260906,'unsupported protocol/seed')
    require(sha(root/'inputs/panel.json') == protocol['panel_sha256'], 'panel drift')
    require(len({p['accession'] for p in panel})==len(panel),'duplicate panel accession')
    for item in panel:
        path=contained(root,item['fasta_path'])
        require(sha(path)==item['fasta_sha256'],'panel FASTA drift')
        lines=path.read_text().splitlines()
        sequence=''.join(''.join(line.split()).upper() for line in lines if not line.startswith('>'))
        require(sum(line.startswith('>') for line in lines)==1,'FASTA must contain one record')
        require(len(sequence)==item['aa_length'] and hashlib.sha256(sequence.encode()).hexdigest()==item['sequence_sha256'],'panel sequence/length mismatch')
    dmd=root/'inputs/NP_000100.3.fasta'
    lines=dmd.read_text().splitlines()
    require(sum(line.startswith('>') for line in lines)==1 and len(''.join(''.join(line.split()) for line in lines if not line.startswith('>')))==3677,'dystrophin input record/length mismatch')
    require(json.dumps(protocol['tasks'],sort_keys=True)==json.dumps(expected_tasks(panel,sha(dmd)),sort_keys=True),'planned task grid/order mismatch')
    require(protocol['openmp']==dict(OMP_DYNAMIC='false',OMP_PROC_BIND='close',OMP_PLACES='cores',OMP_WAIT_POLICY='active',OMP_THREAD_LIMIT='16'),'OpenMP protocol changed')
    # EPYC node-1 adaptation: every wall cap is multiplied by the recorded factor (see epyc-node1-adaptation.json).
    F=read_json(root/'epyc-node1-adaptation.json')['cap_factor']
    require(protocol['address_space_limit_gib']==32 and protocol['workstation_timeout_seconds']==1800*F,'memory/workstation limits changed')
    if 'preflight_protocol_sha256' in protocol:
        require(sha(root/'preflight-protocol.json')==protocol['preflight_protocol_sha256'],'preflight protocol drift')
        preflight=read_json(root/'preflight-protocol.json')
        require(all(protocol[k]==v for k,v in preflight.items() if k not in ('timeout_seconds','workstation_timeout_seconds','affinity')),'undeclared change from preflight protocol')  # EPYC: caps x F and node-1 affinity declared in host_adaptation
        require(bool(protocol.get('host_adaptation')),'missing host adaptation record')
        require(preflight['timeout_seconds']==600 and protocol['timeout_seconds']==900*F and protocol['phase_timeout_seconds']==dict(pilot=600*F,ablation=900*F,profile=900*F,scaling=900*F,workstation=1800*F),'timeout amendment differs from declared caps')
        require(bool(protocol['pilot_amendment']),'missing amendment explanation')
    else:
        require(protocol['timeout_seconds']==600,'unrecorded timeout amendment')
    require({'publication_bench.cc','publication_campaign.py','src/fold_turner.cc','src/main.cc','data/codon_usage_freq_table_human.csv'}<=set(build['source_sha256']),'incomplete source identities')
    for name, digest in build['source_sha256'].items():
        require(sha(contained(root/'source',name)) == digest, f'source drift: {name}')
    require(set(build['binary_sha256'])=={'publication-bench','sparsedesign-square','sparsedesign-packed'},'binary identities incomplete')
    for name,digest in build['binary_sha256'].items():
        if (root/'bin'/name).exists():
            require(sha(root/'bin'/name)==digest,f'present binary drift: {name}')
    records, missing = [], []
    tasks = protocol['tasks']
    planned={t['id'] for t in tasks}
    require(len(planned)==len(tasks),'duplicate task identity')
    # Snapshot completed record paths once; atomically published records have stable outputs.
    paths={p.parent.name:p for p in (root/'timings/runs').glob('*/record.json')}
    require(set(paths)<=planned,'unplanned completed record')
    checked=set()
    for task in tasks:
        if task['fasta'] not in checked:
            require(sha(contained(root,task['fasta'])) == task['fasta_sha256'], 'FASTA drift')
            checked.add(task['fasta'])
        dest = root/'timings/runs'/task['id']
        if task['id'] not in paths:
            missing.append(task['id'])
            continue
        r = read_json(dest/'record.json')
        require(json.dumps(r['task'],sort_keys=True)==json.dumps(task,sort_keys=True), f'task drift: {task["id"]}')
        require(r['status'] in {'ok','failed','timeout','invalid_output'} and isinstance(r['exit_code'],int),'unknown status/exit')
        binary = 'publication-bench' if task['mode'] != 'cli' else f'sparsedesign-{task["layout"]}'
        require(r['binary_sha256'] == build['binary_sha256'][binary], 'binary drift')
        require(r['codon_table_sha256'] == build['source_sha256']['data/codon_usage_freq_table_human.csv'], 'table drift')
        wrappers=build.get('execution_wrappers',{})
        wrapper=wrappers.get('pilot' if task['phase']=='pilot' else 'final_timings')
        if wrapper:
            require(wrapper.startswith('source/'),'collector snapshot must be within frozen source')
            require(wrapper[len('source/'):] in build['source_sha256'],'collector missing from source identity manifest')
            if Path(wrapper).name=='publication_campaign_final.py' or 'collector_sha256' in r:
                require(r.get('collector_sha256')==sha(contained(root,wrapper)),'per-run collector identity mismatch')
                require(r.get('protocol_sha256')==sha(root/'protocol.json'),'per-run protocol identity mismatch')
        files={'stdout.txt','stderr.txt'}|({'time.txt'} if (dest/'time.txt').exists() else set())
        require(set(r['output_sha256'])==files,'required output hash coverage mismatch')
        for name, digest in r['output_sha256'].items():
            require(Path(name).name == name and sha(dest/name) == digest, 'output drift')
        raw_time=(dest/'time.txt').read_text() if 'time.txt' in files else ''
        require(parse_time(raw_time,r['status']=='ok')==r['gnu_time'],'GNU-time dictionary differs from hashed transcript')
        stdout,stderr=(dest/'stdout.txt').read_text(),(dest/'stderr.txt').read_text()
        metrics={}
        if task['mode']=='cli':
            for line in stderr.splitlines():
                if 'dp_cost=' in line:
                    metrics.update(parse_kv(line))
            if 'dp_cost' in metrics:
                metrics['objective_units']=metrics['dp_cost']
        else:
            metrics=parse_kv(stdout)
        require(metrics==r['metrics'],'metrics differ from hashed output transcript')
        cmd=r['command']
        require(isinstance(cmd,list) and len(cmd)>=8 and all(isinstance(x,str) for x in cmd),'invalid command')
        require(cmd[:3]==['/usr/bin/time','-v','-o'] and cmd[4:6]==['taskset','-c'],'command wrapper mismatch')
        require(Path(cmd[3]).parent.name==task['id'] and Path(cmd[3]).name=='time.txt','time destination mismatch')
        require((cmd[6]==','.join(str(16+i) for i in range(task['threads'])) or (task['threads']==1 and cmd[6].isdigit() and 0<=int(cmd[6])<=31)) and Path(cmd[7]).name==binary,'CPU affinity or binary mismatch')  # EPYC: parallel j1 slots may use any CPU (node-local memory)
        tail=cmd[8:]
        if task['mode']=='cli':
            require(len(tail)==7 and tail[:5]==['-l',str(task['lambda']),'-j',str(task['threads']),'-c'] and Path(tail[5]).name=='codon_usage_freq_table_human.csv' and tail[6]=='--diagnostics','CLI arguments mismatch')
        else:
            require(len(tail)==5+int(task['phase']=='profile') and tail[:3]==[task['mode'],str(task['lambda']),str(task['threads'])] and Path(tail[3]).name=='codon_usage_freq_table_human.csv' and Path(tail[4]).name==Path(task['fasta']).name,'harness arguments mismatch')
            require(task['phase']!='profile' or tail[-1]=='profile','missing profile argument')
        timed=r['gnu_time'].get('Command being timed')
        if timed is not None:
            # GNU time displays arguments joined by spaces, without preserving
            # shell quoting. Compare its actual rendering; recorded argv stays typed.
            require(timed=='"'+' '.join(cmd[4:])+'"','GNU-time command mismatch')
        number(r['wall_seconds'],'observer wall',0)
        start,end=(datetime.fromisoformat(r[k]) for k in ('started_utc','ended_utc'))
        require(start.tzinfo is not None and end.tzinfo is not None and (end-start).total_seconds()+.03>=r['wall_seconds'],'event/duration mismatch')
        if 'finalized_utc' in protocol and task['phase']!='pilot':
            require(start>=datetime.fromisoformat(protocol['finalized_utc']),'nonpilot measurement precedes protocol finalization')
        for k in ('pswpin','pswpout'):
            vm_delta(r,k)
        if r['status'] == 'ok':
            number(metrics.get('objective_units'),'objective')
            require(r['exit_code']==0 and r['gnu_time'].get('Exit status')=='0','successful record has nonzero/missing exit')
            require(rss_gib(r) is not None and rss_gib(r)>0 and gnu_seconds(r) is not None,'missing process metrics')
            require(gnu_seconds(r)<=r['wall_seconds']+.03,'GNU wall exceeds observer duration')
            if task['mode'] != 'cli':
                require(not stderr and len(stdout.strip().splitlines())==1 and parse_kv(stdout,True)==metrics,'invalid successful harness transcript')
                keys={'mode','lambda','threads_requested','threads_effective','protein_aa','rna_nt','lattice_nodes','objective_units','solve_seconds','profiled'}
                require(set(metrics)==keys|(set(PHASE_SECONDS)|{'split_visits','split_eligible'} if task['phase']=='profile' else set()),'harness field set mismatch')
                for k,v in [('mode',task['mode']),('lambda',task['lambda']),('threads_requested',task['threads']),('threads_effective',task['threads']),('protein_aa',task['aa_length']),('rna_nt',3*task['aa_length']),('profiled',int(task['phase']=='profile'))]:
                    require(metrics[k]==v,f'harness {k} differs from task')
                number(metrics['solve_seconds'],'solve_seconds',0)
                require(0<metrics['solve_seconds']<=gnu_seconds(r)+.03,'solve duration outside process boundary')
                number(metrics['lattice_nodes'],'lattice nodes',3*task['aa_length']+1)
                if task['phase']=='profile':
                    for k in PHASE_SECONDS:
                        number(metrics[k],k,0)
                    require(all(isinstance(metrics[k],int) and metrics[k]>=0 for k in ('split_visits','split_eligible')) and metrics['split_eligible']<=metrics['split_visits'],'invalid split counters')
                    require(sum(metrics[k] for k in PHASE_SECONDS)<=metrics['solve_seconds']+1e-6,'phase timers exceed enclosing cost call')
            else:
                require(len(stderr.strip().splitlines())==1 and stderr.startswith('diagnostics '),'unexpected successful CLI stderr')
                for k in ('dp_cost','realized_cost','objective_delta','codon_penalty','cai'):
                    number(metrics.get(k),k)
                require(metrics.get('cost_unit')=='centikcal_mol','diagnostic unit mismatch')
                require(abs(metrics['objective_delta'])<=1e-4 and math.isclose(metrics['realized_cost']-metrics['dp_cost'],metrics['objective_delta'],abs_tol=1e-7),'reconstruction discrepancy')
                require(0<metrics['cai']<=1 and metrics['codon_penalty']>=0 and math.isclose(metrics['codon_penalty'],-task['aa_length']*math.log(metrics['cai']),abs_tol=1e-7),'CAI/codon-penalty discrepancy')
                seq=re.findall(r'^mRNA sequence:\s+([ACGU]+)$',stdout,re.M)
                struct=re.findall(r'^mRNA structure:\s+([.()]+)$',stdout,re.M)
                energy=re.findall(r'^mRNA folding free energy: (-?\d+\.\d{2}) kcal/mol; mRNA CAI: (\d+\.\d{3})$',stdout,re.M)
                require(len(seq)==len(struct)==len(energy)==1,'incomplete CLI output')
                require(len(seq[0])==len(struct[0])==3*task['aa_length'],'CLI length mismatch')
                stack=0
                for c in struct[0]:
                    stack+=(c=='(')-(c==')')
                    require(stack>=0,'unbalanced structure')
                require(stack==0,'unbalanced structure')
                require(abs(float(energy[0][0])*100+100*task['lambda']*metrics['codon_penalty']-metrics['realized_cost'])<=.501,'printed MFE/objective disagreement')
                require(abs(float(energy[0][1])-metrics['cai'])<=.000501,'printed CAI disagreement')
        elif r['status']=='timeout':
            cap=protocol.get('phase_timeout_seconds',{}).get(task['phase'],protocol['workstation_timeout_seconds'] if task['phase']=='workstation' else protocol['timeout_seconds'])
            require(r['exit_code']!=0 and r['wall_seconds']>=cap-1,'timeout duration/exit mismatch')
        elif r['status']=='failed':
            require(r['exit_code']!=0,'failed record with zero exit')
        else:
            require(r['exit_code']==0,'invalid-output record with nonzero exit')
        r['_record_sha256']=sha(dest/'record.json')
        records.append(r)
    # EPYC node 1 (2026-09-25): from PARALLEL_START, single-thread runs were deliberately run in
    # parallel (research/epyc-node1/run_parallel.py). Overlap is allowed only between two such runs.
    PARALLEL_START=datetime.fromisoformat('2026-09-25T12:45:00+00:00')
    parallel=lambda r:r['task']['threads']==1 and datetime.fromisoformat(r['started_utc'])>=PARALLEL_START
    intervals=sorted((datetime.fromisoformat(r['started_utc']),r['wall_seconds'],r['task']['id'],parallel(r)) for r in records)
    for i,a in enumerate(intervals):
        for b in intervals[i+1:]:
            if (b[0]-a[0]).total_seconds()+.03>=a[1]:
                break
            require(a[3] and b[3],f'overlapping measurement windows: {a[2]}, {b[2]}')
    require(not require_complete or not missing,f'{len(missing)} planned tasks missing')
    return protocol, build, records, missing


def aggregate_ablation(cases, draws):
    """Shared complete triplets; resample joint indices across all three ratios."""
    results=[]
    for lam,threads in sorted({(c['lambda'],c['threads']) for c in cases}):
        group=[c for c in cases if (c['lambda'],c['threads'])==(lam,threads)]
        full=[c for c in group if c['complete_triplets']==c['expected_repetitions']]
        seed=int(hashlib.sha256(f'20260908:{lam}:{threads}'.encode()).hexdigest()[:16],16)
        rng=random.Random(seed); boots=defaultdict(list)
        if full and draws:
            for _ in range(draws):
                indices=[[rng.randrange(c['expected_repetitions']) for _ in range(c['expected_repetitions'])] for c in full]
                for num,den in RATIOS:
                    name=num+'_over_'+den
                    medians=[statistics.median(c['paired_ratios'][name][i] for i in chosen) for c,chosen in zip(full,indices)]
                    boots[name].append(math.exp(statistics.mean(math.log(v) for v in medians)))
        for num,den in RATIOS:
            name=num+'_over_'+den
            medians=[statistics.median(c['paired_ratios'][name]) for c in full]
            results.append({'lambda':lam,'threads':threads,'ratio':name,'complete_cases':len(full),'expected_cases':len(group),
                            'complete_panel':len(full)==len(group),'included_accessions':[c['accession'] for c in full],
                            'excluded_accessions':[c['accession'] for c in group if c not in full],
                            'geomean_of_paired_medians':math.exp(statistics.mean(math.log(v) for v in medians)) if medians else None,
                            'conditional_repeat_ci95':[quantile(boots[name],.025),quantile(boots[name],.975)] if boots[name] else None,
                            'case_medians':describe(medians),'bootstrap_seed':seed,
                            'estimate_scope':'fixed planned panel' if len(full)==len(group) else 'available fully observed subset; not a full-panel estimate'})
    return results


def analyze(root, require_complete=False, draws=10000):
    require(isinstance(draws,int) and draws>=0,'bootstrap draws must be a nonnegative integer')
    require(sha(Path(__file__))==GENERATOR_SHA256,'analyzer source changed after process startup; rerun with the frozen source')
    root=Path(root)
    protocol,build,records,missing=verified_records(root,require_complete)
    tasks=protocol['tasks']; by_id={r['task']['id']:r for r in records}
    objectives=defaultdict(list);nodes=defaultdict(set); outcomes=[]
    for t in tasks:
        r=by_id.get(t['id'])
        outcomes.append({**{k:t[k] for k in ('id','phase','accession','lambda','threads','mode','layout','repeat')},
                         'status':r['status'] if r else 'missing','exit_code':r['exit_code'] if r else None,
                         'wall_observer_seconds':r['wall_seconds'] if r else None,'wall_gnu_seconds':gnu_seconds(r) if r else None,
                         'rss_gib':rss_gib(r) if r else None,'record_sha256':r['_record_sha256'] if r else None})
        if r and r['status']=='ok':
            objectives[(t['accession'],t['lambda'])].append(r['metrics']['objective_units'])
            if t['mode']!='cli':nodes[t['accession']].add(r['metrics']['lattice_nodes'])
    require(all(len(v)==1 for v in nodes.values()),'same-input lattice size disagreement')
    deltas={str(k):max(v)-min(v) for k,v in objectives.items()}
    require(all(d<=1e-4 for d in deltas.values()),f'objective disagreement: {deltas}')
    completion={}
    for phase in dict.fromkeys(t['phase'] for t in tasks):
        counts=Counter(o['status'] for o in outcomes if o['phase']==phase); total=sum(counts.values())
        completion[phase]={'expected':total,'recorded':total-counts['missing'],**{k:counts[k] for k in ('ok','missing','failed','timeout','invalid_output')}}
    blocks=defaultdict(dict)
    for t in tasks:
        if t['phase']=='ablation':
            blocks[(t['accession'],t['lambda'],t['threads'],t['repeat'])][t['mode']]=by_id.get(t['id'])
    cases=[]
    for key in sorted({k[:3] for k in blocks}):
        reps=sorted(k[3] for k in blocks if k[:3]==key)
        complete=[rep for rep in reps if all(r and r['status']=='ok' for r in blocks[key+(rep,)].values())]
        times,rss,ratios,available_times,statuses={},{},{},{},{}
        for mode in MODES:
            paired=[blocks[key+(rep,)][mode] for rep in complete]
            available=[blocks[key+(rep,)][mode] for rep in reps if blocks[key+(rep,)][mode] and blocks[key+(rep,)][mode]['status']=='ok']
            times[mode]=describe([r['metrics']['solve_seconds'] for r in paired]); rss[mode]=describe([rss_gib(r) for r in paired])
            available_times[mode]=describe([r['metrics']['solve_seconds'] for r in available])
            statuses[mode]=dict(Counter(blocks[key+(rep,)][mode]['status'] if blocks[key+(rep,)][mode] else 'missing' for rep in reps))
        for num,den in RATIOS:
            ratios[num+'_over_'+den]=[blocks[key+(rep,)][num]['metrics']['solve_seconds']/blocks[key+(rep,)][den]['metrics']['solve_seconds'] for rep in complete]
        cases.append({'accession':key[0],'lambda':key[1],'threads':key[2],'expected_repetitions':len(reps),
                      'complete_triplets':len(complete),'paired_repetitions':complete,'outcomes':statuses,'times':times,'rss_gib':rss,
                      'all_successful_unpaired_times':available_times,'paired_ratios':ratios,'ratios':{k:describe(v) for k,v in ratios.items()}})
    profiles=[]; profile_cells=defaultdict(dict)
    for t in tasks:
        if t['phase']=='profile':
            r=by_id.get(t['id']); profiles.append({'task':t,'status':r['status'] if r else 'missing','metrics':r['metrics'] if r else None})
            profile_cells[(t['accession'],t['lambda'],t['threads'])][t['mode']]=r
    comparisons=[]
    for key,group in sorted(profile_cells.items()):
        if all(r and r['status']=='ok' for r in group.values()):
            sparse,scalar,typed=(group[m]['metrics'] for m in ('sparse','dense_scalar','dense_typed'))
            ratio=lambda x,y:x/y if y else None
            comparisons.append({'accession':key[0],'lambda':key[1],'threads':key[2],
                'scalar_over_sparse_split_visits':ratio(scalar['split_visits'],sparse['split_visits']),
                'scalar_over_sparse_eligible_visits':ratio(scalar['split_eligible'],sparse['split_eligible']),
                'scalar_over_sparse_multiloop_seconds':ratio(scalar['multiloop_seconds'],sparse['multiloop_seconds']),
                'typed_over_sparse_multiloop_seconds':ratio(typed['multiloop_seconds'],sparse['multiloop_seconds']),
                'multiloop_fraction_of_profiled_call':{m:group[m]['metrics']['multiloop_seconds']/group[m]['metrics']['solve_seconds'] for m in MODES},
                'residual_profiled_call_seconds':{m:group[m]['metrics']['solve_seconds']-sum(group[m]['metrics'][k] for k in PHASE_SECONDS) for m in MODES}})
    scaling=[];workstation=[]
    for phase,keys in [('scaling',('accession','lambda','threads')),('workstation',('layout',))]:
        planned=defaultdict(list)
        for t in tasks:
            if t['phase']==phase:planned[tuple(t[k] for k in keys)].append(t)
        for key,group in sorted(planned.items()):
            rs=[by_id[t['id']] for t in group if t['id'] in by_id]; ok=[r for r in rs if r['status']=='ok']
            counts=Counter(r['status'] for r in rs); counts['missing']=len(group)-len(rs)
            row={**dict(zip(keys,key)),'expected':len(group),'complete':len(rs)==len(group),'all_successful':len(ok)==len(group),
                 'outcomes':dict(counts),'wall_seconds':describe([gnu_seconds(r) for r in ok]),'rss_gib':describe([rss_gib(r) for r in ok]),
                 'successful_repeats':[r['task']['repeat'] for r in ok]}
            if phase=='scaling':scaling.append(row)
            else:
                row.update(host_swap_in_pages=[vm_delta(r,'pswpin') for r in rs],host_swap_out_pages=[vm_delta(r,'pswpout') for r in rs]);workstation.append(row)
    for row in scaling:
        base=next(r for r in scaling if (r['accession'],r['lambda'],r['threads'])==(row['accession'],row['lambda'],1))
        row['speedup_ratio_of_medians']=base['wall_seconds']['median']/row['wall_seconds']['median'] if base['all_successful'] and row['all_successful'] else None
        row['parallel_efficiency_ratio_of_medians']=row['speedup_ratio_of_medians']/row['threads'] if row['speedup_ratio_of_medians'] else None
    layout_pairs=[]
    for repeat in sorted({t['repeat'] for t in tasks if t['phase']=='workstation'}):
        pair={t['layout']:by_id.get(t['id']) for t in tasks if t['phase']=='workstation' and t['repeat']==repeat}
        if all(r and r['status']=='ok' for r in pair.values()):
            layout_pairs.append({'repeat':repeat,'packed_over_square_wall':gnu_seconds(pair['packed'])/gnu_seconds(pair['square']),
                                 'packed_over_square_peak_rss':rss_gib(pair['packed'])/rss_gib(pair['square'])})
    return {'schema':'sparsedesign-timing-analysis-v2','complete':not missing,'all_successful':not missing and all(r['status']=='ok' for r in records),
            'protocol_sha256':sha(root/'protocol.json'),'preflight_protocol_sha256':protocol.get('preflight_protocol_sha256'),
            'build_sha256':sha(root/'build.json'),'generator_sha256':GENERATOR_SHA256,
            'measurement_source_sha256':{k:build['source_sha256'][k] for k in ('publication_bench.cc','src/fold_turner.cc','src/main.cc')},
            'phase_completion':completion,'missing':missing,'maximum_objective_delta_units':max(deltas.values(),default=None),
            'objective_deltas_units':deltas,'ablation_cases':cases,'ablation_aggregate':aggregate_ablation(cases,draws),
            'profiles':profiles,'profile_comparisons':comparisons,'scaling':scaling,'workstation':workstation,'outcomes':outcomes,
            'workstation_paired':{'complete_pairs':len(layout_pairs),'expected_pairs':len({t['repeat'] for t in tasks if t['phase']=='workstation'}),
                                 'observations':layout_pairs,'packed_over_square_wall':describe([p['packed_over_square_wall'] for p in layout_pairs]),
                                 'packed_over_square_peak_rss':describe([p['packed_over_square_peak_rss'] for p in layout_pairs])},
            'bootstrap':{'draws':draws,'seed':20260908,'estimand':'geometric mean across proteins of medians of paired solve-time ratios; only fully observed three-arm cases',
                'resampling':'joint repetition indices across all three ratios, independently within each protein/lambda/thread case',
                'interval':'linear-interpolated 2.5/97.5 percentile interval; absent when draws=0 or no complete cases',
                'limitation':'repeat variability conditional on explicitly listed cases; incomplete subsets are not full-panel estimates; no population generalization or correction for systematic drift'},
            'measurement_definitions':{
                'cost_seconds':'Harness steady-clock interval around cost call: solver construction, fill and destruction; excludes table/DFA construction, effective-team probe and traceback.',
                'cli_wall_seconds':'GNU-time complete child process: startup, DFA, design, traceback, output and destruction. Python Popen/wait wall is retained separately.',
                'profile_seconds':'Separate instrumented runs. Construction is solver-constructor elapsed time. Closed, multiloop and external sum wavefront phase elapsed timers including barriers; not CPU-seconds. Dense typed closed phase also builds eager branch rows. Phase sum excludes teardown/counter setup/collection and other call overhead.',
                'split_visits':'List elements entered in multiloop scans, including entries rejected by endpoint/start checks. Typed scans can include a terminating future-end entry. Not candidate count Z.',
                'split_eligible':'Scalar/sparse: visited entries with start layer beyond interval start. Typed: excludes direct/same-layer/future-end cases. Pruning comparison uses scalar/sparse in the same right-normal loop.',
                'visit_ratios':'Dense-scalar counter divided by sparse counter for a matching profile case; null when denominator is zero. Not wall speedup or candidate retention fraction.',
                'peak_rss_gib':'GNU maximum resident set size KiB / 2^20; whole-process memory, distinct from the 32GiB virtual address-space limit.',
                'swapping':'Host-wide /proc/vmstat pswpin/pswpout differences. Positive values cannot be attributed to this child; zeros show no host swapping during that window.',
                'partial_outcomes':'Missing, timeout, failed and invalid-output outcomes remain in every planned grid. Timeouts are censored elapsed observations, not completed solve times. Successful-subset memory/timing summaries retain explicit counts.',
                'provenance':'Source/input/output hashes and recorded commands are checked; present binaries are rehashed and absent binaries may be omitted for redistribution. Actual inherited environment, affinity success and exclusive hardware use are not independently authenticated by these records.'}}



def export_tables(summary, out):
    out.mkdir(parents=True,exist_ok=True)
    for name,rows in [('status',summary['outcomes']),('ablation-aggregate',summary['ablation_aggregate']),('profile-comparisons',summary['profile_comparisons']),('scaling',summary['scaling']),('workstation',summary['workstation'])]:
        with (out/(name+'.csv')).open('w',newline='') as stream:
            if rows:
                writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),lineterminator='\n')
                writer.writeheader()
                writer.writerows({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
    end=r'\\'
    def rows_file(name,lines):
        (out/name).write_text('% Generated by analyze_publication_campaign.py; do not edit.\n'+'\n'.join(lines)+'\n')
    def value(x,places=3):
        return f'{x:.{places}f}' if x is not None else '---'
    def tex(text):
        return text.replace('_',r'\_')
    rows_file('status-rows.tex',[f"{phase.capitalize()} & {r['expected']} & {r['ok']} & {r['failed']+r['invalid_output']} & {r['timeout']} & {r['missing']} "+end for phase,r in summary['phase_completion'].items()])
    lines=[]
    for r in summary['ablation_aggregate']:
        label={'dense_typed_over_sparse':'Typed/sparse','dense_scalar_over_sparse':'Scalar/sparse','dense_typed_over_dense_scalar':'Typed/scalar'}[r['ratio']]
        ci=r['conditional_repeat_ci95']; interval=f'[{ci[0]:.3f}, {ci[1]:.3f}]' if ci else '---'
        partial='' if r['complete_panel'] else ' (partial)'
        lines.append(f"{r['lambda']} & {r['threads']} & {label} & {r['complete_cases']}/{r['expected_cases']}{partial} & {value(r['geomean_of_paired_medians'])} & {interval} "+end)
    rows_file('controlled-rows.tex',lines)
    lines=[]
    for lam,threads in sorted({(r['task']['lambda'],r['task']['threads']) for r in summary['profiles']}):
        group=[r for r in summary['profile_comparisons'] if (r['lambda'],r['threads'])==(lam,threads)]
        expected=len({r['task']['accession'] for r in summary['profiles'] if (r['task']['lambda'],r['task']['threads'])==(lam,threads)})
        vals=[]
        for field in ('scalar_over_sparse_eligible_visits','scalar_over_sparse_multiloop_seconds','typed_over_sparse_multiloop_seconds'):
            data=[r[field] for r in group if r[field] is not None]
            vals.append(value(statistics.median(data) if data else None))
        lines.append(f'{lam} & {threads} & {len(group)}/{expected} & '+' & '.join(vals)+' '+end)
    rows_file('profile-rows.tex',lines)
    rows_file('scaling-rows.tex',[f"{tex(r['accession'])} & {r['lambda']} & {r['threads']} & {r['wall_seconds']['n'] if r['wall_seconds'] else 0}/{r['expected']} & {value(r['wall_seconds']['median'] if r['wall_seconds'] else None,2)} & {value(r['speedup_ratio_of_medians'],2)} "+end for r in summary['scaling']])
    rows_file('workstation-rows.tex',[f"{r['layout'].capitalize()} & {r['wall_seconds']['n'] if r['wall_seconds'] else 0}/{r['expected']} & {value(r['wall_seconds']['median'] if r['wall_seconds'] else None,2)} & {value(r['rss_gib']['median'] if r['rss_gib'] else None)} & {r['outcomes'].get('timeout',0)} & {r['outcomes'].get('failed',0)+r['outcomes'].get('invalid_output',0)} & {r['outcomes']['missing']} "+end for r in summary['workstation']])
    claims={'TimingPlannedTasks':sum(r['expected'] for r in summary['phase_completion'].values()),'TimingSuccessfulTasks':sum(r['ok'] for r in summary['phase_completion'].values()),'TimingMissingTasks':len(summary['missing'])}
    rows_file('claims.tex',['\\newcommand{\\'+k+'}{'+str(v)+'}' for k,v in claims.items()])
    (out/'METHODS.md').write_text('\n\n'.join('**'+k.replace('_',' ')+'.** '+v for k,v in summary['measurement_definitions'].items())+'\n\nBootstrap specification:\n\n```json\n'+json.dumps(summary['bootstrap'],indent=2)+'\n```\n')


def plot(summary, panel, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
    lengths = {r['accession']:r['aa_length'] for r in panel['panel']}
    fig, axes = plt.subplots(2,2,figsize=(7.0,4.5),sharex=True,sharey=True)
    for row,lam in enumerate((0,4)):
        for col,threads in enumerate((1,16)):
            ax=axes[row,col]
            for ratio,label,color,marker in [('dense_typed_over_sparse','Typed dense / sparse','#0072B2','o'),
                                            ('dense_scalar_over_sparse','Scalar dense / sparse','#D55E00','s')]:
                cases=sorted([c for c in summary['ablation_cases'] if c['lambda']==lam and c['threads']==threads and c['ratios'][ratio]],key=lambda c:lengths[c['accession']])
                x=[lengths[c['accession']] for c in cases]
                y=[c['ratios'][ratio]['median'] for c in cases]
                lower=[v-c['ratios'][ratio]['min'] for v,c in zip(y,cases)]
                upper=[c['ratios'][ratio]['max']-v for v,c in zip(y,cases)]
                ax.errorbar(x,y,yerr=[lower,upper],fmt=marker,ms=3,capsize=2,color=color,label=label,alpha=.85)
            ax.axhline(1,color='.5',lw=.7,ls='--')
            planned=[c for c in summary['ablation_cases'] if c['lambda']==lam and c['threads']==threads]
            complete=sum(c['complete_triplets']==c['expected_repetitions'] for c in planned)
            ax.set_title(f'λ = {lam}, j = {threads}; {complete}/{len(planned)} complete cases')
            ax.set_xlim(0,max(lengths.values())*1.04)
            if not any(c['complete_triplets'] for c in planned):
                ax.text(.5,.75,'No complete three-arm repetitions',ha='center',transform=ax.transAxes,fontsize=8,color='.4')
            ax.grid(alpha=.18)
            if col==0:ax.set_ylabel('Paired solve-time ratio')
            if row==1:ax.set_xlabel('Protein length (aa)')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False)
    fig.tight_layout(rect=(0,.035,1,.94))
    fig.text(.5,.01,'Points: medians; whiskers: observed ranges of complete three-arm repetitions.',ha='center',fontsize=8)
    fig.savefig(out/'controlled-ablation.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(out/'controlled-ablation.png',dpi=160)
    plt.close(fig)


@contextmanager
def analysis_lock(path,enabled):
    if enabled:
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('a+') as stream:
            fcntl.flock(stream,fcntl.LOCK_EX)
            yield
    else:
        yield


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--require-complete',action='store_true',help='require every planned outcome; recorded failures remain valid outcomes')
    parser.add_argument('--bootstrap-draws',type=int,default=10000)
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--lock',type=Path,help='optional shared compute lock; omit when the calling controller already holds it')
    args=parser.parse_args()
    out=args.output or args.input/'analysis'
    try:
        require(args.bootstrap_draws>=0,'bootstrap draws must be nonnegative')
        require(not out.resolve().is_relative_to((args.input/'timings/runs').resolve()),'output may not overwrite raw run directories')
        with analysis_lock(args.lock,args.lock is not None):
            result=analyze(args.input,args.require_complete,args.bootstrap_draws)
            write_json(out/'summary.json',result)
            export_tables(result,out)
            if args.plot:plot(result,read_json(args.input/'inputs/panel.json'),out)
        print(json.dumps({k:result[k] for k in ('complete','all_successful','phase_completion','maximum_objective_delta_units')},indent=2))
        return 0
    except (ValueError,OSError,KeyError,TypeError,IndexError) as exc:
        print(f'FAIL: {exc}',file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
