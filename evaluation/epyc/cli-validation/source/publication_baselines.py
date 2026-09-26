#!/usr/bin/env python3
"""Pinned, capped, native-serial public baseline evaluations (no vendored sources).

Build unmodified upstream sources using their native build files first. See the
adjacent baseline protocol for source pins, objective units, and interpretation.
Every process is run under GNU time, a virtual-address limit, a wall timeout,
single-CPU affinity, and a shared advisory lock. Existing records are resumed.
"""
import argparse
import csv
import fcntl
import hashlib
import json
import itertools
import io
import math
import os
from pathlib import Path
import platform
import random
import re
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time

PINS = {
    'lineardesign': ('https://github.com/LinearDesignSoftware/LinearDesign.git', 'f0126ca89a8b853088b4bccfd2cc8c378d3678be'),
    'linearcdsfold': ('https://github.com/ablab-nthu/LinearCDSfold.git', 'de1db7e6c8d0a4360bf4f1c9270f11525ad51ff6'),
    'derna': ('https://github.com/elkebir-group/derna.git', 'adecda8277add8c0b75512434bd729f30b8fd7d6'),
}
STANDARD_CODE = dict(zip((''.join(c) for c in itertools.product('UCAG', repeat=3)),
                        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def preserve_or_check(path, value):
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError(f'Resume configuration differs: {path}; use a new output directory')
    else:
        write_json(path, value)


def validate_output(parsed, protein, weights):
    rna, structure = parsed.get('rna', ''), parsed.get('structure', '')
    checks = dict(rna_length=len(rna) == 3*len(protein),
                  structure_length=len(structure) == len(rna),
                  translation=''.join(STANDARD_CODE.get(rna[i:i+3], '?') for i in range(0, len(rna), 3)) == protein)
    stack, paired, balanced = [], True, True
    for i, ch in enumerate(structure):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if not stack:
                balanced = False
                break
            j = stack.pop()
            paired &= i < len(rna) and rna[j]+rna[i] in {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'} and i-j > 3
        elif ch != '.':
            balanced = False
    checks.update(dot_bracket_balanced=balanced and not stack,
                  canonical_pairs_and_minimum_hairpin=paired)
    result = dict(validation_checks=checks, translation_valid=checks['translation'])
    if checks['translation']:
        penalty = -sum(math.log(weights[rna[i:i+3]]) for i in range(0, len(rna), 3))
        result.update(common_codon_penalty=penalty, common_cai=math.exp(-penalty/len(protein)))
    return result

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def parse_output(tool, text):
    if tool == 'sparsedesign':
        tool = 'lineardesign'
    patterns = {
        'lineardesign': (r'mRNA sequence:\s+([ACGU]+)', r'mRNA structure:\s+([().]+)', r'folding free energy:\s+([-\d.]+)'),
        'linearcdsfold': (r'(?:mRNA sequence|RNA sequence|Sequence|sequence):\s*([ACGU]+)', r'(?:mRNA structure|RNA structure|Structure|structure):\s*([().]+)', r'(?:MFE|Minimum Free Energy|[Ff]olding free energy):\s*([-\d.]+)'),
        'derna': (r'zuker (?:cai )?rna:\s*([ACGU]+)[,.]', r'zuker (?:cai )?bp:\s*([().]+)', r'(?:Free Energy|Energy):\s*([-\d.]+)'),
    }
    ans = {}
    for key, pattern in zip(('rna', 'structure', 'native_mfe'), patterns[tool]):
        values = re.findall(pattern, text)
        if values:
            ans[key] = float(values[-1]) if key == 'native_mfe' else values[-1]
    # LinearCDSfold writes bare RNA and dot-bracket lines under result labels.
    if tool == 'linearcdsfold':
        for key, pattern in [('rna', r'^([ACGU]{9,})\s*$'), ('structure', r'^([().]{9,})\s*$')]:
            values = re.findall(pattern, text, re.M)
            if values:
                ans[key] = values[-1]
    return ans

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--panel', type=Path, required=True)
    ap.add_argument('--table', type=Path, required=True)
    ap.add_argument('--sources', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--lock', type=Path, required=True)
    ap.add_argument('--tools', nargs='+', choices=list(PINS)+['sparsedesign'], default=list(PINS))
    ap.add_argument('--sparsedesign-binary', type=Path)
    ap.add_argument('--sparsedesign-source', type=Path)
    ap.add_argument('--stop-after-cap', action='store_true', help='One ascending feasibility run; stop each tool/lambda length chain after its first resource cap')
    ap.add_argument('--feasibility', type=Path, help='Only repeat successful tool/accession/lambda conditions from this complete feasibility directory')
    ap.add_argument('--accessions', nargs='+', help='Restrict to identifiers already in the supplied panel')
    ap.add_argument('--lambdas', nargs='+', type=float, default=[0, 4])
    ap.add_argument('--max-aa', type=int, default=100000)
    ap.add_argument('--min-aa', type=int, default=0)
    ap.add_argument('--repetitions', type=int, default=5)
    ap.add_argument('--timeout', type=float, default=600)
    ap.add_argument('--memory-gib', type=float, default=32)
    ap.add_argument('--cpu', type=int, default=0)
    ap.add_argument('--seed', type=int, default=20260906)
    args = ap.parse_args()
    if len(args.tools) != len(set(args.tools)) or len(args.lambdas) != len(set(args.lambdas)):
        ap.error('Tools and lambdas must not contain duplicates')
    if args.timeout <= 0 or args.memory_gib <= 0 or args.repetitions < 1:
        ap.error('Limits and repetitions must be positive')
    if args.cpu not in os.sched_getaffinity(0):
        ap.error('Selected CPU is outside current affinity')
    if any(x < 0 or not math.isfinite(x) for x in args.lambdas):
        ap.error('Lambdas must be finite and nonnegative')
    if args.stop_after_cap and (args.repetitions != 1 or args.feasibility):
        ap.error('--stop-after-cap requires one repetition and no --feasibility')
    if 'sparsedesign' in args.tools and (not args.sparsedesign_binary or not args.sparsedesign_source):
        ap.error('SparseDesign requires explicit frozen binary and source paths')
    for key in ('panel', 'table', 'sources', 'output', 'lock'):
        setattr(args, key, getattr(args, key).resolve())
    args.output.mkdir(parents=True, exist_ok=True)
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    panel = json.loads(args.panel.read_text())['panel']
    if args.accessions and not set(args.accessions) <= {x['accession'] for x in panel}:
        ap.error('An accession is absent from the supplied panel')
    panel = [x for x in panel if args.min_aa <= x['aa_length'] <= args.max_aa and
             (not args.accessions or x['accession'] in args.accessions)]
    if not panel:
        ap.error('Panel selection is empty')
    provenance = {}
    executables = {'lineardesign': 'bin/LinearDesign_2D', 'linearcdsfold': 'LinearCDSfold', 'derna': 'build/derna'}
    for tool in args.tools:
        if tool == 'sparsedesign':
            args.sparsedesign_binary = args.sparsedesign_binary.resolve()
            args.sparsedesign_source = args.sparsedesign_source.resolve()
            provenance[tool] = dict(url='local frozen study source', commit=None,
                                   executable_sha256=sha(args.sparsedesign_binary),
                                   tracked_file_sha256={str(p.relative_to(args.sparsedesign_source)): sha(p) for p in sorted(args.sparsedesign_source.rglob('*')) if p.is_file()})
            continue
        source = args.sources / tool
        head = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        if head != PINS[tool][1]:
            raise ValueError(f'{tool}: expected {PINS[tool][1]}, got {head}')
        if subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD', '--'], text=True):
            raise ValueError(f'{tool}: tracked source differs from pinned upstream')
        tracked = subprocess.check_output(['git', '-C', str(source), 'ls-files', '-z']).decode().split('\0')
        runtime = {'coding_wheel.txt', 'src/Utils/libraries/LinearDesign_linux64.so', 'src/Utils/libraries/LinearDesign_linux64_old.so'} if tool == 'lineardesign' else set()
        provenance[tool] = dict(url=PINS[tool][0], commit=head, executable_sha256=sha(source / executables[tool]),
                               tracked_file_sha256={p: sha(source/p) for p in tracked if p and (source/p).is_file()},
                               runtime_dependency_sha256={p: sha(source/p) for p in sorted(runtime)})
    # Lossless rearrangement of input decimal frequency strings into DERNA format.
    byaa = {}
    frequencies = {}
    translated = {}
    for row in csv.reader(args.table.open()):
        if len(row) >= 3 and row[0].strip() in {a+b+c for a in 'ACGU' for b in 'ACGU' for c in 'ACGU'}:
            codon, aa, freq = (s.strip() for s in row[:3])
            frequencies[codon] = float(freq)
            translated[codon] = aa
            if aa != '*':
                byaa.setdefault(aa, []).append((codon, freq))
    if translated != STANDARD_CODE or any(v <= 0 for v in frequencies.values()):
        raise ValueError('A complete positive-frequency standard genetic-code table is required')
    maxima = {aa: max(float(f) for _, f in vals) for aa, vals in byaa.items()}
    weights = {c: f/maxima[translated[c]] for c, f in frequencies.items() if translated[c] != '*'}
    derna_table = args.output / 'codon-table-derna.csv'
    fh = io.StringIO(newline='')
    # DERNA checks header token length == 3 without stripping CR. Its parser
    # silently drops the final codon with csv.writer's default CRLF records.
    writer = csv.writer(fh, lineterminator='\n')
    for aa, vals in sorted(byaa.items()):
        writer.writerow([''] + [v[0] for v in vals])
        writer.writerow([aa] + [v[1] for v in vals])
    encoded_table = fh.getvalue().encode()
    if derna_table.exists() and derna_table.read_bytes() != encoded_table:
        raise ValueError(f'Resume table differs: {derna_table}')
    derna_table.write_bytes(encoded_table)
    metadata = dict(schema='sparsedesign-public-baselines-v2',
                    pins=provenance, panel_sha256=sha(args.panel), table_sha256=sha(args.table),
                    panel_selection=[x['accession'] for x in panel], tools=args.tools, lambdas=args.lambdas,
                    stop_after_cap=args.stop_after_cap,
                    derna_table_sha256=sha(derna_table), wrapper_sha256=sha(__file__),
                    python=sys.version, platform=platform.platform(),
                    shared_lock=str(args.lock),
                    memory_limit_kind='RLIMIT_AS virtual address bytes', memory_limit_bytes=int(args.memory_gib*2**30),
                    timeout_seconds=args.timeout, cpu=args.cpu, repetitions=args.repetitions, seed=args.seed,
                    timing_boundary='native CLI process including parsing, allocation, traceback, native output and any native postprocessing')
    feasible = None
    if args.feasibility:
        feasibility_records = [json.loads(p.read_text()) for p in sorted(args.feasibility.glob('*/record.json'))]
        planned = json.loads((args.feasibility/'planned-order.json').read_text())
        if len(feasibility_records) != len(planned):
            raise ValueError('Feasibility grid is incomplete')
        feasible = {(r['accession'], r['lambda_value'], r['tool']) for r in feasibility_records if r['status'] == 'ok'}
        metadata['feasibility_record_sha256'] = {str(p.relative_to(args.feasibility)): sha(p) for p in sorted(args.feasibility.glob('*/record.json'))}
    preserve_or_check(args.output / 'metadata.json', metadata)
    tasks = [(rep, entry, lam, tool) for rep in range(1, args.repetitions+1) for entry in panel for lam in args.lambdas for tool in args.tools]
    if feasible is not None:
        tasks = [(r,e,l,t) for r,e,l,t in tasks if (e['accession'],l,t) in feasible]
    random.Random(args.seed).shuffle(tasks)
    if args.stop_after_cap:
        tasks.sort(key=lambda x: x[1]['aa_length'])
    manifest = [dict(repetition=r, accession=e['accession'], lambda_value=l, tool=t) for r,e,l,t in tasks]
    preserve_or_check(args.output / 'planned-order.json', manifest)
    capped = {}
    for rep, entry, lam, tool in tasks:
        label = f"{entry['accession']}.lambda{lam:g}.{tool}.rep{rep}"
        dest = args.output / label
        if (dest / 'record.json').exists():
            previous = json.loads((dest/'record.json').read_text())
            for key in ('wrapper_sha256', 'table_sha256'):
                if previous[key] != metadata[key]:
                    raise ValueError(f'Resume provenance differs: {dest}/{key}')
            for name, expected in previous['raw_file_sha256'].items():
                if sha(dest/name) != expected:
                    raise ValueError(f'Resume raw digest differs: {dest/name}')
            if previous['status'] in ('timeout', 'allocation_failure'):
                capped[tool, lam] = previous['label']
            continue
        fasta = args.panel.parent.parent / entry['fasta_path']
        if sha(fasta) != entry['fasta_sha256']:
            raise ValueError(f'FASTA digest mismatch: {fasta}')
        protein = ''.join(x.strip() for x in fasta.read_text().splitlines() if not x.startswith('>'))
        if len(protein) != entry['aa_length'] or hashlib.sha256(protein.encode()).hexdigest() != entry['sequence_sha256']:
            raise ValueError(f'Protein digest or length mismatch: {fasta}')
        if dest.exists():
            dest.rename(dest.with_name(dest.name + '.interrupted.' + str(time.time_ns())))
        dest.mkdir(exist_ok=True)
        if args.stop_after_cap and (tool, lam) in capped:
            write_json(dest/'record.json', dict(label=label,tool=tool,accession=entry['accession'],aa_length=entry['aa_length'],
                       lambda_value=lam,repetition=rep,status='not_attempted_after_cap',trigger=capped[tool,lam],
                       wrapper_sha256=metadata['wrapper_sha256'],table_sha256=metadata['table_sha256'],raw_file_sha256={}))
            continue
        source = args.sources / tool
        binary = args.sparsedesign_binary if tool == 'sparsedesign' else source / executables[tool]
        with tempfile.TemporaryDirectory(prefix='sparsedesign-baseline-') as tmp:
            cwd = Path(tmp)
            result_path = dest / 'native-result.txt'
            if tool == 'lineardesign':
                shutil.copyfile(source / 'coding_wheel.txt', cwd / 'coding_wheel.txt')
                # Native upstream binary records a relative DT_NEEDED pathname.
                # Preserve its lookup layout without modifying the binary/source.
                (cwd / 'src').symlink_to(source / 'src', target_is_directory=True)
                cmd = [str(binary), format(lam, 'g'), '1', str(args.table)]
            elif tool == 'linearcdsfold':
                cmd = [str(binary), '-m', 'exact', '-O', 'LD', '-l', format(lam, 'g'), '-c', str(args.table), '-o', str(result_path), '-f', str(dest / 'native-result.csv'), str(fasta)]
            elif tool == 'derna':
                cmd = [str(binary), '-i', str(fasta), '-o', str(result_path), '-m', '1', '-s', '1' if lam == 0 else '2', '-c', str(derna_table)]
                if lam:
                    cmd += ['-l', format(1/(1+100*lam), '.17g')]
            else:
                cmd = [str(binary), '-l', format(lam, 'g'), '-c', str(args.table), '-j', '1', '--diagnostics']
            timed = ['/usr/bin/time', '-f', 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kib=%M\nexit_status=%x', '-o', str(dest/'time.txt'), *cmd]
            env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
            def limits():
                resource.setrlimit(resource.RLIMIT_AS, (int(args.memory_gib*2**30),)*2)
                resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
                os.sched_setaffinity(0, {args.cpu})
            with args.lock.open('a') as lock, fasta.open('rb') as inp, (dest/'stdout.txt').open('wb') as out, (dest/'stderr.txt').open('wb') as err:
                fcntl.flock(lock, fcntl.LOCK_EX)
                started = time.time()
                observer_start = time.monotonic()
                host_load = Path('/proc/loadavg').read_text().strip()
                process = subprocess.Popen(timed, cwd=cwd, env=env, stdin=inp, stdout=out, stderr=err, preexec_fn=limits, start_new_session=True)
                status = 'ok'
                try:
                    rc = process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    status = 'timeout'
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        rc = process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        rc = process.wait()
                elapsed = time.monotonic()-observer_start
                fcntl.flock(lock, fcntl.LOCK_UN)
            raw = (dest/'stdout.txt').read_text(errors='replace')
            if result_path.exists():
                raw += '\n'+result_path.read_text(errors='replace')
            parsed = parse_output(tool, raw)
            if rc and status == 'ok':
                stderr_text = (dest/'stderr.txt').read_text(errors='replace')
                status = ('wrapper_error' if 'error while loading shared libraries' in stderr_text or rc in (126, 127)
                          else 'allocation_failure' if 'bad_alloc' in stderr_text else 'error')
            if status == 'ok' and (len(parsed.get('rna','')) != entry['aa_length']*3 or len(parsed.get('structure','')) != entry['aa_length']*3):
                status = 'invalid_or_unparsed_output'
            verification = validate_output(parsed, protein, weights)
            if status == 'ok' and not all(verification['validation_checks'].values()):
                status = 'invalid_output'
            metrics = {}
            if (dest/'time.txt').exists():
                for line in (dest/'time.txt').read_text().splitlines():
                    if '=' in line:
                        k, v = line.split('=',1)
                        try:
                            metrics[k] = float(v)
                        except ValueError:
                            pass
            generated_bytes = sum(p.stat().st_size for p in cwd.iterdir() if p.is_file() and p.suffix == '.bin')
            record = dict(label=label,tool=tool,accession=entry['accession'],aa_length=entry['aa_length'],lambda_value=lam,repetition=rep,
                          status=status,returncode=rc,elapsed_observer_seconds=elapsed,unix_started=started,command=cmd,
                          source_commit=provenance[tool]['commit'],executable_sha256=provenance[tool]['executable_sha256'],
                          input_fasta_sha256=entry['fasta_sha256'], table_sha256=metadata['table_sha256'],
                          wrapper_sha256=metadata['wrapper_sha256'], host_loadavg_started=host_load,
                          native_mfe_kind='fixed-sequence refold' if tool == 'derna' and lam else 'native reported folding score',
                          native_temporary_bin_bytes=generated_bytes,**metrics,**parsed,**verification)
            record['raw_file_sha256'] = {p.name: sha(p) for p in dest.iterdir() if p.is_file()}
            write_json(dest/'record.json', record)
            if status in ('timeout', 'allocation_failure'):
                capped[tool, lam] = label
            print(json.dumps({k:record.get(k) for k in ('label','status','wall_seconds','max_rss_kib','native_mfe')}),flush=True)
            if status in ('wrapper_error', 'invalid_or_unparsed_output'):
                raise RuntimeError(f'{label}: wrapper or parser failure; inspect retained evidence before continuing')

if __name__ == '__main__':
    main()
