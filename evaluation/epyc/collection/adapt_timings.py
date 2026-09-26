#!/usr/bin/env python3
"""Adapt a freshly built September timing study to run on EPYC NUMA node 1.

Run once, after research/build_publication_campaign.py has created the study
directory and before any measurement. Changes, all recorded in
epyc-node1-adaptation.json:

- collector: CPU list j -> 16..16+j-1 instead of 0..j-1, and every solve
  (and GNU time) gets MPOL_BIND to node 1 via set_mempolicy in the child;
- protocol: every wall cap doubled (EPYC cores are slower than i9 P-cores),
  affinity text updated, host_adaptation recorded; tasks, seed, order,
  OpenMP settings and the 32-GiB address-space cap are unchanged;
- build.json: source hash for the patched collector refreshed;
- analyzer: a copy whose only change is the expected CPU list.
"""
import argparse
import hashlib
import json
from pathlib import Path

CPU_BASE = 16
NODE = 1
CAP_FACTOR = 2


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_once(text, old, new, what):
    if text.count(old) != 1:
        raise SystemExit(f'patch point not found exactly once: {what}')
    return text.replace(old, new)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, required=True)
    ap.add_argument('--analyzer', type=Path, required=True, help='publication/code/research/analyze_publication_campaign.py')
    args = ap.parse_args()
    study = args.study.resolve()
    record = study / 'epyc-node1-adaptation.json'
    if record.exists():
        raise SystemExit('Study already adapted.')
    if (study / 'timings').exists():
        raise SystemExit('Study already has measurements; refusing to change its protocol.')
    before = {name: sha(study / name) for name in
              ('source/publication_campaign_final.py', 'protocol.json', 'build.json')}

    collector = study / 'source/publication_campaign_final.py'
    text = collector.read_text()
    text = replace_once(text, "import argparse\n",
                        "import argparse\nimport ctypes\n", 'imports')
    text = replace_once(text, "MODES = ('dense_typed', 'dense_scalar', 'sparse')\n",
                        "MODES = ('dense_typed', 'dense_scalar', 'sparse')\n"
                        f"CPU_BASE = {CPU_BASE}  # EPYC 7313 NUMA node {NODE}: logical CPUs 16-31\n"
                        f"MEMORY_NODE = {NODE}\n", 'constants')
    text = replace_once(text, "    cpus = list(range(task['threads']))\n",
                        "    cpus = [CPU_BASE + i for i in range(task['threads'])]\n", 'cpu list')
    text = replace_once(text,
                        "        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))\n",
                        "        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))\n"
                        "        mask = ctypes.c_ulong(1 << MEMORY_NODE)  # MPOL_BIND=2, SYS_set_mempolicy=238 (x86_64)\n"
                        "        if ctypes.CDLL(None, use_errno=True).syscall(238, 2, ctypes.byref(mask), 64) != 0:\n"
                        "            os._exit(126)\n", 'mempolicy')
    text = replace_once(text, "    if args.action == 'run':\n",
                        "    if args.action == 'run':\n"
                        "        needed = set(range(CPU_BASE, CPU_BASE + 16))\n"
                        "        if not needed <= os.sched_getaffinity(0):\n"
                        "            raise SystemExit(f'CPUs {min(needed)}-{max(needed)} are not all available')\n", 'cpu check')
    collector.write_text(text)

    protocol_path = study / 'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    original_caps = dict(protocol['phase_timeout_seconds'])
    protocol['phase_timeout_seconds'] = {k: v * CAP_FACTOR for k, v in original_caps.items()}
    protocol['timeout_seconds'] *= CAP_FACTOR
    protocol['workstation_timeout_seconds'] *= CAP_FACTOR
    protocol['affinity'] = (f'nested logical CPUs {CPU_BASE}..{CPU_BASE}+threads-1 on NUMA node {NODE} of a '
                            'dual-socket AMD EPYC 7313 (16 homogeneous cores per node, no SMT); '
                            f'memory bound to node {NODE} with MPOL_BIND')
    protocol['host_adaptation'] = {
        'host': 'AMD EPYC 7313 x2 (arrakis), NUMA node 1',
        'source_protocol_sha256': before['protocol.json'],
        'changes': ['CPU list offset by 16', 'MPOL_BIND to node 1',
                    f'all wall caps multiplied by {CAP_FACTOR} (original {original_caps})'],
        'unchanged': 'tasks, order, seed, repetitions, OpenMP settings, 32-GiB RLIMIT_AS, kernels and sources',
        'compiler_note': 'binaries rebuilt with this host\'s g++; the i9 build used GCC 14.3.0',
    }
    protocol_path.write_text(json.dumps(protocol, indent=2, sort_keys=True) + '\n')

    build_path = study / 'build.json'
    build = json.loads(build_path.read_text())
    build['source_sha256']['publication_campaign_final.py'] = sha(collector)
    build['epyc_node1_adaptation'] = 'see epyc-node1-adaptation.json'
    build_path.write_text(json.dumps(build, indent=2) + '\n')

    analyzer_text = args.analyzer.read_text()
    analyzer_text = replace_once(
        analyzer_text,
        "require(cmd[6]==','.join(str(i) for i in range(task['threads']))",
        f"require(cmd[6]==','.join(str({CPU_BASE}+i) for i in range(task['threads']))",
        'analyzer affinity check')
    analyzer = study / 'analyze_node1.py'
    analyzer.write_text(analyzer_text)

    after = {name: sha(study / name) for name in before}
    record.write_text(json.dumps({
        'cpu_base': CPU_BASE, 'memory_node': NODE, 'cap_factor': CAP_FACTOR,
        'original_caps': original_caps, 'before_sha256': before, 'after_sha256': after,
        'analyzer_source_sha256': sha(args.analyzer), 'analyzer_node1_sha256': sha(analyzer),
    }, indent=2) + '\n')
    print('Adapted', study)


if __name__ == '__main__':
    main()
