#!/usr/bin/env python3
"""Run the frozen CLI verifier with the EPYC-adapted strict record validator."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args, _ = parser.parse_known_args()
    root = Path(__file__).resolve().parents[2]
    source = root / 'publication/code/research'
    sys.path.insert(0, str(source))
    strict = load('analyze_publication_campaign', args.input / 'analyze_node1.py')
    verifier = load('publication_cli_validation', source / 'publication_cli_validation.py')
    status = verifier.main()
    provenance = args.output / 'source'
    provenance.mkdir()
    for path in (Path(__file__), Path(strict.__file__), Path(verifier.__file__),
                 Path(verifier.baseline.__file__), Path(verifier.native.__file__)):
        shutil.copyfile(path, provenance / path.name)
    (args.output / 'invocation.json').write_text(json.dumps(dict(
        argv=sys.argv,
        adaptation='Use analyze_node1.py for EPYC affinity, overlapping permitted j1 windows, and doubled caps.',
        partial_reason='Five excluded pilot records are intentionally missing. All 70 selected CLI records are required to be present.'
    ), indent=2) + '\n')
    summary = json.loads((args.output / 'summary.json').read_text())
    assert summary['planned'] == summary['recorded'] == summary['successful_outputs_checked'] == 70
    assert summary['all_successful_outputs_pass'] and not summary['failed_check_task_ids']
    files = sorted(path for path in args.output.rglob('*') if path.is_file() and path.name != 'SHA256SUMS')
    (args.output / 'SHA256SUMS').write_text(''.join(
        f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(args.output)}\n' for path in files))
    return status


if __name__ == '__main__':
    raise SystemExit(main())
