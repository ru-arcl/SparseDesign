#!/usr/bin/env python3
"""Reproduce the EPYC CLI checks with frozen verification sources; no new designs."""
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'cli-validation/source'))
spec = importlib.util.spec_from_file_location('analyze_publication_campaign', ROOT / 'timings/analyze_node1.py')
strict = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = strict
spec.loader.exec_module(strict)
import publication_cli_validation

if __name__ == '__main__':
    raise SystemExit(publication_cli_validation.main())
