#!/usr/bin/env python3
"""Write/check the active solver source manifest, independent of evaluation snapshots."""
from pathlib import Path
import argparse
import hashlib

ROOT = Path(__file__).resolve().parents[1]
TOP = ['.gitattributes', '.gitignore', 'BENCHMARKS.md', 'CITATION.cff', 'DESIGN.md', 'LICENSE.md', 'LICENSE.docx', 'Makefile', 'NOTICE',
       'PROVENANCE.md', 'README.md', 'RELEASE.md', 'THIRD_PARTY_NOTICES.md']
TREES = ['src', 'test', 'data', 'tools', 'LICENSES']


def render():
    paths = [ROOT / name for name in TOP]
    for name in TREES:
        for path in (ROOT / name).rglob('*'):
            if '__pycache__' in path.parts or path.suffix in {'.pyc', '.pyo', '.o'}:
                continue
            if path.is_symlink():
                raise ValueError('Source manifest does not allow symlinks: ' + str(path))
            if path.is_file():
                paths.append(path)
    result = []
    for path in sorted(paths, key=lambda p: p.relative_to(ROOT).as_posix()):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        result.append(digest + '  ' + path.relative_to(ROOT).as_posix() + '\n')
    return ''.join(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='fail on drift without modifying files')
    args = parser.parse_args()
    expected = render()
    manifest = ROOT / 'SOURCE_SHA256SUMS'
    if args.check:
        if not manifest.exists() or manifest.read_text() != expected:
            raise SystemExit('SOURCE_SHA256SUMS differs; review changes before regenerating')
        print('Active source manifest matches (' + str(len(expected.splitlines())) + ' files)')
    else:
        manifest.write_text(expected)
        print('Wrote ' + str(len(expected.splitlines())) + ' source checksums')


if __name__ == '__main__':
    main()
