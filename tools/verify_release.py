#!/usr/bin/env python3
"""Verify the complete source-and-evidence release before building in the checkout."""
import argparse
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'RELEASE_SHA256SUMS'


def release_files():
    result = []
    for path in ROOT.rglob('*'):
        relative = path.relative_to(ROOT)
        if (path == MANIFEST
                or any(part in {'.git', '__pycache__', '.pytest_cache'}
                       or part.startswith('.venv') for part in relative.parts)
                or path.suffix in {'.pyc', '.pyo'}):
            continue
        if path.is_symlink():
            raise ValueError('Unexpected symbolic link: ' + relative.as_posix())
        if path.is_file():
            with path.open('rb') as stream:
                magic = stream.read(8)
            if magic.startswith((b'\x7fELF', b'MZ', b'!<arch>\n')):
                raise ValueError('Native binary in source release: ' + relative.as_posix()
                                 + '; verify before building, or run make clean first')
            result.append(path)
    return sorted(result, key=lambda path: path.relative_to(ROOT).as_posix())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true',
                        help='refresh after reviewing intentional package changes')
    args = parser.parse_args()
    paths = release_files()
    contents = ''.join(hashlib.sha256(path.read_bytes()).hexdigest() + '  '
                       + path.relative_to(ROOT).as_posix() + '\n' for path in paths)
    if args.write:
        MANIFEST.write_bytes(contents.encode('utf-8'))
    elif not MANIFEST.exists() or MANIFEST.read_text(encoding='utf-8') != contents:
        raise SystemExit('Release contents or membership differ; inspect before using --write')
    print(('Wrote' if args.write else 'Verified') + f' {len(paths)} release files')


if __name__ == '__main__':
    main()
