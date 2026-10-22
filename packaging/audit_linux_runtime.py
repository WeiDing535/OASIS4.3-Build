"""Audit a staged relocatable Linux x86_64 runtime and inventory its bytes.

Run only on a staging copy, never on the original installation. No packages are
installed and no binaries are patched. Non-system dependencies must resolve
inside this prefix without LD_LIBRARY_PATH. This does not grant redistribution.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess

SYSTEM_LIBRARIES = {'libc.so.6', 'libm.so.6', 'libdl.so.2', 'libpthread.so.0',
                    'librt.so.1', 'libresolv.so.2', 'libutil.so.1', 'libmvec.so.1',
                    'ld-linux-x86-64.so.2'}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def parse_ldd(text):
    dependencies = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(('linux-vdso.', 'statically linked')):
            continue
        if 'not found' in line:
            raise ValueError('Missing Linux library: ' + line)
        if '=>' in line:
            name, target = line.split('=>', 1)
            target = target.strip().split(' (', 1)[0]
            if not target.startswith('/'):
                raise ValueError('Unresolved Linux library: ' + line)
            dependencies.append((name.strip(), Path(target)))
        elif line.startswith('/'):
            target = line.split(' (', 1)[0]
            dependencies.append((Path(target).name, Path(target)))
        else:
            raise ValueError('Unrecognized ldd output: ' + line)
    return dependencies


def check_dependencies(root, dependencies):
    root = Path(root).resolve()
    for name, path in dependencies:
        try:
            path.resolve().relative_to(root)
            continue
        except ValueError:
            pass
        if name not in SYSTEM_LIBRARIES:
            raise ValueError('Library requires host installation: ' + name + ' -> ' + str(path))
        # The SONAME alone cannot authorize an arbitrary external prefix.
        if not any(str(path.resolve()).startswith(prefix) for prefix in ('/lib/', '/lib64/', '/usr/lib/', '/usr/lib64/')):
            raise ValueError('System library outside OS library directories: ' + str(path))


def check_rpaths(root, owner, listing):
    root, owner = Path(root).resolve(), Path(owner).resolve()
    for value in re.findall(r'\((?:RPATH|RUNPATH)\).*?\[(.*?)\]', listing):
        for entry in value.split(':'):
            if not entry.startswith(('$ORIGIN', '${ORIGIN}')):
                raise ValueError('Non-relocatable Linux RPATH: ' + entry)
            expanded = entry.replace('${ORIGIN}', str(owner.parent)).replace('$ORIGIN', str(owner.parent))
            if '$' in expanded:
                raise ValueError('Unreviewed Linux loader token: ' + entry)
            try:
                Path(expanded).resolve().relative_to(root)
            except ValueError:
                raise ValueError('Linux RPATH escapes runtime: ' + entry)


def audit(root, kind, allow_manifest=False):
    if platform.system() != 'Linux' or platform.machine().lower() not in ('x86_64', 'amd64'):
        raise RuntimeError('Run this audit on Linux x86_64')
    root = Path(root).resolve()
    required = ['bin/python3'] if kind == 'python' else ['bin/cbuccaneer', 'bin/refmac5',
        'include/environ.def', 'include/default.def', 'lib/data/monomers',
        'lib/data/reference_structures/reference-1tqw.pdb', 'lib/data/reference_structures/reference-1tqw.mtz']
    for name in required:
        if not (root / name).exists():
            raise ValueError('Required runtime component missing: ' + name)
    if not shutil.which('readelf', path='/usr/bin:/bin') or not shutil.which('ldd', path='/usr/bin:/bin'):
        raise RuntimeError('readelf and ldd must already be available on the build host')
    records, elf_checks, versions = [], [], set()
    environment = {'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C'}
    manifest_name = 'python_manifest.json' if kind == 'python' else 'runtime_manifest.json'
    if (root / manifest_name).exists() and not allow_manifest:
        raise FileExistsError('Refusing to overwrite runtime manifest')
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            try:
                path.resolve(strict=True).relative_to(root)
            except (ValueError, OSError, RuntimeError):
                raise ValueError('Broken or escaping runtime symlink: ' + str(path))
        if not path.is_file():
            continue
        with path.open('rb') as handle:
            header = handle.read(20)
        if header.startswith(b'\x7fELF'):
            if len(header) < 20 or header[:6] != b'\x7fELF\x02\x01' or int.from_bytes(header[18:20], 'little') != 62:
                raise ValueError('Expected Linux x86_64 ELF: ' + str(path))
            dynamic = subprocess.run(['readelf', '-d', str(path)], env=environment,
                                     check=True, text=True, capture_output=True).stdout
            check_rpaths(root, path, dynamic)
            dependencies = []
            if '(NEEDED)' in dynamic:
                loaded = subprocess.run(['ldd', str(path)], env=environment,
                                        check=True, text=True, capture_output=True)
                dependencies = parse_ldd(loaded.stdout)
                check_dependencies(root, dependencies)
            symbols = subprocess.run(['readelf', '--version-info', str(path)], env=environment,
                                     check=True, text=True, capture_output=True).stdout
            versions.update(tuple(map(int, v.split('.'))) for v in re.findall(r'GLIBC_(\d+(?:\.\d+)+)', symbols))
            elf_checks.append({'relative_path': str(path.relative_to(root)),
                               'dependencies': [{'name': n, 'path': str(p)} for n, p in dependencies]})
        records.append({'relative_path': str(path.relative_to(root)), 'sha256': digest(path),
                        'bytes': path.stat().st_size})
    return {'platform': 'linux-x86_64', 'runtime_kind': kind, 'files': records,
            'elf_checks': elf_checks, 'glibc_build_host': platform.libc_ver(),
            'glibc_minimum_from_symbols': '.'.join(map(str, max(versions))) if versions else None,
            'public_redistribution_authorized': False,
            'scope': 'Staging inventory and ELF loader audit, not a complete application or scientific validation'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--kind', choices=['python', 'engines'], required=True)
    args = parser.parse_args(argv)
    report = audit(args.root, args.kind)
    name = 'python_manifest.json' if args.kind == 'python' else 'runtime_manifest.json'
    with (args.root / name).open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps({'manifest': str(args.root / name), 'files': len(report['files']),
                      'elf_files': len(report['elf_checks'])}, indent=2))


if __name__ == '__main__':
    main()
