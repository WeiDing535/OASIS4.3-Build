#!/usr/bin/env python3
"""Remove components the SAD pipeline never loads from a built bundle.

Every rule here is a directory or glob that belongs to a GUI, plotting or
validation tool that the bundled engines (phaser, refmac5, cbuccaneer, cparrot,
RESOLVE through oasis-density-python) never import. Anything the pipeline needs
is kept, and package_standalone re-runs its smoke tests after pruning.
"""
import hashlib, json, os, shutil, sys
from pathlib import Path

SITE_PACKAGES_27 = 'runtime/density/conda_base/lib/python2.7/site-packages'
PYTHON_312 = 'python/lib/python3.12'

# (path or glob relative to the bundle root, why)
DROP = [
    # GUI and plotting stacks inside the RESOLVE python 2.7 environment.
    (SITE_PACKAGES_27 + '/bokeh', 'plotting'),
    (SITE_PACKAGES_27 + '/PySide2', 'Qt bindings'),
    (SITE_PACKAGES_27 + '/wx*', 'GUI toolkit'),
    (SITE_PACKAGES_27 + '/matplotlib*', 'plotting'),
    (SITE_PACKAGES_27 + '/pandas', 'data frames'),
    (SITE_PACKAGES_27 + '/sklearn', 'machine learning'),
    (SITE_PACKAGES_27 + '/babel', 'localisation'),
    (SITE_PACKAGES_27 + '/IPython', 'interactive shell'),
    (SITE_PACKAGES_27 + '/notebook', 'browser notebook'),
    (SITE_PACKAGES_27 + '/jupyter*', 'browser notebook'),
    (SITE_PACKAGES_27 + '/ensurepip', 'bundled pip'),
    (SITE_PACKAGES_27 + '/plat-mac', 'legacy platform'),
    ('runtime/density/conda_base/lib/libQt5*', 'Qt runtime for the dropped GUI stack'),
    ('runtime/density/conda_base/lib/libicu*', 'Qt dependency'),
    ('runtime/density/conda_base/plugins', 'Qt plugins'),
    ('runtime/density/conda_base/translations', 'Qt translations'),
    # Installer-only packages in the bundled python 3.12.
    (PYTHON_312 + '/site-packages/pip', 'installer'),
    (PYTHON_312 + '/site-packages/setuptools', 'installer'),
    (PYTHON_312 + '/site-packages/pygments', 'syntax highlighting'),
    (PYTHON_312 + '/site-packages/_pytest', 'test runner'),
    (PYTHON_312 + '/site-packages/pytest', 'test runner'),
    (PYTHON_312 + '/site-packages/pluggy', 'test runner'),
    (PYTHON_312 + '/site-packages/iniconfig', 'test runner'),
    (PYTHON_312 + '/idlelib', 'IDLE'),
    (PYTHON_312 + '/ensurepip', 'bundled pip'),
    (PYTHON_312 + '/lib2to3', '2to3 converter'),
    (PYTHON_312 + '/turtledemo', 'turtle demos'),
    (PYTHON_312 + '/pydoc_data', 'pydoc text'),
    (PYTHON_312 + '/test', 'test suite'),
    # Sugar and nucleotide validation data, used by nucleofind/privateer only.
    ('runtime/lib/data/nucleofind_models', 'nucleofind models'),
    ('runtime/lib/data/privateer_data', 'privateer data'),
]

def _expand(root, pattern):
    root = Path(root)
    if any(ch in pattern for ch in '*?['):
        return sorted(root.glob(pattern))
    path = root / pattern
    return [path] if path.exists() else []

def prune(bundle, keep_reference=True):
    root = Path(bundle).expanduser().resolve()
    removed = []
    for pattern, why in DROP:
        for path in _expand(root, pattern):
            if not path.exists():
                continue
            size = sum(p.stat().st_size for p in path.rglob('*') if p.is_file()) if path.is_dir() else path.stat().st_size
            shutil.rmtree(path) if path.is_dir() else path.unlink()
            removed.append(dict(path=str(path.relative_to(root)), bytes=size, why=why))
    # RESOLVE reference structures: the GUI default points at reference-1tqw.mtz.
    refs = root / 'runtime/lib/data/reference_structures'
    if keep_reference and refs.is_dir():
        for path in sorted(refs.iterdir()):
            if path.name.startswith('reference-1tqw.'):
                continue
            size = sum(p.stat().st_size for p in path.rglob('*') if p.is_file()) if path.is_dir() else path.stat().st_size
            shutil.rmtree(path) if path.is_dir() else path.unlink()
            removed.append(dict(path=str(path.relative_to(root)), bytes=size, why='unused RESOLVE reference'))
    # Byte-code caches are regenerated on first import.
    caches = [p for p in root.rglob('__pycache__') if p.is_dir()]
    for path in caches:
        size = sum(p.stat().st_size for p in path.rglob('*') if p.is_file())
        shutil.rmtree(path)
        removed.append(dict(path=str(path.relative_to(root)), bytes=size, why='byte-code cache'))
    return removed

# Payloads that ship with the RESOLVE environment but no bundled program loads.
DROP_LATE = [
    ('runtime/density/conda_base/lib/libLLVM-9.dylib', 'LLVM (no llvmlite/numba present)'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/sphinx', 'documentation builder'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/OpenGL', 'OpenGL bindings'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/py2app', 'app bundler'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/pip', 'installer'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/pygments', 'syntax highlighting'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/google', 'cloud client stubs'),
    ('runtime/density/conda_base/lib/python2.7/site-packages/setuptools', 'installer'),
    # segments/ is NOT droppable: RESOLVE reads rho.list, Richardson_rho.dat and
    # side_chains_rho.dat from it whenever the pipeline has to discriminate the
    # hand, so dropping it made every auto-hand run die at the P+ stage.
    ('runtime/density/modules/solve_resolve/ext_ref_files/segments_em', 'electron-microscopy segment library'),
    ('runtime/density/modules/solve_resolve/ext_ref_files/segments_rich_small', 'alternative RESOLVE segment library'),
    ('runtime/density/modules/solve_resolve/ext_ref_files/patterns/patterns_txt', 'text dump of the pattern library'),
]


def prune_late(root):
    removed = []
    for pattern, why in DROP_LATE:
        for path in _expand(root, pattern):
            if not path.exists():
                continue
            size = sum(p.stat().st_size for p in path.rglob('*') if p.is_file()) if path.is_dir() else path.stat().st_size
            shutil.rmtree(path) if path.is_dir() else path.unlink()
            removed.append(dict(path=str(path.relative_to(root)), bytes=size, why=why))
    return removed


def dedupe(root, min_size=200_000):
    """Replace byte-identical payloads with relative symlinks.

    The conda environment ships eight identical BLAS/LAPACK libraries and repeats
    several data files; linking them keeps every embedded path working while
    removing the copies.
    """
    root = Path(root)
    groups = {}
    for dirpath, dirs, files in os.walk(root):
        for name in files:
            path = Path(dirpath) / name
            if path.is_symlink():
                continue
            try:
                size = path.stat().st_size
            except OSError:
                continue
            if size >= min_size:
                groups.setdefault(size, []).append(path)
    saved = 0
    linked = 0
    for size, paths in groups.items():
        if len(paths) < 2:
            continue
        buckets = {}
        for path in paths:
            buckets.setdefault(_digest(path), []).append(path)
        for digest, same in buckets.items():
            if len(same) < 2:
                continue
            canonical = same[0]
            for duplicate in same[1:]:
                relative = os.path.relpath(canonical, duplicate.parent)
                duplicate.unlink()
                os.symlink(relative, duplicate)
                saved += size
                linked += 1
    return saved, linked


# Monomer dictionaries for the residues a SAD phasing run actually builds, plus
# the common ions and cryo additives found in protein crystals.
MONOMER_KEEP = set('''
ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
MSE UNK HOH DOD
NA CL K MG CA ZN MN FE CU NI CD HG CO SR CS BR IOD F
SO4 PO4 NO3 NH4 ACT FMT ACE EDO GOL PEG PG4 MPD IPA CIT TRS MES BME DTT
ADP ATP GDP GTP AMP GMP FAD NAD NAP HEM SAH SAM
BPB BOG LDA DDM OLC PLM MYR STE
'''.split())


def prune_monomers(root):
    """Keep only protein/water/ion dictionaries; the shard index stays intact."""
    base = Path(root) / 'runtime/lib/data/monomers'
    if not base.is_dir():
        return []
    removed = []
    for shard in sorted(base.iterdir()):
        if not shard.is_dir() or len(shard.name) != 1:
            continue
        for path in sorted(shard.glob('*.cif')):
            if path.stem.upper() in MONOMER_KEEP:
                continue
            size = path.stat().st_size
            path.unlink()
            removed.append(dict(path=str(path.relative_to(root)), bytes=size,
                                why='monomer dictionary not used by SAD building'))
    return removed


def _digest(path):
    handle = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            handle.update(block)
    return handle.hexdigest()


def refresh_manifest(root, name='runtime_manifest.json'):
    """Drop manifest rows for files the prune removed; kept files are unchanged."""
    root = Path(root)
    manifest = root / name
    if not manifest.is_file():
        return 0, 0
    data = json.loads(manifest.read_text())
    kept = [row for row in data['files'] if (root / row['relative_path']).is_file()]
    dropped = len(data['files']) - len(kept)
    if dropped:
        data['files'] = kept
        manifest.write_text(json.dumps(data, indent=1) + '\n')
        checksum = manifest.with_name(name.replace('.json', '.sha256'))
        if checksum.is_file() or name == 'runtime_manifest.json':
            Path(str(manifest).replace('.json', '.sha256')).write_text(
                _digest(manifest) + '  ' + name + '\n')
    return len(kept), dropped


def main():
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle')
    parser.add_argument('--json')
    parser.add_argument('--keep-all-monomers', action='store_true',
                        help='retain every ligand dictionary (adds ~670 MB)')
    parser.add_argument('--keep-optional', action='store_true',
                        help='retain GUI/plotting payloads and keep duplicate files')
    args = parser.parse_args()
    removed = prune(args.bundle)
    if not args.keep_all_monomers:
        removed += prune_monomers(args.bundle)
    if not args.keep_optional:
        removed += prune_late(args.bundle)
    saved, linked = (0, 0) if args.keep_optional else dedupe(args.bundle)
    total = sum(item['bytes'] for item in removed) + saved
    print('removed %d items, %.1f MB (dedupe: %d links, %.1f MB)' % (
        len(removed), sum(item['bytes'] for item in removed) / 1e6, linked, saved / 1e6))
    for name in ('runtime/runtime_manifest.json', 'python/python_manifest.json'):
        kept, dropped = refresh_manifest(Path(args.bundle) / Path(name).parent, Path(name).name)
        if kept or dropped:
            print('%s: kept %d rows, dropped %d' % (name, kept, dropped))
    if args.json:
        Path(args.json).write_text(json.dumps(removed, indent=2) + '\n')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
