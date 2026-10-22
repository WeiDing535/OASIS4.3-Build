"""Locate a relocatable external-engine runtime without changing science settings."""
import json
import os
import tempfile
from pathlib import Path


def _candidates():
    source = Path(__file__).resolve()
    # Source checkout: <app>/src/phasebuild. Portable app: <bundle>/app/src/phasebuild.
    yield source.parents[2] / 'runtime'
    yield source.parents[3] / 'runtime'
    if os.environ.get('PHASEBUILD_STANDALONE') != '1':
        if os.environ.get('CCP4'):
            yield Path(os.environ['CCP4'])
        # Compatibility discovery for existing local installations only. Portable
        # launchers disable this fallback, so missing bundled files fail clearly.
        yield Path('/Applications/ccp4-9')


def default_runtime():
    explicit = os.environ.get('PHASEBUILD_RUNTIME')
    if explicit:
        return str(Path(explicit).expanduser().resolve())
    for path in _candidates():
        if (path / 'bin').is_dir():
            return str(path.resolve())
    return ''


def resolve_runtime(explicit=None, required=()):
    value = explicit if explicit is not None else default_runtime()
    if not value:
        raise ValueError('External engine runtime unavailable; set --runtime or PHASEBUILD_RUNTIME')
    root = Path(value).expanduser().resolve()
    if not root.is_dir():
        raise ValueError('Runtime directory does not exist: ' + str(root))
    for name in required:
        if name not in ('refmac5', 'cbuccaneer'):
            raise ValueError('Unknown runtime executable: ' + name)
        binary = root / 'bin' / name
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError('Runtime executable unavailable: ' + str(binary))
    return root


def runtime_environment(root, inherited=None, scratch=None):
    root = resolve_runtime(root)
    environment = dict(os.environ if inherited is None else inherited)
    for name in list(environment):
        if name.startswith(('DYLD_', 'CCP4')) or name in (
                'CLIBD', 'CLIBD_MON', 'CBIN', 'CEXAM', 'CINCL', 'CLIB',
                'PYTHONHOME', 'PYTHONPATH', 'BINSORT_SCR'):
            environment.pop(name, None)
    # CCP4's Fortran startup requires an existing writable scratch directory.
    # Replace a host installation's setting, and keep per-job scratch outside
    # the bundle so the application also works from a read-only installation.
    scratch = Path(scratch if scratch is not None else tempfile.gettempdir()).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    if not scratch.is_dir() or not os.access(scratch, os.W_OK):
        raise ValueError('Runtime scratch directory is not writable: ' + str(scratch))
    environment.update(CCP4=str(root), CINCL=str(root / 'include'), CCP4_SCR=str(scratch), CLIBD=str(root / 'lib/data'),
                       CLIBD_MON=str(root / 'lib/data/monomers') + '/')
    return environment


def runtime_info(root=None):
    root = resolve_runtime(root)
    manifest = root / 'runtime_manifest.json'
    result = {'root': str(root), 'external_engines': ['Buccaneer', 'REFMAC'],
              'bundled_manifest': str(manifest) if manifest.is_file() else None,
              'standalone_mode': os.environ.get('PHASEBUILD_STANDALONE') == '1'}
    if manifest.is_file():
        metadata = json.loads(manifest.read_text())
        result['distribution_status'] = metadata.get('distribution_status', 'local-evaluation')
    return result
