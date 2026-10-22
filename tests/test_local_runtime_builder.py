import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest

PATH = Path(__file__).resolve().parents[1] / 'packaging/build_local_runtime.py'
SPEC = importlib.util.spec_from_file_location('phasebuild_local_runtime_builder', PATH)
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)


def test_otool_install_id_is_not_a_dependency():
    identity = '/obsolete/build/lib/libclipper-core.2.dylib'
    def run(command, **kwargs):
        if command[1] == '-D':
            return SimpleNamespace(stdout='file:\n' + identity + '\n')
        return SimpleNamespace(stdout='file:\n\t' + identity + ' (compatibility version3)\n' +
            '\t@loader_path/libfftw.2.dylib (compatibility version3)\n' +
            '\t/usr/lib/libSystem.B.dylib (compatibility version1)\n')
    names, evidence = BUILDER.dependency_names(Path('/local/libclipper-core.2.dylib'), run=run)
    assert names == ['@loader_path/libfftw.2.dylib', '/usr/lib/libSystem.B.dylib']
    assert evidence['dylib_install_ids_excluded'] == [identity]


def test_resolver_preserves_loader_relative_layout_and_rejects_host_fallback(tmp_path):
    source = tmp_path / 'source'
    (source / 'lib').mkdir(parents=True)
    (source / 'lib/libfixture.dylib').write_bytes(b'fixture')
    assert BUILDER.resolve_dependency('@loader_path/../lib/libfixture.dylib', source / 'bin/tool', source) == source / 'lib/libfixture.dylib'
    assert BUILDER.resolve_dependency('/usr/lib/libSystem.B.dylib', source / 'bin/tool', source) is None
    with pytest.raises(RuntimeError, match='absolute dependency'):
        BUILDER.resolve_dependency('/host/libfixture.dylib', source / 'bin/tool', source)


def test_copy_verifies_all_bytes_and_dereferences_file_symlinks(tmp_path):
    source = tmp_path / 'original'
    source.write_bytes(bytes(range(256)) * 16)
    alias = tmp_path / 'alias'
    alias.symlink_to(source.name)
    bundle = tmp_path / 'runtime'
    record = BUILDER.copy_verified(alias, bundle / 'lib/alias', 'dynamic_library', bundle, 'fixture_license')
    assert record['relative_path'] == record['path'] == 'lib/alias'
    assert record['sha256'] == record['source_sha256']
    assert record['source_file_symlink'] == 'original'
    assert not (bundle / 'lib/alias').is_symlink()
    assert source.read_bytes() == (bundle / 'lib/alias').read_bytes()


def test_ccp4fyp_definitions_are_required_runtime_files(tmp_path):
    source, bundle = tmp_path / 'source', tmp_path / 'runtime'
    (source / 'include').mkdir(parents=True)
    (source / 'include/environ.def').write_text('HKLIN=in.mtz\n')
    with pytest.raises(RuntimeError, match='default.def'):
        BUILDER.startup_definition_plan(source, bundle)
    (source / 'include/default.def').write_text('ATOMSF=atomsf.lib\n')
    plan = BUILDER.startup_definition_plan(source, bundle)
    assert [item[1].relative_to(bundle).as_posix() for item in plan] == [
        'include/environ.def', 'include/default.def']
    records = [BUILDER.copy_verified(original, destination, role, bundle, group)
               for original, destination, role, group in plan]
    assert all(record['sha256'] == record['source_sha256'] for record in records)


def test_startup_probe_uses_application_environment_without_host_scratch(tmp_path, monkeypatch):
    runtime = tmp_path / 'runtime'
    runtime.mkdir()
    monkeypatch.setenv('CCP4_SCR', '/host-only/scratch')
    monkeypatch.setenv('HOME', '/host-only/home')
    scratch = tmp_path / 'job/scratch'
    environment, source = BUILDER.probe_runtime_environment(runtime, scratch)
    assert environment['CCP4_SCR'] == str(scratch.resolve())
    assert scratch.is_dir()
    assert 'HOME' not in environment
    assert environment['CINCL'] == str(runtime / 'include')
    assert source['sha256'] == BUILDER.sha256(PATH.parents[1] / 'src/phasebuild/runtime.py')
