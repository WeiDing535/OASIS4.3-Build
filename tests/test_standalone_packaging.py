"""Reject damaged or host-dependent runtime copies before packaging."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

SPEC = importlib.util.spec_from_file_location('standalone_builder', Path(__file__).resolve().parents[1] / 'packaging/build_standalone.py')
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)


def fixture_runtime(root):
    root.mkdir()
    (root / 'component').write_bytes(b'known-runtime-bytes')
    manifest = {'files': [{'relative_path': 'component', 'sha256': hashlib.sha256(b'known-runtime-bytes').hexdigest()}]}
    (root / 'manifest.json').write_text(json.dumps(manifest))
    return root


def test_runtime_verifier_detects_corruption(tmp_path):
    root = fixture_runtime(tmp_path / 'runtime')
    assert BUILDER.verify_runtime(root, 'manifest.json')['files_checked'] == 1
    (root / 'component').write_bytes(b'changed')
    with pytest.raises(ValueError, match='hash differs'):
        BUILDER.verify_runtime(root, 'manifest.json')


def test_runtime_verifier_refuses_unlisted_escape_link(tmp_path):
    root = fixture_runtime(tmp_path / 'runtime')
    (tmp_path / 'host-file').write_text('host dependency')
    (root / 'escape').symlink_to(tmp_path / 'host-file')
    with pytest.raises(ValueError, match='symlink escapes'):
        BUILDER.verify_runtime(root, 'manifest.json')


def test_existing_bundle_is_never_overwritten(tmp_path):
    out = tmp_path / 'existing';out.mkdir()
    with pytest.raises(FileExistsError, match='overwrite bundle'):
        BUILDER.main(['--python-source', str(tmp_path / 'none'), '--runtime-source', str(tmp_path / 'none'), '--output', str(out)])
    assert list(out.iterdir()) == []


def test_linux_platform_and_engine_architecture_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(BUILDER.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(BUILDER.platform, 'machine', lambda: 'x86_64')
    assert BUILDER.target_platform() == 'linux-x86_64'
    root = tmp_path / 'runtime';(root / 'bin').mkdir(parents=True)
    header = bytearray(20);header[:6] = b'\x7fELF\x02\x01';header[18:20] = (62).to_bytes(2, 'little')
    for name in ['cbuccaneer', 'refmac5']:(root / 'bin' / name).write_bytes(header)
    BUILDER.verify_linux_engines(root)
    (root / 'bin/refmac5').write_bytes(b'macOS binary')
    with pytest.raises(ValueError, match='Linux x86_64 ELF'):
        BUILDER.verify_linux_engines(root)


def test_linux_tar_preserves_executable_and_relative_symlink(tmp_path):
    import tarfile
    root = tmp_path / 'bundle with spaces';root.mkdir()
    (root / 'launcher').write_text('#!/bin/sh\nexit 0\n');(root / 'launcher').chmod(0o755)
    (root / 'alias').symlink_to('launcher')
    target = tmp_path / 'bundle.tar.gz';BUILDER.write_tar(root, target)
    with tarfile.open(target) as archive:
        assert archive.getmember(root.name + '/launcher').mode & 0o111
        link = archive.getmember(root.name + '/alias')
        assert link.issym() and link.linkname == 'launcher'
    with pytest.raises(FileExistsError):BUILDER.write_tar(root, target)


def test_doctor_cannot_silently_use_host_python(tmp_path):
    root = tmp_path / 'bundle'
    report = {'python': str(root / 'python/bin/python3'), 'python_prefix': str(root / 'python'),
              'packages': {'gemmi': {'path': str(root / 'python/lib/gemmi.so')}},
              'runtime': {'root': str(root / 'runtime')}}
    BUILDER.verify_doctor_paths(report, root)
    report['packages']['gemmi']['path'] = '/host/python/gemmi.so'
    with pytest.raises(ValueError, match='outside its prefix'):
        BUILDER.verify_doctor_paths(report, root)
