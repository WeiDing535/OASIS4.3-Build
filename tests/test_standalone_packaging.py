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
