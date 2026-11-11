"""Local portable-builder safety and runtime contract."""
import importlib.util
import json
from pathlib import Path
import pytest

SPEC=importlib.util.spec_from_file_location('free_packager',Path(__file__).resolve().parents[1]/'packaging/build_installation_free.py')
BUILDER=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(BUILDER)


def test_existing_destination_is_refused_without_touching_it(tmp_path):
    out=tmp_path/'existing';out.mkdir();marker=out/'keep.txt';marker.write_text('preserve')
    with pytest.raises(FileExistsError):BUILDER.main(['--python-source','missing','--runtime-source','missing','--output',str(out)])
    assert marker.read_text()=='preserve' and list(out.iterdir())==[marker]


def test_copy_manifest_rejects_escaping_symlink(tmp_path):
    root=tmp_path/'bundle';root.mkdir();outside=tmp_path/'outside';outside.write_text('private')
    (root/'link').symlink_to(outside)
    with pytest.raises(ValueError,match='escaping'):BUILDER.manifest(root)


def test_launcher_enforces_no_suite_discovery_and_uses_relative_components():
    script=BUILDER.LAUNCHER
    assert 'OASIS_INSTALLATION_FREE=1' in script and 'PHASEBUILD_STANDALONE=1' in script
    assert '"$base/python/bin/python3"' in script and 'PHASEBUILD_RUNTIME="$base/runtime"' in script
    assert '/Applications/' not in script and 'unset PHENIX' in script


def test_reused_runtime_manifest_does_not_hash_its_previous_self(tmp_path):
    root=tmp_path/'python';root.mkdir()
    (root/'python_manifest.json').write_text('{"files":[]}')
    (root/'library.so').write_bytes(b'unchanged library')
    rows=BUILDER.manifest(root,exclude={'python_manifest.json'})
    assert [r['relative_path'] for r in rows]==['library.so']
    assert rows[0]['sha256']==BUILDER.sha(root/'library.so')


def test_bundled_engine_list_comes_from_the_reused_runtime_manifest(tmp_path):
    root=tmp_path/'bundle';(root/'runtime').mkdir(parents=True)
    rows=[{'relative_path':'bin/'+name} for name in ('cbuccaneer','refmac5','cparrot')]
    (root/'runtime/runtime_manifest.json').write_text(json.dumps({'dependency_graph':rows}))
    assert BUILDER.bundled_engines(root)==['cparrot','cbuccaneer','refmac5']
    assert BUILDER.ENGINE_LABELS['cparrot']=='Parrot'
    rows.append({'relative_path':'bin/phaser'})
    (root/'runtime/runtime_manifest.json').write_text(json.dumps({'dependency_graph':rows}))
    assert BUILDER.bundled_engines(root)==['phaser','cparrot','cbuccaneer','refmac5']
