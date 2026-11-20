"""Source packages and portable launchers include the integrated phase entry."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import pytest

ROOT=Path(__file__).resolve().parents[1]


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_focused_source_archive_ships_the_desktop_window_and_refuses_overwrite(tmp_path):
    package=load('oasis43_source_package',ROOT/'package_oasis.py')
    target=tmp_path/'source.zip';package.main(['--output',str(target)])
    before=target.read_bytes()
    with zipfile.ZipFile(target) as archive:
        prefix='OASIS4.3-focused-development/'
        manifest=json.loads(archive.read(prefix+'PACKAGE_MANIFEST.json'))
        records={item['path']:item['sha256'] for item in manifest['files']}
        for name in ['run_desktop.py','run_sad_to_model.py','src/phasebuild/desktop.py',
                     'src/phasebuild/jobs.py','src/phasebuild/oasis_workflow.py']:
            assert name in records
        # The retired browser front end and its assets are not shipped at all.
        for name in records:
            assert not name.endswith(('run_gui.py','package_gui.py','phasebuild/gui.py','phasebuild/gui_workflow.py'))
            assert 'gui_assets' not in name and 'oasis_assets' not in name
        for name,digest in records.items():
            assert hashlib.sha256(archive.read(prefix+name)).hexdigest()==digest
        assert not any('/runtime/bin/' in n or '/cases/' in n for n in archive.namelist())
        assert archive.getinfo(prefix+'OASIS4.3.command').external_attr>>16&0o111
    with pytest.raises(FileExistsError):package.main(['--output',str(target)])
    assert target.read_bytes()==before


def test_portable_launcher_dispatches_both_hl_modes_and_records_current_sources(tmp_path):
    builder=load('phasebuild_app_assembler',ROOT/'packaging/assemble_bundle.py')
    bundle=tmp_path/'bundle with spaces';builder.main(['--output',str(bundle)])
    manifest=json.loads((bundle/'app/app_manifest.json').read_text())
    files={f['relative_path']:f for f in manifest['files']}
    assert 'run_build_from_phases.py' in files
    for name,record in files.items():
        assert hashlib.sha256((bundle/'app'/name).read_bytes()).hexdigest()==record['sha256']
        assert (bundle/'app'/name).read_bytes()==(ROOT/name).read_bytes()
    # Exercise the real launcher/parser only; no engine or structure job runs.
    (bundle/'python/bin').mkdir(parents=True)
    (bundle/'python/bin/python3').symlink_to(sys.executable)
    for mode in ['build-hl','hybrid-hl']:
        process=subprocess.run([str(bundle/'phasebuild'),mode,'--help'],capture_output=True,text=True)
        assert process.returncode==0,process.stderr
        assert '--free-value' in process.stdout
        assert ('--model' in process.stdout)==(mode=='hybrid-hl')
    before=(bundle/'phasebuild').read_bytes()
    with pytest.raises(FileExistsError):builder.main(['--output',str(bundle)])
    assert (bundle/'phasebuild').read_bytes()==before


def test_linux_launcher_runs_from_space_path_with_host_runtime_disabled(tmp_path):
    import os
    builder = load('phasebuild_linux_assembler', ROOT / 'packaging/assemble_bundle.py')
    bundle = tmp_path / 'Linux bundle with spaces'
    builder.main(['--platform', 'linux-x86_64', '--output', str(bundle)])
    assert (bundle / 'OASIS4.3-Build.sh').is_file()
    (bundle / 'python/bin').mkdir(parents=True)
    (bundle / 'python/bin/python3').symlink_to(sys.executable)
    environment = dict(os.environ, CCP4='/missing/host/ccp4', PHASEBUILD_RUNTIME='/missing/host/runtime')
    process = subprocess.run([str(bundle / 'phasebuild'), '--help'], capture_output=True, text=True, env=environment)
    assert process.returncode == 0 and 'linux-x86_64' in process.stdout
    process = subprocess.run([str(bundle / 'OASIS4.3-Build.sh'), '--help'], capture_output=True, text=True, env=environment)
    assert process.returncode == 0 and '--exit-on-finish' in process.stdout
