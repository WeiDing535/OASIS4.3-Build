"""Scientific input validation and reference packaging contracts for desktop jobs."""
import hashlib
import json
from pathlib import Path
import numpy as np
import gemmi
import pytest
from oasis_py.reference_data import bundled_reference,read_reference,reference_info
from oasis_py.desktop_job import build_sad_command
from oasis_py.cli import main
from test_science import dataset


def job(tmp_path):
    data=dataset('P 1',extent=4)
    m=gemmi.Mtz(with_base=True);m.cell=data.cell;m.spacegroup=data.spacegroup;m.add_dataset('test')
    for label,kind in [('FP','F'),('SIGFP','Q'),('DANO','D'),('SIGDANO','Q')]:m.add_column(label,kind)
    m.set_data(np.column_stack([data.hkl,data.f,data.sigf,np.ones(len(data.f)),np.ones(len(data.f))]).astype('float32'))
    src=tmp_path/'input with space.mtz';m.write_to_file(str(src))
    from oasis_py.data import Atoms
    from oasis_py.substructure import write_atoms
    heavy=tmp_path/'heavy.pdb';write_atoms(heavy,Atoms([[.2,.3,.4]],['Se'],[1.],[20.]),data)
    return dict(input=str(src),heavy=str(heavy),output=str(tmp_path/'output.mtz'),
        fpp='Se=3.8',composition='C:100,N:30,O:40,Se:1',solvent='.5',
        cycles=5,dm_cycles=5,labels={})


@pytest.mark.reference_data
def test_bundled_reference_hash_provenance_and_columns():
    m=read_reference();info=reference_info()
    assert m.nreflections==328518 and info['pdb_id']=='1TQW'
    assert info['coefficients_changed'] is False
    assert hashlib.sha256(bundled_reference().read_bytes()).hexdigest()==info['mtz_sha256']
    assert np.isfinite(info['resolution_limit'])


def test_reference_validation_does_not_accept_experimental_amplitudes_as_map(tmp_path):
    spec=job(tmp_path)
    with pytest.raises(ValueError,match='amplitude.*phase'):read_reference(spec['input'])


def test_desktop_command_is_shell_free_and_uses_bundled_reference(tmp_path):
    spec=job(tmp_path);cmd=build_sad_command(spec)
    assert cmd[cmd.index('phase')+1]==spec['input']
    assert cmd[cmd.index('--reference-mtz')+1]==str(bundled_reference())
    assert cmd[cmd.index('--dm-engine')+1]=='mlhl-best'
    assert '--composition' in cmd and '--residues' not in cmd
    assert not any('parrot' in x.lower() for x in cmd)


@pytest.mark.parametrize('key,value',[
    ('solvent','1.1'),('solvent','nan'),('fpp','Se=-1'),('fpp','Hg=3'),
    ('composition','C:-1'),('cycles',0),('dm_cycles',101),('dmin','nan'),
])
def test_desktop_rejects_invalid_scientific_inputs_before_execution(tmp_path,key,value):
    spec=job(tmp_path);spec[key]=value
    with pytest.raises((ValueError,KeyError)):build_sad_command(spec)
    assert not Path(spec['output']).exists()


def test_desktop_rejects_existing_metadata_and_preserves_it(tmp_path):
    spec=job(tmp_path);sidecar=Path(spec['output']).with_suffix('.json');sidecar.write_text('keep')
    with pytest.raises(ValueError,match='文件已存在'):build_sad_command(spec)
    assert sidecar.read_text()=='keep'


def test_pplus_uses_validated_legacy_dm_and_dm_can_be_disabled(tmp_path):
    spec=job(tmp_path);spec['method']='pplus';cmd=build_sad_command(spec)
    assert cmd[cmd.index('--dm-engine')+1]=='legacy'
    spec['density_modification']=False;spec['solvent']='';cmd=build_sad_command(spec)
    assert '--native-dm' not in cmd


def test_reference_info_cli_runs_without_target_data(capsys):
    assert main(['reference-info'])==0
    assert json.loads(capsys.readouterr().out)['coefficients_changed'] is False
