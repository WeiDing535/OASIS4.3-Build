import numpy as np
import pytest
import gemmi
from oasis_py.cli import main
from oasis_py.data import Reflections,Atoms,parse_composition,read_atoms,write_mtz,read_mtz
from test_science import dataset


def test_reject_duplicate_hkl_and_invalid_sigmas():
    d=dataset()
    with pytest.raises(ValueError,match='Duplicate'):
        Reflections([[1,0,0],[1,0,0]],[2,3],[1,1],d.cell,d.spacegroup)
    with pytest.raises(ValueError,match='positive amplitude uncertainties'):
        Reflections([[1,0,0]],[2],[0],d.cell,d.spacegroup)


def test_composition_validates_and_accumulates_elements():
    assert parse_composition('C:10,C:5,S:1')=={'C':15,'S':1}
    with pytest.raises(ValueError):parse_composition('X:3')
    with pytest.raises(ValueError):parse_composition('C:-1')
    with pytest.raises(ValueError):parse_composition()


def test_explicit_free_exclusion(tmp_path):
    d=dataset();path=tmp_path/'in.mtz'
    write_mtz(path,d,np.zeros(len(d.f)),np.ones(len(d.f)))
    mtz=gemmi.read_mtz_file(str(path));base=mtz.array.copy()
    mtz.add_column('FREE','I');flags=np.arange(len(base))%10
    mtz.set_data(np.column_stack([base,flags]).astype('float32'));mtz.write_to_file(str(path))
    kept,meta=read_mtz(path,free_label='FREE',free_value=0)
    assert len(kept.f)==int(np.sum(flags!=0))
    with pytest.raises(ValueError,match='explicitly'):read_mtz(path,free_label='FREE')


def test_amplitude_reader_rejects_intensities(tmp_path):
    d=dataset();m=gemmi.Mtz(with_base=True);m.spacegroup=d.spacegroup;m.set_cell_for_all(d.cell)
    m.add_dataset('I');m.add_column('I','J');m.add_column('SIGI','Q')
    m.set_data(np.column_stack([d.hkl,d.f,d.sigf]).astype('float32'));p=tmp_path/'i.mtz';m.write_to_file(str(p))
    with pytest.raises(ValueError,match='amplitudes, not intensities'):
        read_mtz(p,labels={'f':'I','sigf':'SIGI'})


def test_cli_never_overwrites_input_or_existing_output(tmp_path,capsys):
    d=dataset();p=tmp_path/'original.mtz';write_mtz(p,d,np.zeros(len(d.f)),np.ones(len(d.f)))
    before=p.read_bytes()
    result=main(['phase',str(p),'--out',str(p),'--residues','50'])
    assert result==2 and p.read_bytes()==before
    assert 'Output exists' in capsys.readouterr().err


def test_model_cell_mismatch_is_not_silently_ignored(tmp_path):
    from oasis_py.substructure import write_atoms
    d=dataset();a=Atoms([[.2,.3,.4]],['S'],[1],[10]);p=tmp_path/'s.pdb';write_atoms(p,a,d)
    other=dataset();other.cell=gemmi.UnitCell(30,30,30,90,90,90)
    with pytest.raises(ValueError,match='unit cells differ'):read_atoms(p,other)


def test_missing_sad_columns_have_actionable_error(tmp_path):
    d=dataset();p=tmp_path/'native.mtz';write_mtz(p,d,np.zeros(len(d.f)),np.ones(len(d.f)))
    with pytest.raises(ValueError,match='Missing dano'):read_mtz(p,'sad')
