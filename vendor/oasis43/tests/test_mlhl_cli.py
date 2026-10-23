import numpy as np
import pytest
import gemmi
from oasis_py.cli import main
from oasis_py.data import write_mtz
from test_science import dataset


@pytest.mark.parametrize('extra',[
    [],
    ['--native-dm','--solvent-fraction','.5','--reference-mtz','reference.mtz','--sad-method','pplus'],
    ['--native-dm','--solvent-fraction','.5','--reference-mtz','reference.mtz','--dm-anisotropy'],
])
def test_unvalidated_mlhl_combinations_fail_before_phasing(tmp_path,capsys,extra):
    data=dataset();data.dano=np.zeros(len(data.f));data.sigdano=np.ones(len(data.f))
    path=tmp_path/'input.mtz';write_mtz(path,data,np.zeros(len(data.f)),np.ones(len(data.f)))
    m=gemmi.read_mtz_file(str(path));base=m.array.copy()
    m.add_column('DANO','D');m.add_column('SIGDANO','Q')
    m.set_data(np.column_stack([base,np.zeros(len(base)),np.ones(len(base))]).astype('float32'))
    m.write_to_file(str(path))
    code=main(['phase',str(path),'--mode','sad','--sad-method','likelihood',
        '--dm-engine','mlhl-best','--out',str(tmp_path/'out.mtz'),*extra])
    assert code==2
    assert 'mlhl-best' in capsys.readouterr().err.lower()
    assert not (tmp_path/'out.mtz').exists()


def test_mlhl_public_api_rejects_invalid_cycle_and_volume_settings():
    from oasis_py.density_mlhl import density_modify_mlhl
    with pytest.raises(ValueError,match='cycles or solvent fraction'):
        density_modify_mlhl(None,None,.5,None,cycles=0)
    with pytest.raises(ValueError,match='cycles or solvent fraction'):
        density_modify_mlhl(None,None,1.1,None,cycles=1)
