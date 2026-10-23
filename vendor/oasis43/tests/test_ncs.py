import numpy as np
import gemmi,pytest
from oasis_py.native_density import metric_filter,smooth
from oasis_py.ncs import prepare_ncs,apply_ncs,resample,validate_operators


def test_known_ncs_reduces_independent_noise_in_duplicate_copies():
    shape=(32,32,32);rng=np.random.default_rng(237)
    cell=gemmi.UnitCell(64,64,64,90,90,90);kernel=metric_filter(cell,shape,4.)
    half=smooth(rng.normal(size=shape),metric_filter(cell,shape,1.))[:16]
    clean=np.concatenate([half,half],axis=0);rho=clean+rng.normal(size=shape)*clean.std()*.2
    grid=gemmi.FloatGrid(*shape);grid.unit_cell=cell;grid.spacegroup=gemmi.SpaceGroup('P 1');grid.array[:]=rho
    ops=[{'rotation':np.eye(3).tolist(),'translation':[32.,0,0]}]
    plans,diagnostics=prepare_ncs(grid,np.ones(shape),kernel,ops)
    averaged=apply_ncs(rho,cell,plans)
    assert np.mean((averaged-clean)**2)<.75*np.mean((rho-clean)**2)
    assert max(v['active_fraction'] for v in diagnostics)>.5


def test_wrong_ncs_is_rejected_by_local_density_correlation():
    shape=(32,32,32);rng=np.random.default_rng(238);cell=gemmi.UnitCell(64,64,64,90,90,90)
    rho=rng.normal(size=shape);grid=gemmi.FloatGrid(*shape);grid.unit_cell=cell;grid.spacegroup=gemmi.SpaceGroup('P 1');grid.array[:]=rho
    plans,diagnostics=prepare_ncs(grid,np.ones(shape),metric_filter(cell,shape,4.),[{'rotation':np.eye(3).tolist(),'translation':[17.,13.,7.]}])
    assert max(v['active_fraction'] for v in diagnostics)<.001
    np.testing.assert_array_equal(apply_ncs(rho,cell,[]),rho)


def test_skew_cell_fractional_translation_resamples_correctly():
    cell=gemmi.UnitCell(30,35,40,90,100,110);shape=(30,35,40)
    values=np.indices(shape)[0].astype(float)
    t=np.array(cell.orth.mat)@np.array([1/30,0,0])
    result=resample(values,cell,np.eye(3),t)
    np.testing.assert_allclose(result,np.roll(values,-1,axis=0),atol=1e-10)


def test_bad_or_identity_ncs_operator_is_not_silently_accepted():
    with pytest.raises(ValueError,match='Identity'):validate_operators([{'rotation':np.eye(3).tolist(),'translation':[0,0,0]}])
    with pytest.raises(ValueError,match='proper'):validate_operators([{'rotation':np.diag([1,1,-1]).tolist(),'translation':[1,0,0]}])
