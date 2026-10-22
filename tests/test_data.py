import gemmi
import numpy as np
import pytest
from phasebuild.data import Density


def make_mtz(path, fom, zero=False):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz.set_cell_for_all(gemmi.UnitCell(12, 12, 12, 90, 90, 90))
    mtz.add_dataset('test')
    for label, kind in [('FWT','F'), ('PHWT','P'), ('FOM','W')]:
        mtz.add_column(label, kind)
    mtz.set_data(np.array([[1,0,0,0 if zero else 10,0,fom],
                           [0,1,0,0 if zero else 8,45,fom],
                           [0,0,1,0 if zero else 6,90,fom]], dtype=np.float32))
    mtz.write_to_file(str(path))


def test_map_coefficients_not_weighted_twice_and_periodic_sampling(tmp_path):
    a, b = tmp_path / 'a.mtz', tmp_path / 'b.mtz'
    make_mtz(a, .1)
    make_mtz(b, 1.)
    x, y = Density(a), Density(b)
    np.testing.assert_allclose(np.array(x.grid), np.array(y.grid))
    p = np.array([[2.,3.,4.], [14.,3.,4.]])
    assert x.sample(p)[0] == pytest.approx(x.sample(p)[1], abs=1e-6)


def test_zero_variance_map_rejected(tmp_path):
    path = tmp_path / 'empty.mtz'
    make_mtz(path, 1., zero=True)
    with pytest.raises(ValueError, match='variance'):
        Density(path)
