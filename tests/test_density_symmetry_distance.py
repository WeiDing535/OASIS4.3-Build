"""Regression: a map-derived UnitCell must retain crystallographic images."""
import gemmi
import numpy as np
import pytest
from phasebuild.data import Density


@pytest.mark.parametrize('spacegroup',['P 1 21 1','C 1 2 1','P 21 21 21'])
def test_density_distance_includes_screw_rotation_and_centering_images(tmp_path,spacegroup):
    mtz=gemmi.Mtz(with_base=True)
    mtz.spacegroup=gemmi.find_spacegroup_by_name(spacegroup)
    mtz.set_cell_for_all(gemmi.UnitCell(40.,44.,48.,90.,90.,90.))
    mtz.add_dataset('synthetic')
    mtz.add_column('FWT','F');mtz.add_column('PHWT','P')
    hkl=gemmi.make_miller_array(mtz.cell,mtz.spacegroup,4.)
    random=np.random.RandomState(1729)
    mtz.set_data(np.column_stack([hkl,random.uniform(1.,20.,len(hkl)),
                                 random.uniform(0.,360.,len(hkl))]).astype(np.float32))
    path=tmp_path/'non_p1_map.mtz';mtz.write_to_file(str(path))
    density=Density(path)
    assert len(density.cell.images)==len(density.sg.operations())-1
    original=gemmi.Fractional(.11,.19,.23)
    for operation in list(density.sg.operations())[1:]:
        image=gemmi.Fractional(*operation.apply_to_xyz(list(original)))
        a=density.cell.orthogonalize(original);b=density.cell.orthogonalize(image)
        b.x+=.1
        assert density.distance(np.array(list(a)),np.array(list(b)))==pytest.approx(.1,abs=1e-8)
        # Explicitly show that this pair cannot be explained by translations.
        p1_cell=gemmi.UnitCell(*density.cell.parameters)
        assert p1_cell.find_nearest_image(a,b).dist()>1.
