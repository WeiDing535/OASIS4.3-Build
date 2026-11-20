import sys
from pathlib import Path
import numpy as np
import pytest
from oasis_py.data import Atoms,structure_factors
from test_science import dataset

from support.scoring import metrics


def test_allowed_origin_translation_does_not_become_a_false_failure():
    data=dataset('P 21 21 21',extent=4)
    atoms=Atoms([[.13,.27,.39],[.41,.12,.23]],['S','C'],[1,.7],[10,15])
    truth=structure_factors(data,atoms)
    shifted=Atoms(atoms.fractional+[0,.5,0],atoms.elements,atoms.occupancy,atoms.b_iso)
    calculated=structure_factors(data,shifted)
    np.testing.assert_allclose(abs(calculated),abs(truth),atol=1e-12)
    phase=np.angle(calculated);fom=np.ones(len(phase))
    assert metrics(data,phase,fom,truth)['phase_mae_deg']>60
    aligned=metrics(data,phase,fom,truth,origin_shift=[0,.5,0])
    assert aligned['phase_mae_deg']<1e-10
    assert aligned['fourier_coefficient_cc']==pytest.approx(1.)


def test_scoring_rejects_translation_that_changes_symmetry_setting():
    data=dataset('P 21 21 21')
    with pytest.raises(ValueError,match='space-group setting'):
        metrics(data,np.zeros(len(data.f)),np.ones(len(data.f)),
                np.ones(len(data.f),complex),origin_shift=[0,.25,0])


def test_polar_axis_continuous_origin_is_allowed():
    data=dataset('P 21');shift=[0,.12345,0]
    phase=2*np.pi*(data.hkl@shift)
    score=metrics(data,phase,np.ones(len(phase)),np.ones(len(phase),complex),origin_shift=shift)
    assert score['phase_mae_deg']==0
