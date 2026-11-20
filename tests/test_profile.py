import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from phasebuild.profile import DensityProfile, patches


class GaussianDensity:
    def __init__(self, atoms):
        self.atoms = atoms

    def sample(self, xyz):
        delta = xyz[..., None, :] - self.atoms
        return np.exp(-np.sum(delta**2, axis=-1) / .5).sum(axis=-1)


def test_local_profile_is_rotation_translation_invariant_and_direction_sensitive(tmp_path):
    xyz = np.array([[-1.3,.6,0.], [0.,0.,0.], [1.5,0.,0.], [2.,1.,0.]])
    offsets = np.array(np.meshgrid(np.arange(-2., 2.1, 1.), np.arange(-2., 2.1, 1.),
                                   np.arange(-2., 2.1, 1.), indexing='ij')).reshape(3,-1).T
    original = GaussianDensity(xyz)
    target = patches(original, xyz, offsets)
    profile = DensityProfile(offsets, target, np.ones(len(offsets)), {'purpose':'test'})
    rotation = Rotation.from_euler('xyz', [32,67,15], degrees=True).as_matrix()
    shifted = xyz @ rotation.T + [10., -12., 4.]
    assert profile.score(original, xyz) == pytest.approx(profile.score(GaussianDensity(shifted), shifted), abs=1e-6)
    wrong = xyz.copy()
    wrong[[0,2]] = wrong[[2,0]]
    assert profile.score(original, xyz) > profile.score(original, wrong) + .5
    path = tmp_path / 'target.npz'
    profile.save(path)
    restored = DensityProfile.load(path)
    assert restored.score(original, xyz) == pytest.approx(profile.score(original, xyz))


def test_invalid_profile_rejected():
    with pytest.raises(ValueError, match='Invalid density-profile'):
        DensityProfile(np.zeros((3,3)), np.zeros(3), -np.ones(3), {})
