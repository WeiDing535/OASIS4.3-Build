import numpy as np
from phasebuild.residue_logistic import (ResidueLogisticProfile, fit_residue_logistic,
                                         perturb_backbones)


def test_balanced_logistic_learns_residue_signals_and_roundtrips(tmp_path):
    random = np.random.RandomState(4)
    labels = np.repeat(np.arange(20), 6)
    values = np.eye(20)[labels]+random.normal(scale=.03, size=(len(labels), 20))
    profile = fit_residue_logistic(values, labels, np.zeros((20, 3)))
    assert np.mean(profile.scores_from_features(values).argmax(axis=1) == labels) == 1.
    assert profile.metadata['optimizer_success']
    path = tmp_path/'target.npz'; profile.save(path)
    loaded = ResidueLogisticProfile.load(path)
    np.testing.assert_array_equal(loaded.scores_from_features(values), profile.scores_from_features(values))


def test_score_magnitude_tracks_evidence_and_is_not_forced_to_unit_variance():
    weights = np.arange(20).reshape(1, 20).astype(float)
    profile = ResidueLogisticProfile(np.zeros((1, 3)), np.zeros(1), np.ones(1), weights, np.zeros(20), {})
    scores = profile.scores_from_features(np.array([[0.], [.2], [1.]]))
    np.testing.assert_array_equal(scores[0], np.zeros(20))
    np.testing.assert_allclose(scores[1], .2*scores[2])
    assert scores[1].std() < scores[2].std()


def test_local_augmentation_preserves_intraresidue_geometry_and_is_reproducible():
    xyz = np.random.RandomState(5).normal(size=(12, 4, 3))
    actual = perturb_backbones(xyz, np.random.RandomState(3))
    other = perturb_backbones(xyz, np.random.RandomState(3))
    np.testing.assert_array_equal(actual, other)
    original_dist = np.linalg.norm(xyz[:, :, None]-xyz[:, None, :], axis=-1)
    actual_dist = np.linalg.norm(actual[:, :, None]-actual[:, None, :], axis=-1)
    np.testing.assert_allclose(actual_dist, original_dist, atol=1e-12)
    assert not np.allclose(actual, xyz)
