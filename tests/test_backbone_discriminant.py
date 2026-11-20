import numpy as np
import pytest
from phasebuild.backbone_discriminant import (BackboneDiscriminant, fit_discriminant,
    discrimination_summary, supported_segments)


def test_logistic_model_learns_decoy_separation_and_roundtrips(tmp_path):
    rng = np.random.RandomState(51)
    offsets = np.zeros((2, 3))
    values = rng.normal(size=(120, 14))
    labels = (values[:, 0] + values[:, 1] > 0).astype(int)
    model = fit_discriminant(values, labels, offsets, {'purpose': 'synthetic test'}, regularization=2.)
    scores = model.score_features(values)
    assert discrimination_summary(labels, scores)['synthetic_auc'] > .97
    model.save(tmp_path / 'model.npz')
    restored = BackboneDiscriminant.load(tmp_path / 'model.npz')
    np.testing.assert_allclose(restored.score_features(values), scores)


def test_pruning_joins_only_bounded_internal_gaps():
    assert supported_segments([0., 2., 3., 0., 2., 3., 0.], 1., 4, 1) == [(1, 6)]
    assert supported_segments([2., 3., 0., 0., 2., 3.], 1., 4, 1) == []
    assert supported_segments([2., 3., 0., 0., 2., 3.], 1., 2, 1) == [(0, 2), (4, 6)]


def test_discrimination_auc_ties_and_fixed_thresholds():
    result = discrimination_summary([1, 0, 1, 0], [2., 0., 1., 1.], {'fixed': 1.5})
    assert result['synthetic_auc'] == pytest.approx(.875)
    assert result['thresholds']['fixed']['positive_retention'] == .5
    assert result['thresholds']['fixed']['negative_retention'] == 0.


def test_model_rejects_nonfinite_parameters():
    with pytest.raises(ValueError, match='Invalid discriminant'):
        BackboneDiscriminant(np.zeros((1, 3)), np.zeros(11), np.ones(11), np.zeros(11), np.nan, {})
