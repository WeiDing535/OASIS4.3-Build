"""Independent fixed-Q uniform-phase/noise calibration tests."""
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.integrate import quad

from phasebuild.pplus_scale import (
    fit_phase_marginal_scale, score_phase_marginal,
    PhaseMarginalQuadratureError,
)


def generated(alpha, seed=807231, n=1800):
    rng = np.random.default_rng(seed)
    q = rng.uniform(.5, 1.8, n)
    sigma = rng.uniform(.18, .45, n)
    d = 2 * alpha * q * np.cos(rng.uniform(0, 2*np.pi, n)) + rng.normal(0, sigma)
    return SimpleNamespace(dano=d, sigdano=sigma, centric=np.zeros(n, bool)), q


@pytest.mark.parametrize('true_alpha', [.08, .35, 1.1])
def test_uniform_phase_noise_scale_recovery(true_alpha):
    data, q = generated(true_alpha)
    alpha, meta = fit_phase_marginal_scale(data, q)
    # Finite-sample recovery, rather than an implementation-mirroring identity.
    assert alpha == pytest.approx(true_alpha, rel=.18, abs=.018)
    assert meta['train_loglikelihood'] >= meta['null_loglikelihood']
    assert meta['quadrature_levels'][-1]['converged']
    assert meta['reported_sigmaD_range'] == [float(data.sigdano.min()), float(data.sigdano.max())]


def test_null_and_unidentifiable_zero_substructure_are_explicit():
    data, q = generated(.5, n=60)
    data.dano[:] = 0
    alpha, meta = fit_phase_marginal_scale(data, q)
    assert alpha == 0 and meta['null_selected'] and meta['identifiable']
    flat, degenerate = fit_phase_marginal_scale(data, q * 0)
    assert flat == 0 and not degenerate['identifiable']
    assert degenerate['loglikelihood_gain_vs_null'] == 0


def test_zero_Q_rows_and_missing_centric_placeholders_are_handled_explicitly():
    data, q = generated(.4, n=90)
    q[:3] = 0
    data.centric[-1] = True
    data.dano[-1] = data.sigdano[-1] = q[-1] = np.nan
    alpha, meta = fit_phase_marginal_scale(data, q)
    assert np.isfinite(alpha) and meta['selected_zero_Q_rows'] == 3
    assert meta['selected_rows'] == 89
    values, scored = score_phase_marginal(data, q, alpha,
                                          samples=meta['samples_used'])
    assert np.isnan(values[-1]) and scored['selected_rows'] == 89
    expected = .5*(data.dano[:3]/data.sigdano[:3])**2 + np.log(data.sigdano[:3]*np.sqrt(2*np.pi))
    np.testing.assert_allclose(values[:3], expected, atol=1e-12)


def test_strong_signal_uses_only_training_quadrature_doubling():
    rng = np.random.default_rng(909871)
    q = rng.uniform(.6, 1.6, 180)
    sigma = np.full(len(q), .015)
    d = 2*.55*q*np.cos(rng.uniform(0, 2*np.pi, len(q))) + rng.normal(0, sigma)
    data = SimpleNamespace(dano=d, sigdano=sigma, centric=np.zeros(len(q), bool))
    alpha, meta = fit_phase_marginal_scale(data, q, samples=32)
    assert alpha == pytest.approx(.55, rel=.12)
    levels = meta['quadrature_levels']
    assert len(levels) > 1 and levels[-1]['converged']
    assert all(b['samples'] == 2*a['samples'] for a, b in zip(levels, levels[1:]))


def test_amplitude_units_and_substructure_phase_do_not_change_fit():
    data, q = generated(.45, n=350)
    a, _ = fit_phase_marginal_scale(data, q)
    units = 37.2
    scaled = SimpleNamespace(dano=data.dano*units, sigdano=data.sigdano*units,
                             centric=data.centric.copy())
    b, _ = fit_phase_marginal_scale(scaled, q*units)
    c, _ = fit_phase_marginal_scale(data, q*np.exp(1j*np.arange(len(q))))
    assert b == pytest.approx(a, rel=2e-6, abs=1e-8)
    assert c == pytest.approx(a, rel=2e-6, abs=1e-8)
    nll, _ = score_phase_marginal(data, q, a)
    other, _ = score_phase_marginal(scaled, q*units, a)
    np.testing.assert_allclose(other-nll, np.log(units), atol=1e-11)


def test_training_mask_and_centric_placeholders_cannot_affect_fit():
    data, q = generated(.55, n=200)
    data.centric[-20:] = True
    fit = np.zeros(200, bool); fit[:140] = True; fit[-20:] = True
    before = data.sigdano.copy()
    a, first = fit_phase_marginal_scale(data, q, fit_mask=fit)
    data.sigdano[-20:] = 1e80
    data.dano[140:180] = 1e50
    b, second = fit_phase_marginal_scale(data, q, fit_mask=fit)
    assert a == b and first == second
    np.testing.assert_array_equal(data.sigdano[:180], before[:180])
    assert first['selected_rows'] == 140
    assert first['requested_centric_rows_excluded'] == 20


def test_predictive_row_mask_is_explicit_and_alpha_is_not_refitted():
    data, q = generated(.4, n=30)
    data.centric[-1] = True
    mask = np.arange(30) % 2 == 0
    values, meta = score_phase_marginal(data, q, .4, mask)
    use = mask & ~data.centric
    assert np.isnan(values[~use]).all()
    assert meta['mean_nll'] == pytest.approx(values[use].mean())
    assert meta['alpha'] == .4 and not meta['alpha_refitted']


def test_gaussian_phase_integral_matches_independent_adaptive_quadrature():
    data = SimpleNamespace(dano=np.array([-.7, .2, 1.1]),
        sigdano=np.array([.23, .4, .17]), centric=np.zeros(3, bool))
    q = np.array([.6, .9, 1.2]); alpha = .8
    nll, _ = score_phase_marginal(data, q, alpha, samples=512)
    expected = []
    for d, sigma, magnitude in zip(data.dano, data.sigdano, q):
        density = quad(lambda theta: np.exp(-.5*((d-2*alpha*magnitude*np.cos(theta))/sigma)**2)
            / (sigma*np.sqrt(2*np.pi)*np.pi), 0, np.pi,
            epsabs=1e-12, epsrel=1e-12)[0]
        expected.append(-np.log(density))
    np.testing.assert_allclose(nll, expected, atol=2e-12, rtol=0)


def test_underresolved_prediction_raises_without_adapting_to_validation():
    data = SimpleNamespace(dano=np.array([.123]), sigdano=np.array([.0002]),
                           centric=np.zeros(1, bool))
    with pytest.raises(PhaseMarginalQuadratureError) as caught:
        score_phase_marginal(data, np.array([1.]), 1., samples=32)
    assert not caught.value.diagnostics['converged']


@pytest.mark.parametrize('change', ['sigma_zero', 'd_nan', 'q_nan', 'bad_mask', 'empty_mask', 'bad_shape'])
def test_invalid_selected_data_and_masks_fail(change):
    data, q = generated(.3, n=20); mask = None
    if change == 'sigma_zero': data.sigdano[0] = 0
    elif change == 'd_nan': data.dano[0] = np.nan
    elif change == 'q_nan': q[0] = np.nan
    elif change == 'bad_mask': mask = np.ones(20, int)
    elif change == 'empty_mask': mask = np.zeros(20, bool)
    elif change == 'bad_shape': q = q[:-1]
    with pytest.raises(ValueError):
        fit_phase_marginal_scale(data, q, fit_mask=mask)


@pytest.mark.parametrize('alpha', [-1., np.nan, np.inf, True, 1j, '.3'])
def test_invalid_fixed_alpha_is_rejected(alpha):
    data, q = generated(.2, n=10)
    with pytest.raises(ValueError): score_phase_marginal(data, q, alpha)
