"""Fixed-substructure scale from a noise-integrated small-signal SAD model.

DANO = 2 * alpha * |Q| * cos(theta) + Gaussian(reported_sigmaD), with
uniform theta. This is not a joint Bijvoet likelihood, a measurement of
f-double-prime, or a criterion for ranking sites or hands. No phases, protein
models or free-set information are accepted or inspected.
"""
from __future__ import annotations

import hashlib
from numbers import Real
import numpy as np
from scipy.optimize import minimize_scalar


_MAX_FIT_SAMPLES = 2048
_MEAN_QUADRATURE_TOL = 2e-5
_MAX_QUADRATURE_TOL = 2e-3
_LOG_SQRT_2PI = .5 * np.log(2 * np.pi)


class PhaseMarginalQuadratureError(ValueError):
    """The predeclared quadrature sequence did not resolve the likelihood."""

    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics


def _samples(value):
    if type(value) is not int or not 32 <= value <= _MAX_FIT_SAMPLES:
        raise ValueError('samples must be an integer from 32 to 2048')
    return value


def _selected(data, imaginary_sf, mask):
    if data.dano is None or data.sigdano is None:
        raise ValueError('Measured DANO and reported sigmaD are required')
    d = np.asarray(data.dano, dtype=float)
    sigma = np.asarray(data.sigdano, dtype=float)
    sf = np.asarray(imaginary_sf)
    centric = np.asarray(data.centric)
    if (d.ndim != 1 or sigma.shape != d.shape or sf.shape != d.shape
            or centric.shape != d.shape or centric.dtype != np.bool_):
        raise ValueError('DANO, sigmaD, substructure and centric shapes differ')
    requested = np.ones(len(d), dtype=bool) if mask is None else np.asarray(mask)
    if requested.shape != d.shape or requested.dtype != np.bool_:
        raise ValueError('mask must be a boolean vector matching the data')
    use = requested & ~centric
    if not use.any():
        raise ValueError('At least one selected acentric observed pair is required')
    d, sigma, sf = d[use], sigma[use], sf[use]
    if (not np.isfinite(d).all() or not np.isfinite(sigma).all()
            or not np.isfinite(sf).all() or np.any(sigma <= 0)):
        raise ValueError('Selected acentric DANO/Q must be finite and sigmaD positive')
    q = np.abs(sf).astype(float)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        normalized = np.r_[d / sigma, q / sigma]
    if not np.isfinite(q).all() or not np.isfinite(normalized).all():
        raise ValueError('Selected observations exceed the supported numerical range')
    meta = {
        'total_rows': len(use), 'requested_rows': int(requested.sum()),
        'selected_rows': int(use.sum()),
        'requested_centric_rows_excluded': int(np.sum(requested & centric)),
        'selected_zero_Q_rows': int(np.sum(q == 0)),
        'selection_sha256': hashlib.sha256(use.tobytes()).hexdigest(),
        'reported_sigmaD_range': [float(sigma.min()), float(sigma.max())],
        'DANO_range': [float(d.min()), float(d.max())],
        'Q_magnitude_range': [float(q.min()), float(q.max())],
        'sigma_policy': 'Unmodified reported sigmaD; no AOE shrinkage or fitted noise inflation',
    }
    return use, d, sigma, q, meta


def _nll(d, sigma, q, alpha, samples):
    """Uniform-theta Gaussian mixture, symmetry folded to theta in (0, pi)."""
    if alpha == 0 or not np.any(q):
        return .5 * (d / sigma) ** 2 + np.log(sigma) + _LOG_SQRT_2PI
    cosine = np.cos(np.pi * (np.arange(samples) + .5) / samples)
    result = np.empty(len(d))
    block = max(1, min(256, 1048576 // samples))
    for start in range(0, len(d), block):
        stop = start + block
        with np.errstate(over='ignore', invalid='ignore'):
            z = ((d[start:stop, None]
                  - (2 * alpha * q[start:stop, None]) * cosine)
                 / sigma[start:stop, None])
            exponent = -.5 * z * z
        largest = np.max(exponent, axis=1)
        log_mean = np.full(len(largest), -np.inf)
        finite = np.isfinite(largest)
        log_mean[finite] = (largest[finite]
            + np.log(np.exp(exponent[finite] - largest[finite, None]).sum(axis=1))
            - np.log(samples))
        result[start:stop] = (np.log(sigma[start:stop]) + _LOG_SQRT_2PI
                              - log_mean)
    return result


def _quadrature_check(d, sigma, q, alpha, samples):
    coarse = _nll(d, sigma, q, alpha, samples)
    fine = _nll(d, sigma, q, alpha, 2 * samples)
    error = np.abs(coarse - fine)
    maximum, mean = float(error.max()), float(error.mean())
    return coarse, {
        'fit_samples': samples, 'check_samples': 2 * samples,
        'max_abs_row_loglikelihood_difference': maximum,
        'mean_abs_row_loglikelihood_difference': mean,
        'max_tolerance': _MAX_QUADRATURE_TOL,
        'mean_tolerance': _MEAN_QUADRATURE_TOL,
        'converged': bool(np.isfinite(error).all()
            and maximum <= _MAX_QUADRATURE_TOL
            and mean <= _MEAN_QUADRATURE_TOL),
    }


def _rms(x):
    maximum = float(np.max(np.abs(x)))
    return maximum * float(np.sqrt(np.mean((x / maximum) ** 2))) if maximum else 0.


def _fit_level(d, sigma, q, samples, previous=None):
    # Work in dimensionless beta, so changing all amplitude units cannot
    # alter the search grid or its absolute stopping tolerance.
    qref, dref = _rms(q), max(_rms(d), _rms(sigma))
    factor = dref / (2 * qref)
    positive = q > 0
    ratio = ((np.abs(d[positive] / dref) + 6 * sigma[positive] / dref)
             / (q[positive] / qref))
    upper = max(32., 4 * float(np.quantile(ratio, .95)))
    if not np.isfinite(upper) or not np.isfinite(factor) or factor <= 0:
        raise ValueError('Scale search range is not numerically representable')
    evaluations = 0

    def objective(beta):
        nonlocal evaluations
        evaluations += 1
        return float(np.mean(_nll(d, sigma, q, factor * beta, samples)))

    null = objective(0.)
    if not np.isfinite(null):
        raise ValueError('Null likelihood exceeds the supported numerical range')
    expansions = 0
    while True:
        grid = np.unique(np.r_[np.geomspace(1e-7, upper, 97),
            np.linspace(0, min(4., upper), 33)[1:],
            [] if previous is None else [previous / factor]])
        values = np.array([objective(x) for x in grid])
        candidates = [(null, 0.)] + list(zip(values, grid))
        minima = [i for i in range(1, len(grid) - 1)
                  if values[i] <= values[i-1] and values[i] <= values[i+1]]
        # Refine several basins, rather than assuming a globally unimodal
        # objective. Quadrature checks separately reject unresolved spikes.
        for i in sorted(minima, key=lambda j: values[j])[:8]:
            opt = minimize_scalar(lambda x: objective(float(np.exp(x))),
                bounds=(np.log(grid[i-1]), np.log(grid[i+1])), method='bounded',
                options={'xatol': 1e-9, 'maxiter': 80})
            if opt.success and np.isfinite(opt.fun):
                candidates.append((float(opt.fun), float(np.exp(opt.x))))
        value, beta = min(candidates)
        if value >= null - 1e-10:
            value, beta = null, 0.
        if beta < .98 * upper:
            break
        expansions += 1
        if expansions > 3:
            raise ValueError('Scale optimum remains at the expanded search boundary')
        upper *= 4
    return float(factor * beta), {
        'samples': samples, 'objective_evaluations': evaluations,
        'alpha_search_bounds': [0., float(factor * upper)],
        'dimensionless_beta_search_bounds': [0., float(upper)],
        'search_boundary_expansions': expansions,
        'search': 'Explicit null plus logarithmic/linear grids and up to eight refined basins',
    }


def fit_phase_marginal_scale(data, imaginary_sf, *, fit_mask=None, samples=256):
    """Return effective alpha and training-only diagnostics for fixed |Q|.

    Predeclared quadrature doubling may refit the same training objective,
    through at most 2048 nodes; the final check may use 4096. Validation data
    never choose alpha or numerical precision. Unresolved integration raises
    PhaseMarginalQuadratureError with all attempted-level diagnostics.
    """
    level = _samples(samples)
    _, d, sigma, q, meta = _selected(data, imaginary_sf, fit_mask)
    levels, previous = [], None
    if not np.any(q):
        alpha, search = 0., {'samples': level, 'objective_evaluations': 0,
                            'alpha_search_bounds': [0., 0.]}
    while True:
        if np.any(q):
            alpha, search = _fit_level(d, sigma, q, level, previous)
        nll, check = _quadrature_check(d, sigma, q, alpha, level)
        levels.append({**search, **check, 'alpha': alpha})
        if check['converged']:
            break
        if level >= _MAX_FIT_SAMPLES:
            raise PhaseMarginalQuadratureError(
                'Phase-marginal training quadrature did not converge through 2048 nodes', levels)
        previous, level = alpha, min(2 * level, _MAX_FIT_SAMPLES)
    null = _nll(d, sigma, q, 0., level)
    meta.update(method='fixed_HA_uniform_phase_Gaussian_DANO_marginal_v1',
        anomalous_scale=alpha, samples_requested=samples, samples_used=level,
        quadrature_levels=levels, train_loglikelihood=float(-nll.sum()),
        train_mean_nll=float(nll.mean()), null_loglikelihood=float(-null.sum()),
        null_mean_nll=float(null.mean()),
        loglikelihood_gain_vs_null=float((null - nll).sum()),
        identifiable=bool(np.any(q)), null_selected=bool(alpha == 0),
        quadrature='Equal-weight midpoint theta nodes in (0,pi); symmetry-folded uniform full circle',
        interpretation='Effective amplitude scale under a small-signal approximation; not f-double-prime or a site/hand ranking score')
    return alpha, meta


def score_phase_marginal(data, imaginary_sf, alpha, mask=None, *, samples=256):
    """Predict per-row/mean NLL at fixed alpha; excluded rows are NaN.

    This never refits alpha or increases samples from validation scores. A
    fixed next-grid accuracy check raises explicitly if the requested grid
    is insufficient. Include masks in any interpretation of returned means.
    """
    level = _samples(samples)
    if (isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, Real)
            or not np.isfinite(alpha) or alpha < 0):
        raise ValueError('alpha must be a finite nonnegative scalar')
    use, d, sigma, q, meta = _selected(data, imaginary_sf, mask)
    values, check = _quadrature_check(d, sigma, q, float(alpha), level)
    if not check['converged']:
        raise PhaseMarginalQuadratureError(
            'Requested predictive quadrature is insufficient; alpha was not refitted', check)
    result = np.full(len(use), np.nan)
    result[use] = values
    meta.update(alpha=float(alpha), samples_used=level, quadrature_check=check,
                mean_nll=float(values.mean()), loglikelihood=float(-values.sum()),
                output='Original row order; unselected and centric rows are NaN',
                alpha_refitted=False)
    return result, meta
