"""Opt-in original-likelihood-compatible, KL-bounded P+ branch updates.

This is a conditional-shape two-state research approximation. It is not a
replacement observation likelihood or an exact continuous SAD posterior.
"""
import numpy as np
from scipy.special import expit, logsumexp, xlogy
from oasis_py.clipper_likelihood import reflection_classes
from .pplus_conditional import (ConditionalPplusResult, prepare_conditional_branches,
                                iterate_conditional_pplus)


def _mass_array(log_mass):
    lm = np.asarray(log_mass, dtype=float)
    if (lm.ndim != 2 or lm.shape[1] != 2 or not len(lm) or np.isnan(lm).any()
            or np.isposinf(lm).any() or np.any(np.all(np.isneginf(lm), axis=1))
            or np.any(abs(logsumexp(lm, axis=1)) > 1e-10)):
        raise ValueError('Normalized, nonempty two-branch log masses are required')
    return lm


def _kl(p, lm):
    with np.errstate(invalid='ignore'):
        a = np.where(p == 0, 0, p * lm[:, 0])
        b = np.where(p == 1, 0, (1-p) * lm[:, 1])
    return np.maximum(0, xlogy(p, p) + xlogy(1-p, 1-p) - a - b)


def bernoulli_kl(probability, log_mass):
    """KL(q || original branch mass), with exact zero-mass conventions."""
    lm = _mass_array(log_mass)
    p = np.asarray(probability, dtype=float)
    if p.shape != (len(lm),) or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('Branch probabilities must be finite and in [0,1]')
    return _kl(p, lm)


def bounded_branch_mass(log_mass, score, budget):
    """Move toward the variational update, bounded relative to ORIGINAL mass.

    Bisection on the probability segment uses convexity of Bernoulli KL.
    budget=0 is an explicit mathematical no-update control; None is unbounded.
    """
    lm = _mass_array(log_mass)
    s = np.asarray(score, dtype=float)
    if s.shape != (len(lm),) or not np.isfinite(s).all():
        raise ValueError('Invalid P+ branch field')
    if budget is not None and (isinstance(budget, bool) or not np.isfinite(budget) or budget < 0):
        raise ValueError('Invalid probability KL budget')
    base = np.clip(np.exp(lm[:, 0]), 0, 1)
    if budget == 0:
        return base.copy()
    proposed = expit(lm[:, 0] - lm[:, 1] + s)
    deterministic = np.isneginf(lm).any(axis=1)
    proposed[deterministic] = base[deterministic]
    if budget is None:
        return proposed
    bad = _kl(proposed, lm) > budget
    lo, hi = np.zeros(len(lm)), np.ones(len(lm))
    for _ in range(60):
        mid = (lo + hi) / 2
        p = base + mid * (proposed - base)
        over = _kl(p, lm) > budget
        hi = np.where(over, mid, hi)
        lo = np.where(over, lo, mid)
    return np.where(bad, base + lo * (proposed - base), proposed)


def iterate_trusted_pplus(data, branches, relations, *, kl_budget=.0005,
                          cycles=100, damping=.25, tolerance=1e-6):
    """Update P+ using triplet fields with frozen conditional branch shapes.

    The supplied original observation likelihood is kept fixed, never
    accumulated. No additional Sim field is added: the demonstrated SEC17
    frontend already models substructure scattering in its SAD likelihood.
    A positive budget is mandatory; callers must audit actual update size.
    """
    if (isinstance(kl_budget, bool) or kl_budget is None or not np.isfinite(kl_budget)
            or kl_budget <= 0 or type(cycles) is not int or cycles < 2
            or isinstance(damping, bool) or not np.isfinite(damping) or not 0 < damping <= 1
            or not np.isfinite(tolerance) or tolerance <= 0):
        raise ValueError('Invalid trusted P+ controls; positive KL budget required')
    # Reuse existing crystal identity, centric-state and shape validation.
    iterate_conditional_pplus(data, branches, np.zeros(len(data.f), complex),
                              relations, cycles=1)
    p0 = np.exp(branches.log_mass[:, 0])
    p = p0.copy()
    history = []
    for cycle in range(cycles):
        previous = p.copy()
        first = p * branches.first[:, 0] + (1-p) * branches.first[:, 1]
        field = np.asarray(relations.field(np.angle(first), abs(first)), dtype=complex)
        if field.shape != (len(p),) or not np.isfinite(field).all():
            raise ValueError('Invalid relation field')
        score = np.real(field * np.conj(branches.first[:, 0] - branches.first[:, 1]))
        proposed = bounded_branch_mass(branches.log_mass, score, kl_budget)
        p = (1-damping) * p + damping * proposed
        delta = float(abs(p-previous).max())
        kl = _kl(p, branches.log_mass)
        history.append({'iteration': cycle+1, 'maximum_probability_change': delta,
                        'mean_KL_from_original': float(kl.mean()),
                        'max_KL_from_original': float(kl.max())})
        if cycle >= 2 and delta < tolerance:
            break
    first = p * branches.first[:, 0] + (1-p) * branches.first[:, 1]
    second = p * branches.second[:, 0] + (1-p) * branches.second[:, 1]
    phase, fom = np.angle(first) % (2*np.pi), np.minimum(abs(first), 1)
    _, centric, allowed = reflection_classes(data)
    signed = np.real(first[centric] / np.exp(1j*allowed[centric]))
    phase[centric] = (allowed[centric] + np.where(signed < 0, np.pi, 0)) % (2*np.pi)
    metadata = {'method': 'original-observation-compatible-KL-bounded-two-state-Pplus',
                'approximation': 'Frozen conditional shapes, variational branch masses; not exact continuous posterior',
                'budget_nats_per_reflection': float(kl_budget), 'damping': float(damping),
                'iterations': len(history), 'converged': delta < tolerance,
                'observation_likelihood_accumulated': False,
                'original_observation_HL_modified': False, 'posterior_HL_exported': False,
                'additional_Sim_field_added': False, 'role': 'Map-seed moments only',
                'mean_abs_probability_change': float(abs(p-p0).mean()),
                'rows_probability_change_over_0p01': int((abs(p-p0) > .01).sum())}
    return ConditionalPplusResult(phase, fom, p, first, second, history, metadata)


def phase_trusted_pplus(data, observation_hl, reference, relations, *, samples=1440, **kwargs):
    """Construct both branch shapes and their P+ update from the SAME HL."""
    branches = prepare_conditional_branches(data, observation_hl, reference, samples=samples)
    return iterate_trusted_pplus(data, branches, relations, **kwargs)
