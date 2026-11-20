"""Separate SAD observation likelihood from phase-improvement posterior fields.

The upstream numerical result remains unchanged for density modification.
Only this observation likelihood is exported as the fixed REFMAC restraint.
"""
import numpy as np
from oasis_py.data import structure_factors
from oasis_py.phasing import PhaseResult, sad_hl
from oasis_py.sad import posterior_moments


def observation_prior(data, result, heavy, fpp, method):
    if method == 'likelihood':
        if result.metadata.get('method') != 'sad_likelihood':
            raise ValueError('Observation prior requires matching SAD likelihood metadata')
        alpha = result.metadata['anomalous_scale']
        error = result.metadata['model_error_fraction']
        if not np.isfinite(alpha) or alpha < 0 or not np.isfinite(error) or not 0 <= error <= 1:
            raise ValueError('Invalid saved SAD observation scale/error')
        # Reuse the fitted upstream scale; never refit against a protein model.
        sf = alpha * structure_factors(data, heavy, fpp)
        sigma = np.hypot(data.sigdano, error * np.abs(sf))
        hl = sad_hl(np.angle(sf), np.abs(sf), data.dano, sigma)
        details = dict(anomalous_scale=alpha, model_error_fraction=error,
                       definition='Upstream SAD observation likelihood without final triplet/partial-structure field')
    elif method == 'pplus':
        # P+ already returns its error-calibrated observation likelihood as HL.
        hl = np.array(result.hl, dtype=float, copy=True)
        details = dict(error_scale=result.metadata['error_scale'],
                       definition='Upstream P+ error-calibrated SAD observation likelihood')
    else:
        raise ValueError('Unsupported SAD observation-prior method')
    phase, fom = posterior_moments(data, hl)
    return PhaseResult(phase, fom, None, hl, dict(
        role='Fixed SAD observation prior for refinement; not the phase-improvement posterior',
        phase_method=method, method='sad_observation_prior',
        structural_triplet_field_included=False, protein_model_used=False, **details))
