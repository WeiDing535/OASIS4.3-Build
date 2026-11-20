"""Seed-only P+ ablation using an observed anomalous-amplitude scale.

Only the doublet geometry changes. The existing calibration's attenuation,
partial-structure field, triplet graph and experimental HL remain fixed.
This is not a jointly recalibrated SAD likelihood or a measurement of f''.
"""
import hashlib
import numpy as np
from oasis_py.data import structure_factors
from oasis_py.phasing import PhaseResult, iterate_pplus
from oasis_py.sad import fit_anomalous_scale
from .pplus_experimental import PplusCalibration, _calibration_identity


def phase_pplus_observed_scale(data, composition, heavy, fpp, *,
                               calibration_context, centric_singletons=None,
                               cycles=5, aoe=.5, kmin=.03,
                               max_neighbors=600, max_per_target=96):
    """Change doublet angles using a work-observation-only scale estimate.

    The caller must supply the same work-only reflections, heavy substructure,
    composition and settings used to prepare ``calibration_context``. This
    function does not select/exclude a free set: that is an input contract.
    Its returned HL is the frozen legacy observation model, whereas phase,
    FOM and PPLUS describe the separate map seed. They must not be treated
    as moments of the returned HL or multiplied into that likelihood.
    """
    context = calibration_context
    if not isinstance(context, PplusCalibration):
        raise ValueError('A frozen PplusCalibration context is required')
    if type(cycles) is not int or cycles < 1:
        raise ValueError('P+ cycles must be a positive integer')
    if data.dano is None or data.sigdano is None:
        raise ValueError('P+ observed-scale seed requires SAD observations')
    n = len(data.f)
    single = (np.zeros(n, dtype=bool) if centric_singletons is None
              else np.asarray(centric_singletons, dtype=bool))
    if single.shape != (n,) or np.any(single & ~data.centric) or np.all(single):
        raise ValueError('Invalid centric singleton mask')
    if (not fpp or any(el not in fpp or not np.isfinite(fpp[el]) or fpp[el] <= 0
                       for el in heavy.elements)):
        raise ValueError('Positive finite f-double-prime required')
    identity = _calibration_identity(data, composition, heavy, fpp, single,
                                     aoe, kmin, max_neighbors, max_per_target)
    if context.identity != identity:
        raise ValueError('Frozen P+ calibration does not match input/settings')
    if (context.error.shape != (n,) or context.sim.shape != (n,)
            or context.observations.shape != (n, 4)
            or context.graph.nreflections != n):
        raise ValueError('Frozen P+ calibration has inconsistent arrays')

    raw = structure_factors(data, heavy, fpp)
    alpha, fit = fit_anomalous_scale(data, raw)
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError('No positive observed anomalous scale; seed ablation unavailable')
    sf = alpha * raw
    reference = np.angle(sf)
    ratio = data.dano / (2 * np.maximum(np.abs(sf), 1e-8))
    difference = np.arccos(np.clip(ratio, -1, 1))
    difference[data.centric] = np.pi / 2
    phase, fom, probability, history = iterate_pplus(
        reference, difference, context.error, context.sim, context.graph,
        cycles=cycles)
    observations = context.observations.copy()
    prior_digest = hashlib.sha256(data.hkl.tobytes() + observations.tobytes()).hexdigest()
    metadata = {
        **context.metadata,
        'method': 'pplus_observed_scale_geometry_seed_experimental',
        'role': 'Map seed only; fixed legacy observation HL retained separately',
        'calibration_context_sha256': context.identity,
        'anomalous_scale_fit': fit,
        'geometry_anomalous_scale': float(alpha),
        'scale_interpretation': 'Effective anomalous amplitude scale; not measured f-double-prime',
        'seed_only_geometry_ablation': True,
        'fully_recalibrated_likelihood': False,
        'cosine_rank': False,
        'error_weights_reused': True,
        'partial_structure_field_reused': True,
        'triplet_graph_reused': True,
        'observation_HL_unchanged': True,
        'fixed_observation_HL_sha256': prior_digest,
        'hl_definition': 'Exact copy of frozen legacy calibrated SAD observation model; not map-seed moments',
        'centric_singletons_used_in_scale_fit': False,
        'scale_fit_reflection_selection': 'Working acentric observations only; centric rows excluded',
        'geometry_clipped_doublet_fraction_acentric': float(np.mean(np.abs(ratio[~data.centric]) > 1)),
        'cycles': history,
    }
    return PhaseResult(phase, fom, probability, observations, metadata)
