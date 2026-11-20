"""Opt-in observation-consistent P+ seeds under a small-signal SAD model.

Scale fitting and branch shapes use the same reported DANO uncertainties.
There is no AOE sigma contraction, clipped geometric doublet, reference
protein, or external phasing result. This remains an approximate SAD model,
not a joint Bijvoet likelihood or a validated default implementation.
"""
from dataclasses import dataclass
import hashlib
import json
import numpy as np
from oasis_py.data import structure_factors
from oasis_py.phasing import PhaseResult, sad_hl, _sim_field
from oasis_py.relations import build_relations
from oasis_py.sad import posterior_moments
from .pplus_experimental import prepare_pplus_normalization, _normalization_identity
from .pplus_scale import fit_phase_marginal_scale
from .pplus_conditional import prepare_conditional_branches, iterate_conditional_pplus


def _data_identity(data):
    digest=hashlib.sha256()
    for value in [data.hkl,data.f,data.sigf,data.dano,data.sigdano]:
        if value is None:raise ValueError('Measured SAD differences and uncertainties required')
        array=np.ascontiguousarray(value)
        digest.update(str(array.dtype).encode());digest.update(str(array.shape).encode())
        digest.update(array.tobytes())
    digest.update(json.dumps({'cell':list(data.cell.parameters),
        'symmetry':[op.triplet() for op in data.spacegroup.operations()]},sort_keys=True).encode())
    return digest.hexdigest()


@dataclass(frozen=True)
class ObservationPplusContext:
    identity: str
    data_identity: str
    reference: np.ndarray
    observations: np.ndarray
    sim: np.ndarray
    graph: object
    branches: object
    metadata: dict


def _context_identity(data_identity,reference,observations,sim,graph,branches,metadata):
    digest=hashlib.sha256(data_identity.encode())
    arrays=[reference,observations,sim,branches.log_mass,branches.first,branches.second]
    arrays += [getattr(graph,name) for name in
        ['target','left','right','left_sign','right_sign','shift','kappa']]
    for value in arrays:
        array=np.ascontiguousarray(value)
        digest.update(str(array.dtype).encode());digest.update(str(array.shape).encode())
        digest.update(array.tobytes())
    digest.update(json.dumps(metadata,sort_keys=True).encode())
    return digest.hexdigest()


def _check_context(data,context):
    if (not isinstance(context,ObservationPplusContext)
            or context.data_identity!=_data_identity(data)
            or context.identity!=_context_identity(context.data_identity,
                context.reference,context.observations,context.sim,
                context.graph,context.branches,context.metadata)):
        raise ValueError('Frozen P+ context does not match observations or numerical state')


def prepare_observation_pplus(data,composition,heavy,fpp,*,centric_singletons=None,
        normalization_context=None,fit_mask=None,scale_samples=256,
        phase_samples=720,kmin=.03,max_neighbors=600,max_per_target=96):
    """Freeze experimental constraints and conditional branches before updates.

    A fit_mask can hold out working observations for conditional scale
    diagnostics. Existing HA coordinates may already have seen those rows;
    such a diagnostic is not end-to-end independent cross-validation.
    """
    data_identity=_data_identity(data)
    n=len(data.f)
    single=np.zeros(n,bool) if centric_singletons is None else np.asarray(centric_singletons,dtype=bool)
    if single.shape!=(n,) or np.any(single&~data.centric):
        raise ValueError('Invalid centric singleton mask')
    if (not fpp or any(el not in fpp or not np.isfinite(fpp[el]) or fpp[el]<=0
                      for el in heavy.elements)):
        raise ValueError('Positive finite f-double-prime required')
    normalization=normalization_context or prepare_pplus_normalization(data,composition)
    if normalization.identity!=_normalization_identity(data,composition):
        raise ValueError('Frozen normalization does not match observations/composition')
    imaginary=structure_factors(data,heavy,fpp)
    alpha,fit=fit_phase_marginal_scale(data,imaginary,fit_mask=fit_mask,samples=scale_samples)
    reference=np.angle(imaginary)
    observations=sad_hl(reference,alpha*np.abs(imaginary),data.dano,data.sigdano)
    observations[single]=0.
    sim=_sim_field(data,heavy,composition,normalization.expected,normalization.statistics['scale'])
    graph=build_relations(data,normalization.e,normalization.kappa,kmin,max_neighbors,max_per_target)
    branches=prepare_conditional_branches(data,observations,reference,samples=phase_samples)
    for array in [reference,observations,sim]:array.setflags(write=False)
    for name in ['target','left','right','left_sign','right_sign','shift','kappa']:
        getattr(graph,name).setflags(write=False)
    metadata={'method':'observation_consistent_conditional_pplus_experimental',
        'status':'opt_in_research_not_validated_default',
        'observation_model':'DANO=2*alpha*abs(Q)*cos(phi-reference)+Gaussian(reported_sigmaD)',
        'uncertainty_policy':'Reported sigmaD unchanged; no AOE rescaling or sigma contraction',
        'anomalous_scale':float(alpha),'anomalous_scale_fit':fit,
        'scale_interpretation':'Effective scale for a fixed HA model; not measured f-double-prime or a site-ranking score',
        'centric_singletons_without_anomalous_likelihood':int(single.sum()),
        'partial_structure_field':'Existing Wilson-normalized Sim field; no f-prime/ADP uncertainty recalibration',
        'joint_bijvoet_likelihood':False,'fully_calibrated_SAD_model':False,
        'normalization':dict(normalization.statistics),'normalization_context_sha256':normalization.identity,
        'relations':graph.statistics,'phase_samples':phase_samples,
        'scale_samples_requested':scale_samples,
        'fit_mask_sha256':None if fit_mask is None else
            hashlib.sha256(np.asarray(fit_mask,dtype=bool).tobytes()).hexdigest(),
        'conditional_branch_shape':'Frozen experimental half-circle distributions; variational two-state weight update',
        'protein_reference_used':False,'external_phases_used':False}
    identity=_context_identity(data_identity,reference,observations,sim,graph,branches,metadata)
    return ObservationPplusContext(identity,data_identity,reference,
        observations,sim,graph,branches,metadata)


def observation_pplus_prior(data,context):
    """Export only the frozen experimental HL, with its own integrated moments."""
    _check_context(data,context)
    phase,fom=posterior_moments(data,context.observations,
        samples=context.metadata['phase_samples'])
    return PhaseResult(phase,fom,None,context.observations.copy(),{
        **context.metadata,'role':'Pure small-signal observation prior; no Sim/triplet posterior',
        'context_sha256':context.identity})


def phase_observation_pplus(data,context,*,cycles=5,triplet_weight=1.,damping=1.):
    """Update P+ branch masses while retaining the experimental branch shapes.

    The returned seed carries no HL. Persist the experimental prior through
    observation_pplus_prior, so an HL-only downstream reader cannot silently
    substitute the experimental prior for this updated map seed.
    The second seed moment is reported in metadata only through diagnostics;
    this API does not pretend that four HL exactly represent its branch mixture.
    """
    _check_context(data,context)
    result=iterate_conditional_pplus(data,context.branches,context.sim,context.graph,
        cycles=cycles,triplet_weight=triplet_weight,damping=damping)
    metadata={**context.metadata,'role':'P+ map seed; experimental HL retained separately',
        'context_sha256':context.identity,'cycles':result.history,
        'triplet_weight':float(triplet_weight),'damping':float(damping),
        'hl_definition':'No HL in this map seed; export the experimental prior separately',
        'second_seed_moment_mean_magnitude':float(np.mean(np.abs(result.second))),
        'iteration':result.metadata}
    return PhaseResult(result.phase,result.fom,result.probability,None,metadata)
