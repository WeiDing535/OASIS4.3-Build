"""Explicit research handoff of a phase distribution to density modification.

P+ PHIB/FOM and its observation-only HL are different distributions. Choose
which to pass explicitly, then use this same persisted prior in either engine.
The fixed observation likelihood used by REFMAC remains a separate artifact.
"""
import numpy as np
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py.native_density import initial_prior
from oasis_py.hl_constraints import canonicalize_centric_hl
from .pplus_experimental import mtz_safe_bimodal_prior


def prepare_dm_prior(data,initial,*,representation):
    """Return consistent, float32-safe PHIB/FOM/HL for an explicit distribution.

    observation_hl intentionally excludes P+ branch updates. phase_posterior
    retains consistent posterior HL, or matches P+ PHIB/FOM by a von Mises
    distribution (centric rows are discrete); the latter loses doublet shape.
    This changes only the DM handoff, never the original observation artifact.
    """
    if representation=='observation_hl':
        if initial.hl is None:raise ValueError('Observation HL required')
        hl=np.array(initial.hl,dtype=float,copy=True);kind='original_observation_HL_no_Pplus_branch_updates';mismatch=None
    elif representation=='phase_posterior':
        hl,kind,mismatch=initial_prior(data,initial)
    else:raise ValueError('Choose observation_hl or phase_posterior explicitly')
    hl,canonical=canonicalize_centric_hl(data,hl)
    phase,fom=posterior_moments(data,hl)
    meta={'method':'explicit_dm_prior_research','role':'DM input distribution; not the REFMAC observation prior',
          'requested_representation':representation,'effective_distribution':kind,
          'source_method':initial.metadata.get('method'),'source_HL_moment_mismatch':mismatch,
          'centric_canonicalization':canonical,'source_observation_artifact_modified':False}
    result=PhaseResult(phase,fom,None,hl,meta)
    return mtz_safe_bimodal_prior(data,result)
