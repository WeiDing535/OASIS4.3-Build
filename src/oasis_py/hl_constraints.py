# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Exact symmetry constraints on HL representations; no empirical correction."""
import numpy as np
from .clipper_likelihood import reflection_classes


def canonicalize_centric_hl(data,hl):
    """Remove unidentifiable centric HL terms, preserving both phase probabilities.

    On allowed phases theta0 and theta0+pi, C cos(2theta)+D sin(2theta)
    is constant and the first-harmonic component perpendicular to theta0
    is zero. Keeping u=A cos(theta0)+B sin(theta0) as the only coefficient
    preserves every centric likelihood ratio exactly, up to roundoff.

    This can prevent older integrators from underflowing when their global
    normalization bound is dominated by a forbidden phase direction. It
    neither calibrates FOM nor changes acentric coefficients or phase data.
    """
    a=np.array(hl,dtype=float,copy=True)
    if a.shape!=(len(data.f),4) or not np.all(np.isfinite(a)):
        raise ValueError('Expected finite HL coefficients with shape (n,4)')
    _,centric,theta=reflection_classes(data)
    c=np.cos(theta[centric]);s=np.sin(theta[centric])
    u=a[centric,0]*c+a[centric,1]*s
    a[centric,0]=u*c;a[centric,1]=u*s;a[centric,2:]=0.
    return a,dict(method='exact_centric_HL_equivalent_representation',
        centric_rows=int(centric.sum()),acentric_rows_unchanged=True,
        phase_likelihood_ratios_preserved=True,empirical_FOM_correction=False,
        target_reference_used=False)
