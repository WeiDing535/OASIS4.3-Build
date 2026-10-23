# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental alternating fixed-prior SAD density modification.

A and B are initial P+ / continuous SAD solutions on the SAME reflection
rows. Each DM round consumes the preceding map but replaces DM evidence and
uses exactly one fixed initial prior. This is not alternating triplet updates,
not a joint independent-data likelihood, and not a Parrot implementation.
"""
from .native_density import _density_modify


def density_modify_alternating(data, first, second, solvent_fraction,
                              reference_mtz, cycles=20, pattern='AB',
                              seed=20261016, max_grid_points=8_000_000,
                              cycle_callback=None):
    """Run A-only, B-only, AB or BA with identical DM round budgets.

    ``first`` and ``second`` must use data.hkl order. The first active prior
    also sets the starting phase and generic histogram simulation. A and AB
    therefore share initialization; B and BA share initialization. A callback
    receives copies for logging, and cannot change the running state. Return
    the FINAL round, never select the best round against a target model.
    """
    if pattern not in ('A','B','AB','BA'):
        raise ValueError('pattern must be A, B, AB or BA')
    return _density_modify(data,first,solvent_fraction,reference_mtz,
        cycles=cycles,seed=seed,max_grid_points=max_grid_points,
        _alternate_initial=second,_pattern=pattern,_cycle_callback=cycle_callback)
