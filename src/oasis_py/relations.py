# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Symmetry-aware three-phase relations without packed Miller indices."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class Relations:
    target: np.ndarray
    left: np.ndarray
    right: np.ndarray
    left_sign: np.ndarray
    right_sign: np.ndarray
    shift: np.ndarray
    kappa: np.ndarray
    nreflections: int
    statistics: dict

    def field(self,phase,fom):
        phase=np.asarray(phase);fom=np.asarray(fom)
        if phase.shape!=(self.nreflections,) or fom.shape!=phase.shape:
            raise ValueError('Phase/FOM shape mismatch')
        angle=(self.left_sign*phase[self.left]+self.right_sign*phase[self.right]+self.shift)
        weight=self.kappa*fom[self.left]*fom[self.right]
        real=np.bincount(self.target,weights=weight*np.cos(angle),minlength=self.nreflections)
        imag=np.bincount(self.target,weights=weight*np.sin(angle),minlength=self.nreflections)
        return real+1j*imag


def build_relations(data,e,kappa_scale,kmin=.03,max_neighbors=1000,
                    max_per_target=256,max_relations=2_000_000):
    """Build h=h'+h'' using symmetry and Friedel equivalents.

    Neighbor/per-target caps are explicit approximations. Metadata reports
    whether they removed candidates. max_neighbors=0/max_per_target=0 means
    no cap. A hard global limit raises rather than silently truncating.
    """
    e=np.asarray(e,dtype=float)
    n=len(data.f)
    if e.shape!=(n,) or np.any(e<0) or not np.all(np.isfinite(e)):
        raise ValueError('Invalid normalized amplitudes')
    if kmin<0 or kappa_scale<=0 or min(max_neighbors,max_per_target)<0:
        raise ValueError('Invalid relation thresholds')
    lookup={}
    for i,h in enumerate(data.hkl):
        for op in data.spacegroup.operations():
            equiv=tuple(op.apply_to_hkl(h.tolist()))
            shift=op.phase_shift(h.tolist())
            lookup.setdefault(equiv,(i,1,shift))
            lookup.setdefault(tuple(-x for x in equiv),(i,-1,-shift))
    keys=sorted(lookup,key=lambda h:(-e[lookup[h][0]],h))
    expanded=len(keys)
    capped_neighbors=bool(max_neighbors and len(keys)>max_neighbors)
    if max_neighbors: keys=keys[:max_neighbors]
    records=[];discarded=0
    for i,h in enumerate(data.hkl):
        hi,hj,hk=map(int,h)
        local=[]
        for key in keys:
            other=(hi-key[0],hj-key[1],hk-key[2])
            if other not in lookup:continue
            j,sj,tj=lookup[key];k,sk,tk=lookup[other]
            if i==j or i==k or j==k:continue
            kap=kappa_scale*e[i]*e[j]*e[k]
            if kap<kmin:continue
            local.append((i,j,k,sj,sk,tj+tk,kap))
        if max_per_target and len(local)>max_per_target:
            local.sort(key=lambda r:-r[-1]);discarded+=len(local)-max_per_target
            local=local[:max_per_target]
        records.extend(local)
        if len(records)>max_relations:
            raise ValueError('Relation memory budget exceeded; increase budget or set explicit caps')
    arr=np.asarray(records,dtype=float).reshape(-1,7)
    ints=[arr[:,i].astype(np.int64) for i in range(5)]
    stats={'relations':len(records),'expanded_hkl':expanded,'neighbors_used':len(keys),
           'neighbor_cap_applied':capped_neighbors,'discarded_by_per_target_cap':discarded,
           'kmin':kmin, 'approximate':capped_neighbors or discarded>0}
    return Relations(*ints,arr[:,5],arr[:,6],n,stats)
