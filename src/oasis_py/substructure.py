# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental sparse-atom dual-space search (not a port of SAPI)."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import gemmi
from scipy.ndimage import maximum_filter
from .data import Atoms,structure_factors


@dataclass
class SubstructureResult:
    atoms: Atoms
    correlation: float
    trials: list[dict]
    metadata: dict


def coefficient_map(data,values,sample_rate=3.):
    asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,
                           np.ascontiguousarray(data.hkl,dtype=np.int32),
                           np.ascontiguousarray(values,dtype=np.complex64))
    return asu.transform_f_phi_to_map(sample_rate=sample_rate)


def pick_peaks(grid,count,min_distance=1.5):
    """Pick ASU-distinct maxima using periodic crystallographic distances."""
    if count<1 or min_distance<0:raise ValueError('Invalid peak-search controls')
    rho=np.asarray(grid.array)
    maxima=(rho==maximum_filter(rho,size=3,mode='wrap')) & (rho>0)
    coords=np.argwhere(maxima)
    if not len(coords):raise ValueError('No positive map maxima')
    values=rho[tuple(coords.T)]
    # Examine more than count peaks because symmetry creates whole orbits.
    order=np.argsort(-values,kind='stable')
    selected=[];heights=[]
    for row in order:
        frac=coords[row]/np.asarray(rho.shape)
        duplicate=False
        for previous in selected:
            for op in grid.spacegroup.operations():
                equivalent=op.apply_to_xyz(previous.tolist())
                distance=grid.unit_cell.find_nearest_pbc_image(
                    gemmi.Fractional(*frac),gemmi.Fractional(*equivalent)).dist()
                if distance<min_distance:
                    duplicate=True;break
            if duplicate:break
        if not duplicate:
            selected.append(frac);heights.append(float(values[row]))
            if len(selected)==count:break
    if len(selected)<count:
        raise ValueError(f'Only {len(selected)} separated ASU peaks for {count} requested sites')
    return np.asarray(selected),np.asarray(heights)


def _correlation(a,b):
    if np.std(a)<1e-12 or np.std(b)<1e-12:return 0.
    return float(np.corrcoef(a,b)[0,1])


def search_native(data,nsites,element='S',trials=8,cycles=50,seed=1,
                  min_distance=1.5,sample_rate=3.,max_grid_points=8_000_000):
    """Sparse peak-projection search using |DANO| as an FA proxy.

    The anomalous difference is NOT an exact substructure amplitude. No
    claim of SHELXD/HySS parity is made. Both hands must be evaluated by
    downstream experimental phasing/density modification. This routine is
    useful for algorithm experiments and exact-amplitude synthetic tests.
    """
    if data.dano is None:raise ValueError('Substructure search needs DANO')
    if min(nsites,trials,cycles)<1:raise ValueError('Positive sites/trials/cycles required')
    if gemmi.Element(element).atomic_number==0:raise ValueError('Unknown heavy element')
    # Remove radial attenuation by equal-population shell normalization.
    proxy=np.abs(data.dano).astype(float)
    if np.max(proxy)<=0:raise ValueError('No anomalous signal')
    for ix in np.array_split(np.argsort(data.stol2),min(12,max(1,len(proxy)//50))):
        proxy[ix]/=max(np.sqrt(np.mean(proxy[ix]**2)),1e-12)
    asu=gemmi.ComplexAsuData(data.cell,data.spacegroup,data.hkl,
                           np.ascontiguousarray(proxy,dtype=np.complex64))
    # Ask Gemmi for symmetry-compatible grid sizes before allocating a map.
    size=asu.get_size_for_hkl(sample_rate=sample_rate)
    if np.prod(size)>max_grid_points:
        raise ValueError(f'FFT map requires {np.prod(size)} points; exceeds explicit memory budget')
    rng=np.random.default_rng(seed);best=None;results=[]
    for trial in range(trials):
        atoms=Atoms(rng.random((nsites,3)),[element]*nsites,np.ones(nsites),np.zeros(nsites))
        best_trial=-np.inf
        for cycle in range(cycles):
            calculated=structure_factors(data,atoms,chunk=256)
            values=proxy*np.exp(1j*np.angle(calculated))
            grid=coefficient_map(data,values,sample_rate)
            try:
                xyz,heights=pick_peaks(grid,nsites,min_distance)
            except ValueError:
                break
            # Unit occupancies avoid fitting an FA proxy as if it were an
            # experimentally measured substructure structure factor.
            atoms=Atoms(xyz,[element]*nsites,np.ones(nsites),np.zeros(nsites))
            calculated=structure_factors(data,atoms,chunk=256)
            # Match the radial normalization used for the target proxy.
            fc=np.abs(calculated)
            for ix in np.array_split(np.argsort(data.stol2),min(12,max(1,len(proxy)//50))):
                fc[ix]/=max(np.sqrt(np.mean(fc[ix]**2)),1e-12)
            cc=_correlation(proxy,fc)
            if cc>best_trial:best_trial=cc
            if best is None or cc>best[0]:best=(cc,atoms,trial+1,cycle+1)
        results.append({'trial':trial+1,'best_proxy_correlation':float(best_trial) if np.isfinite(best_trial) else None})
    if best is None:raise ValueError('All substructure starts failed to yield a valid peak set')
    cc,atoms,trial,cycle=best
    return SubstructureResult(atoms,cc,results,
        {'status':'experimental_not_SAPI_equivalent','algorithm':'sparse_peak_projection',
         'amplitude_proxy':'absolute_Bijvoet_difference','best_trial':trial,'best_cycle':cycle,
         'seed':seed,'grid_size':list(size),'hand_resolved':False})


def invert_hand(atoms,spacegroup):
    """Apply crystallographic change-of-hand operation, including origin shift."""
    op=spacegroup.change_of_hand_op()
    other=Atoms([op.apply_to_xyz(p.tolist()) for p in atoms.fractional],
                list(atoms.elements),atoms.occupancy.copy(),atoms.b_iso.copy())
    ops=spacegroup.operations();ops.change_basis_forward(op)
    other_sg=gemmi.find_spacegroup_by_ops(ops)
    if other_sg is None:raise ValueError('Could not identify changed-hand space group')
    return other,other_sg


def write_atoms(path,atoms,data,spacegroup=None):
    st=gemmi.Structure();st.cell=data.cell;st.spacegroup_hm=(spacegroup or data.spacegroup).hm
    model=gemmi.Model('1');chain=gemmi.Chain('H')
    for i,(xyz,el,occ,b) in enumerate(zip(atoms.fractional,atoms.elements,atoms.occupancy,atoms.b_iso),1):
        residue=gemmi.Residue();residue.name=el.upper();residue.seqid=gemmi.SeqId(i,' ')
        atom=gemmi.Atom();atom.name=el.upper();atom.element=gemmi.Element(el)
        atom.pos=data.cell.orthogonalize(gemmi.Fractional(*xyz));atom.occ=float(occ);atom.b_iso=float(b)
        residue.add_atom(atom);chain.add_residue(residue)
    model.add_chain(chain);st.add_model(model);st.write_pdb(str(path))
