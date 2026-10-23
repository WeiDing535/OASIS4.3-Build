# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Experimental local weighted NCS averaging with explicit Cartesian operators.

Operators must come from heavy atoms or an allowed partial model, never the
complete target structure. Automatic heavy-atom operator discovery is not yet
implemented: a missing/unsupported operator means NCS remains OFF.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from .native_density import smooth


def validate_operators(operators):
    result=[]
    if not isinstance(operators,(list,tuple)):raise ValueError('NCS operators must be a list')
    if any(not isinstance(op,dict) or 'rotation' not in op or 'translation' not in op for op in operators):raise ValueError('NCS operators require rotation and translation')
    if len(operators)>4:raise ValueError('At most four NCS operators supported')
    for op in operators:
        rotation=np.asarray(op['rotation'],dtype=float);translation=np.asarray(op['translation'],dtype=float)
        if rotation.shape!=(3,3) or translation.shape!=(3,) or not np.isfinite(rotation).all() or not np.isfinite(translation).all():raise ValueError('Invalid NCS operator dimensions or values')
        if not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-5) or not np.isclose(np.linalg.det(rotation),1,atol=1e-5):raise ValueError('NCS rotation must be proper and orthogonal')
        if np.allclose(rotation,np.eye(3),atol=1e-5) and np.linalg.norm(translation)<.1:raise ValueError('Identity NCS operator provides no information')
        result.append({'rotation':rotation,'translation':translation})
    return result


def resample(values,cell,rotation,translation):
    """Sample rho(R x + t), with R/t Cartesian and crystallographic PBC."""
    shape=np.array(values.shape);frac=np.indices(values.shape,dtype=float).reshape(3,-1)/shape[:,None]
    orth=np.array(cell.orth.mat);inv=np.array(cell.frac.mat)
    transformed=(inv@rotation@orth)@frac+(inv@translation)[:,None]
    # grid-wrap interpolates across the actual periodic voxel boundary.
    coords=(transformed%1)*shape[:,None]
    return map_coordinates(values,coords,order=1,mode='grid-wrap',prefilter=False).reshape(values.shape)


def prepare_ncs(grid,mask,kernel,operators,threshold=.7):
    """Only average locally correlated ordered density; fixed weights per cycle."""
    if not 0<threshold<1:raise ValueError('NCS correlation threshold must lie in (0,1)')
    ops=validate_operators(operators);rho=np.array(grid.array,dtype=float)
    a=smooth(rho,kernel);va=np.maximum(smooth(rho*rho,kernel)-a*a,0)
    plans=[];diagnostics=[]
    for op in ops:
        rotation=op['rotation'];translation=op['translation']
        for r,t in [(rotation,translation),(rotation.T,-rotation.T@translation)]:
            other=resample(rho,grid.unit_cell,r,t);b=smooth(other,kernel)
            vb=np.maximum(smooth(other*other,kernel)-b*b,0)
            correlation=(smooth(rho*other,kernel)-a*b)/np.sqrt(np.maximum(va*vb,1e-30))
            target_mask=resample(mask,grid.unit_cell,r,t)
            weight=np.clip((correlation-threshold)/(1-threshold),0,1)*mask*target_mask
            plans.append((r,t,weight))
            diagnostics.append({'active_fraction':float(np.mean(weight>.05)),'mean_weight':float(weight.mean())})
    return plans,diagnostics


def apply_ncs(rho,cell,plans):
    values=rho.copy();weights=np.ones_like(rho)
    for r,t,w in plans:
        values+=w*resample(rho,cell,r,t);weights+=w
    return values/weights
