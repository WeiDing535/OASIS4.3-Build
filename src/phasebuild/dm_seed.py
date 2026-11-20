"""Explicit optional precision contract for the map-only DM seed.

The original P+ result and SAD HL remain unchanged. This limits insignificant
floating-point differences crossing discontinuous external DM operations; it
does not establish insensitivity to larger, scientifically meaningful changes.
"""
import numpy as np


def canonicalize(phase_degrees,fom):
    phase=np.asarray(phase_degrees,dtype=float);weight=np.asarray(fom,dtype=float)
    if phase.ndim!=1 or phase.shape!=weight.shape or not np.isfinite(phase).all() or not np.isfinite(weight).all() or np.any((weight<0)|(weight>1)):
        raise ValueError('DM seed needs matching finite phase/FOM vectors and FOM in [0,1]')
    rounded_phase=np.round(phase,3).astype(np.float32)
    rounded_fom=np.round(weight,6).astype(np.float32)
    return rounded_phase,rounded_fom,dict(enabled=True,phase_decimal_places=3,fom_decimal_places=6,
        max_phase_change_degrees=float(np.max(np.abs(rounded_phase.astype(float)-phase))) if phase.size else 0.,
        max_fom_change=float(np.max(np.abs(rounded_fom.astype(float)-weight))) if weight.size else 0.,
        scope='Map-only seed precision, original P+ phase file and raw SAD HL unchanged; no Rfree or reference-model input')
