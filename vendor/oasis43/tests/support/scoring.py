import numpy as np

def metrics(data,phase,fom,truth,origin_shift=None):
    """Evaluate phases, with an optional explicitly established origin shift.

    origin_shift is model-origin minus reference-origin in fractional units.
    It must come from independently matched coordinates/provenance, never
    selection of whichever shift gives the best reference score. This does
    not search origins or hands and never changes solver output.
    """
    if origin_shift is not None:
        shift=np.asarray(origin_shift,dtype=float)
        if shift.shape!=(3,) or not np.all(np.isfinite(shift)):
            raise ValueError('Invalid fractional origin shift')
        for op in data.spacegroup.operations():
            delta=(np.asarray(op.apply_to_xyz(shift.tolist()))
                   -np.asarray(op.apply_to_xyz([0.,0.,0.]))-shift)
            if np.max(np.abs((delta+.5)%1-.5))>1e-7:
                raise ValueError('Origin translation changes the space-group setting')
        phase=np.asarray(phase)-2*np.pi*(data.hkl@shift)
    delta=np.angle(np.exp(1j*(phase-np.angle(truth))))
    good=np.abs(truth)>1e-7
    weight=data.f[good]**2
    # Equivalent same-amplitude reciprocal-space correlation. No origin/hand search.
    m=fom[good];d=delta[good]
    corr=np.sum(weight*m*np.cos(d))/np.sqrt(np.sum(weight*m*m)*np.sum(weight))
    amp=np.abs(truth);scale=np.dot(amp,data.f)/np.dot(amp,amp)
    return dict(phase_mae_deg=float(np.mean(np.abs(d))*180/np.pi),
        weighted_phase_mae_deg=float(np.average(np.abs(d),weights=weight)*180/np.pi),
        phase_cosine=float(np.mean(np.cos(d))),fourier_coefficient_cc=float(corr),
        mean_fom=float(np.mean(fom)),reference_amplitude_r=float(np.sum(np.abs(data.f-scale*amp))/np.sum(data.f)))
