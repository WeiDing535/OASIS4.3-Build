# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Explicit adapter to external CCP4 Parrot density modification."""
from pathlib import Path
import json,subprocess,time
import gemmi,numpy as np
from .data import read_mtz,write_mtz
from .sad import posterior_moments


def density_modify_parrot(executable,hklin,out_dir,solvent_fraction,reference_mtz,
                           reference_pdb,cycles=20,representation='hl',timeout=300,
                           canonical_centric=False):
    """Use a generic reference histogram; never accept a target complete model.

    Input uses the native FP/SIGFP plus HLA..HLD or PHIB/FOM schema.
    The generic reference structure must be unrelated to the target.
    """
    if not 0<solvent_fraction<1 or cycles<1 or timeout<=0:
        raise ValueError('Invalid solvent fraction, cycle count or timeout')
    if representation not in {'hl','phifom'}:raise ValueError('representation must be hl or phifom')
    if canonical_centric and representation!='hl':
        raise ValueError('Centric HL canonicalization requires representation=hl')
    for path in [executable,hklin,reference_mtz,reference_pdb]:
        if not Path(path).is_file():raise ValueError(f'File not found: {path}')
    out=Path(out_dir).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Use an empty density-modification output directory')
    data,input_filter=read_mtz(hklin,'dmr',{'f':'FP','sigf':'SIGFP'})
    out.mkdir(parents=True,exist_ok=True);raw=out/'parrot.mtz'
    effective_input=Path(hklin).resolve();canonical_meta={'enabled':False}
    if canonical_centric:
        from .hl_constraints import canonicalize_centric_hl
        original=gemmi.read_mtz_file(str(hklin));original.ensure_asu()
        lookup={tuple(h):i for i,h in enumerate(original.make_miller_array())}
        order=[lookup[tuple(h)] for h in data.hkl]
        try:prior=np.column_stack([original.column_with_label('HL'+k).array[order] for k in 'ABCD'])
        except AttributeError:raise ValueError('Centric canonicalization requires HLA, HLB, HLC, HLD')
        canonical,canonical_meta=canonicalize_centric_hl(data,prior)
        phase,fom=posterior_moments(data,canonical)
        effective_input=out/'canonical_input.mtz';write_mtz(effective_input,data,phase,fom,canonical)
        canonical_meta['enabled']=True
    command=[str(Path(executable).resolve()),'-mtzin-ref',str(Path(reference_mtz).resolve()),
        '-pdbin-ref',str(Path(reference_pdb).resolve()),
        '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
        '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
        '-mtzin',str(effective_input),'-colin-fo','/*/*/[FP,SIGFP]',
        '-mtzout',str(raw),'-colout','DM','-solvent-flatten','-histogram-match',
        '-solvent-content',str(solvent_fraction),'-cycles',str(cycles)]
    command+=['-colin-hl','/*/*/[HLA,HLB,HLC,HLD]'] if representation=='hl' else ['-colin-phifom','/*/*/[PHIB,FOM]']
    start=time.perf_counter()
    with (out/'parrot.log').open('w') as log:
        try:run=subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,timeout=timeout,check=False)
        except subprocess.TimeoutExpired:raise RuntimeError(f'Parrot exceeded {timeout}s; see {out}/parrot.log')
    if run.returncode or not raw.exists():raise RuntimeError(f'Parrot failed with exit {run.returncode}; see {out}/parrot.log')
    mtz=gemmi.read_mtz_file(str(raw));mtz.ensure_asu()
    if mtz.spacegroup.operations()!=data.spacegroup.operations() or not np.allclose(mtz.cell.parameters,data.cell.parameters,rtol=1e-5,atol=1e-3):
        raise RuntimeError('Parrot output changed cell or space-group setting')
    lookup={tuple(h):i for i,h in enumerate(mtz.make_miller_array())}
    try:
        ix=[lookup[tuple(h)] for h in data.hkl]
        hl=np.column_stack([mtz.column_with_label('DM.ABCD.'+k).array[ix] for k in 'ABCD']).astype(float)
    except (KeyError,AttributeError):raise RuntimeError('Missing Parrot output reflections or HL coefficients')
    phase,fom=posterior_moments(data,hl)
    write_mtz(out/'phases.mtz',data,phase,fom,hl)
    meta={'backend':'CCP4 Parrot external density modification','command':command,
          'solvent_fraction':solvent_fraction,'cycles':cycles,'input_representation':representation,
          'input_filter':input_filter,'target_complete_model_used':False,
          'mtz':str(out/'phases.mtz'),'raw_mtz':str(raw),
          'mean_fom':float(fom.mean()),'seconds':time.perf_counter()-start}
    if canonical_centric:meta['centric_HL_canonicalization']=canonical_meta
    (out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta
