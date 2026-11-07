"""CCP4 Parrot density-modification adapter with anisotropic scaling correction.

This is an optional external engine: it never runs inside the bundled runtime and
never receives SHELXE phases, reference coordinates or free-set reflections. Only
the persisted original observation likelihood (HLA..HLD) and, in later feedback
cycles, a work-only partial model's Fcalc are supplied. The generated map is a
work-only phase estimate; the free set is excluded before Parrot is called.
"""
from pathlib import Path
import os
import shutil
import subprocess
import gemmi
import numpy as np

DEFAULT_REFERENCE='lib/data/reference_structures/reference-1tqw'

CANDIDATES=('cparrot',)


def locate(explicit=None, runtime=None):
    """Return a usable cparrot path, else raise with a clear message."""
    if explicit:
        path=Path(explicit).expanduser()
        if path.is_file() and os.access(path,os.X_OK):return str(path)
        raise ValueError('Configured Parrot executable is not runnable: '+str(path))
    if runtime:
        candidate=Path(runtime)/'bin'/'cparrot'
        if candidate.is_file() and os.access(candidate,os.X_OK):return str(candidate)
    found=shutil.which('cparrot')
    if found:return found
    for base in ('/Applications/ccp4-9','/opt/ccp4','/usr/local/ccp4'):
        candidate=Path(base)/'bin'/'cparrot'
        if candidate.is_file() and os.access(candidate,os.X_OK):return str(candidate)
    raise ValueError('Parrot (cparrot) was not found; set --dm-executable or use the internal engine')


def _environment(executable, scratch):
    binary=Path(executable).resolve()
    base=binary.parent.parent
    env=dict(os.environ)
    env.update({'PATH':str(binary.parent)+':'+env.get('PATH','/usr/bin:/bin'),
                'CCP4':str(base),'CCP4_SCR':str(scratch),
                'CLIBD':str(base/'lib/data'),'CLIBD_MON':str(base/'lib/data/monomers'),
                'SYMINFO':str(base/'lib/data/syminfo.lib')})
    Path(scratch).mkdir(parents=True,exist_ok=True)
    return env


def reference_files(executable, runtime=None):
    """Generic unrelated reference histogram (1tqw), never a target structure."""
    if runtime:
        stem=Path(runtime)/DEFAULT_REFERENCE
        if stem.with_suffix('.mtz').is_file() and stem.with_suffix('.pdb').is_file():
            return stem.with_suffix('.mtz'),stem.with_suffix('.pdb')
    base=Path(executable).resolve().parent.parent
    for stem in (base/DEFAULT_REFERENCE,Path('/Applications/ccp4-9')/DEFAULT_REFERENCE):
        if stem.with_suffix('.mtz').is_file():return stem.with_suffix('.mtz'),stem.with_suffix('.pdb')
    raise ValueError('Generic Parrot reference histogram is unavailable')


def _hl_columns(path):
    m=gemmi.read_mtz_file(str(path));labels=m.column_labels()
    for prefix in ('DM','PIRATE','PARROT'):
        quartet=[prefix+'.ABCD.'+k for k in 'ABCD']
        if all(c in labels for c in quartet):return quartet
    quartet=['HL'+k for k in 'ABCD']
    if not all(c in labels for c in quartet):raise ValueError('No HL quartet in '+str(path))
    return quartet


def to_phifom(source,destination):
    """Convert a Parrot HL quartet to PHIB(deg)/FOM using von Mises moments."""
    from scipy.special import i0e,i1e
    m=gemmi.read_mtz_file(str(source));m.ensure_asu();a=np.array(m,copy=True)
    labels=m.column_labels();quartet=_hl_columns(source)
    ix=[labels.index(c) for c in quartet]
    values=a[:,ix]
    if not np.isfinite(values).all():raise ValueError('Non-finite HL coefficients in '+str(source))
    kappa=np.hypot(values[:,0],values[:,1])
    phase=np.degrees(np.arctan2(values[:,1],values[:,0]))%360.
    fom=np.where(kappa>30,np.ones_like(kappa),i1e(np.clip(kappa,0,30))/np.maximum(i0e(np.clip(kappa,0,30)),1e-300))
    out=gemmi.Mtz(with_base=True);out.spacegroup=m.spacegroup;out.set_cell_for_all(m.cell);out.add_dataset('parrot_dm')
    for label,kind in [('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHIB','P'),('FOM','W'),('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')]:
        out.add_column(label,kind)
    fp=np.asarray(m.column_with_label('FP')) if 'FP' in labels else np.asarray(m.column_with_label('F'))
    sig=np.asarray(m.column_with_label('SIGFP')) if 'SIGFP' in labels else np.asarray(m.column_with_label('SIGF'))
    out.set_data(np.column_stack([a[:,:3],fp,sig,fp*fom,phase,phase,fom,values]).astype(np.float32))
    out.write_to_file(str(destination))
    return {'mean_fom':float(fom.mean()),'columns':quartet}


def run(work_mtz, hl_mtz, output_mtz, solvent_fraction, *, cycles=20, mask_radius=3.0,
        anisotropy=True, model_pdb=None, executable=None, runtime=None, metadata=None,
        scratch=None, timeout=3600, phenix_bin=None):
    """Run Parrot on the work-only data and write the pipeline map schema.

    ``work_mtz`` must contain FP/SIGFP and FreeR_flag; free reflections are never
    passed to Parrot. ``model_pdb`` is an optional work-only partial model whose
    Fcalc is combined with the fixed experimental likelihood.
    """
    if not 0<solvent_fraction<1:raise ValueError('Solvent fraction must lie between 0 and 1')
    executable=locate(executable,runtime)
    work=gemmi.read_mtz_file(str(work_mtz));work.ensure_asu();labels=work.column_labels()
    fp_label='FP' if 'FP' in labels else 'F'
    sig_label='SIGFP' if 'SIGFP' in labels else 'SIGF'
    hl_source=hl_mtz or work_mtz
    quartet=_hl_columns(hl_source)
    # stage work-only input carrying the experimental HL quartet
    src=gemmi.read_mtz_file(str(hl_source));src.ensure_asu();src_labels=src.column_labels()
    lookup={tuple(int(v) for v in row[:3]):row for row in np.array(src)}
    rows=[]
    for row in np.array(work):
        h=tuple(int(v) for v in row[:3]);other=lookup.get(h)
        if other is None:continue
        rows.append(list(h)+[row[labels.index(fp_label)],row[labels.index(sig_label)]]+[other[src_labels.index(c)] for c in quartet])
    if not rows:raise ValueError('No work reflections shared with the HL prior')
    staged=Path(output_mtz).parent/'parrot_input.mtz'
    m=gemmi.Mtz(with_base=True);m.spacegroup=work.spacegroup;m.set_cell_for_all(work.cell);m.add_dataset('parrot_work')
    m.add_column('FP','F');m.add_column('SIGFP','Q')
    for k in 'ABCD':m.add_column('HL'+k,'A')
    m.set_data(np.asarray(rows,np.float32));m.write_to_file(str(staged))
    reference_mtz,reference_pdb=reference_files(executable,runtime)
    raw=Path(output_mtz).parent/'parrot_raw.mtz'
    command=[executable,'-mtzin-ref',str(reference_mtz),'-pdbin-ref',str(reference_pdb),
        '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
        '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
        '-mtzin',str(staged),'-colin-fo','/*/*/[FP,SIGFP]','-colin-hl','/*/*/[HLA,HLB,HLC,HLD]',
        '-solvent-content','%.10g'%solvent_fraction,'-solvent-mask-filter-radius','%.4g'%mask_radius,
        '-cycles',str(int(cycles)),'-solvent-flatten','-histogram-match',
        '-mtzout',str(raw),'-colout','DM']
    if anisotropy:command.append('-anisotropy-correction')
    if model_pdb is not None:
        if not Path(model_pdb).is_file():raise ValueError('Feedback model is missing')
        command+=['-colin-fc','/*/*/[FC,PHIC]']
        work_with_fc=Path(output_mtz).parent/'parrot_with_model.mtz'
        _attach_fc(staged,model_pdb,work_with_fc,phenix_bin)
        command[command.index('-mtzin')+1]=str(work_with_fc)
    log_path=Path(output_mtz).parent/'parrot.log'
    with log_path.open('w') as log:
        completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,
            env=_environment(executable,scratch or Path(output_mtz).parent/'parrot_scratch'),timeout=timeout)
    if completed.returncode or not raw.is_file():
        raise RuntimeError('Parrot failed with exit %s; see %s'%(completed.returncode,log_path))
    info=to_phifom(raw,output_mtz)
    return {'engine':'CCP4 Parrot','executable':executable,'solvent_fraction':float(solvent_fraction),
            'mask_filter_radius':float(mask_radius),'anisotropy_correction':bool(anisotropy),
            'cycles':int(cycles),'partial_model_used':model_pdb is not None,
            'work_reflections':len(rows),'free_reflections_passed':0,
            'phase_source':'original observation HL' if model_pdb is None else 'original HL combined with work-only model Fcalc',
            **(metadata or {}), **info}


def _attach_fc(staged, model_pdb, destination, phenix_bin=None):
    """Attach bulk-solvent-corrected FC/PHIC to the staged MTZ.

    A Phenix ``fmodel`` pass is preferred when an explicit binary directory is
    supplied: it applies bulk-solvent scaling without moving the model, which is
    the numerically stable option for a partial trace. Otherwise a single REFMAC
    cycle provides the same columns with the bundled runtime.
    """
    if phenix_bin:
        return _attach_fc_phenix(staged,model_pdb,destination,phenix_bin)
    from .refinement import refmac_refine
    source=gemmi.read_mtz_file(str(staged));labels=source.column_labels()
    fp_label='FP' if 'FP' in labels else 'F'
    sig_label='SIGFP' if 'SIGFP' in labels else 'SIGF'
    work=Path(destination).parent/('refmac_fc_'+Path(destination).stem)
    if work.exists():raise ValueError('REFMAC Fcalc workspace already exists: '+str(work))
    native=Path(destination).parent/('native_'+Path(destination).stem+'.mtz')
    m=gemmi.read_mtz_file(str(staged));copy=gemmi.Mtz(with_base=True)
    copy.spacegroup=m.spacegroup;copy.set_cell_for_all(m.cell);copy.add_dataset('native')
    for label,kind in [('F','F'),('SIGF','Q')]:copy.add_column(label,kind)
    copy.add_column('FreeR_flag','I')
    data=np.column_stack([np.array(m)[:,:3],np.array(m)[:,labels.index(fp_label)],np.array(m)[:,labels.index(sig_label)],
        (np.arange(m.nreflections)%20!=0).astype(float)])
    copy.set_data(data.astype(np.float32));copy.write_to_file(str(native))
    refined=refmac_refine(native,model_pdb,work,cycles=1)
    out=gemmi.read_mtz_file(str(work/'refined.mtz'))
    out_labels=out.column_labels()
    lookup={tuple(int(v) for v in row[:3]):row for row in np.array(out)}
    fcs=np.zeros(source.nreflections);phis=np.zeros(source.nreflections)
    for i,row in enumerate(np.array(source)):
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        fcs[i]=match[out_labels.index('FC')];phis[i]=match[out_labels.index('PHIC')]
    base=np.array(source,copy=True)
    source.add_column('FC','F');source.add_column('PHIC','P')
    source.set_data(np.column_stack([base,fcs,phis]).astype(np.float32))
    source.write_to_file(str(destination))
    return refined


def _attach_fc_phenix(staged, model_pdb, destination, phenix_bin):
    """FC/PHIC from phenix.fmodel; no coordinate movement, bulk solvent included."""
    binary=Path(phenix_bin)/'phenix.fmodel'
    if not binary.is_file():raise ValueError('phenix.fmodel is missing in '+str(phenix_bin))
    source=gemmi.read_mtz_file(str(staged));labels=source.column_labels()
    base=np.array(source,copy=True)
    fc_mtz=Path(destination).with_suffix('.fmodel.mtz')
    completed=subprocess.run([str(binary),str(model_pdb),'high_resolution=%.4f'%source.resolution_high(),
                   'output.file_name='+str(fc_mtz)],
                   stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    if completed.returncode or not fc_mtz.is_file():
        raise RuntimeError('phenix.fmodel Fcalc pass failed for '+str(model_pdb))
    model_mtz=gemmi.read_mtz_file(str(fc_mtz));ml=model_mtz.column_labels();ma=np.array(model_mtz)
    lookup={tuple(int(v) for v in row[:3]):row for row in ma}
    fcs=np.zeros(source.nreflections);phis=np.zeros(source.nreflections)
    for i,row in enumerate(base):
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        fcs[i]=match[ml.index('FMODEL')];phis[i]=match[ml.index('PHIFMODEL')]
    source.add_column('FC','F');source.add_column('PHIC','P')
    source.set_data(np.column_stack([base,fcs,phis]).astype(np.float32))
    source.write_to_file(str(destination))
    return {'engine':'phenix.fmodel','binary':str(binary)}
