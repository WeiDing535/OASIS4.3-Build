"""Work-only trace/phase-feedback loop on top of external density modification.

Each cycle: density-modify the fixed original observation likelihood (optionally
combined with the previous cycle's work-only partial model) -> trace a partial
model from that map -> attach the model's bulk-solvent-corrected Fcalc to the same
fixed likelihood -> repeat. Nothing from the reference structure, SHELXE, or the
free set enters the loop; model acceptance is fixed in advance and never uses R-free.
"""
from pathlib import Path
import json
import shutil
import subprocess
import gemmi
import numpy as np

from .data import fingerprint
from .runtime import runtime_environment


def _cbuccaneer(runtime):
    path=Path(runtime)/'bin'/'cbuccaneer'
    if not path.is_file():raise ValueError('Bundled cbuccaneer is missing')
    return str(path)


def trace(work_mtz, phase_mtz, sequence, runtime, output_dir, *, resolution=2.0,
          cycles=5, jobs=1, timeout=3600, progress=print):
    """Trace a partial model from a work-only map and report what was built."""
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    reference=Path(runtime)/'lib/data/reference_structures/reference-1tqw'
    command=[_cbuccaneer(runtime),'-mtzin',str(work_mtz),
        '-colin-fo','/*/*/[FP,SIGFP]','-colin-phifom','/*/*/[PHIB,FOM]',
        '-colin-free','/*/*/[FreeR_flag]','-seqin',str(sequence),
        '-mtzin-ref',str(reference.with_suffix('.mtz')),'-pdbin-ref',str(reference.with_suffix('.pdb')),
        '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
        '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
        '-find','-grow','-join','-link','-sequence','-correct','-filter','-prune','-rebuild',
        '-cycles',str(int(cycles)),'-jobs',str(int(jobs)),'-resolution','%.3f'%resolution,
        '-new-residue-name','UNK','-pdbout',str(out/'model.pdb')]
    env=runtime_environment(Path(runtime),scratch=out/'scratch')
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    with (out/'buccaneer.log').open('w') as log:
        completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,
                                 env=env,timeout=timeout)
    if completed.returncode or not (out/'model.pdb').is_file():
        raise RuntimeError('Partial-model tracing failed; see '+str(out/'buccaneer.log'))
    residues=0
    for line in (out/'buccaneer.log').read_text(errors='ignore').splitlines():
        if 'residues were built' in line:
            try:residues=int(line.split()[0])
            except ValueError:pass
    fragments=0
    for line in (out/'buccaneer.log').read_text(errors='ignore').splitlines():
        if 'residues were built' in line and 'fragments' in line:
            try:fragments=int(line.split('in')[1].split('fragments')[0])
            except (ValueError,IndexError):pass
    return {'model':str(out/'model.pdb'),'residues':residues,'fragments':fragments,
            'resolution':float(resolution),'cycles':int(cycles),
            'command_sha256':fingerprint(out/'model.pdb')['sha256']}


def merge_phase_map(work_mtz, dm_mtz, destination, *, free_value=0):
    """Write the pipeline's 03_dm schema (F,SIGF,FWT,PHWT,PHIB,FOM,HL) work-only."""
    work=gemmi.read_mtz_file(str(work_mtz));work.ensure_asu()
    dm=gemmi.read_mtz_file(str(dm_mtz));dm.ensure_asu()
    work_labels=work.column_labels();dm_labels=dm.column_labels()
    fp_label,sig_label=_amp_labels(work_labels)
    lookup={tuple(int(v) for v in row[:3]):row for row in np.array(dm)}
    quartet=['HL'+k for k in 'ABCD']
    rows=[]
    for row in np.array(work):
        if row[work_labels.index('FreeR_flag')]==free_value:continue
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        f=row[work_labels.index(fp_label)];sig=row[work_labels.index(sig_label)]
        values=[match[dm_labels.index(c)] for c in quartet]
        if not np.isfinite(values).all():continue
        kappa=float(np.hypot(values[0],values[1]))
        if kappa<=0:fom=0.0
        elif kappa>30:fom=1.0
        else:
            from scipy.special import i0e,i1e
            fom=float(i1e(kappa)/max(i0e(kappa),1e-300))
        phase=float(np.degrees(np.arctan2(values[1],values[0]))%360.)
        rows.append([*row[:3],f,sig,f*fom,phase,phase,fom,*values])
    m=gemmi.Mtz(with_base=True);m.spacegroup=work.spacegroup;m.set_cell_for_all(work.cell);m.add_dataset('dm')
    for label,kind in [('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('PHIB','P'),('FOM','W'),
                       ('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')]:
        m.add_column(label,kind)
    if not rows:raise ValueError('No work reflections survived the DM merge')
    m.set_data(np.asarray(rows,np.float32));m.write_to_file(str(destination))
    return {'work_rows':len(rows)}


def run(work_mtz, prior_mtz, sequence, runtime, output, *, solvent_fraction, cycles=20,
        feedback_cycles=1, trace_resolution=2.0, mask_radius=3.0, anisotropy=True,
        executable=None, trace_cycles=5, progress=print, phenix_bin=None, fcalc='fc_all'):
    """Iterate external DM and work-only tracing; persist every cycle."""
    from . import parrot_dm
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    history=[]
    model=None
    merged=out/'work_only.mtz'
    merge_work(work_mtz,merged)
    dm_path=None
    for cycle in range(max(1,feedback_cycles+1)):
        folder=out/('cycle_%02d'%(cycle+1));folder.mkdir(exist_ok=True)
        dm_path=folder/'dm.mtz'
        info=parrot_dm.run(merged,prior_mtz,dm_path,solvent_fraction,cycles=cycles,
            mask_radius=mask_radius,anisotropy=anisotropy,model_pdb=model,executable=executable,
            runtime=runtime,scratch=folder/'scratch',phenix_bin=phenix_bin,fcalc=fcalc)
        row={'cycle':cycle+1,'dm':str(dm_path),'dm_info':info}
        if cycle<feedback_cycles:
            trace_in=folder/'trace_input.mtz'
            phase_input(merged,dm_path,trace_in)
            traced=trace(trace_in,dm_path,sequence,runtime,folder/'trace',
                         resolution=trace_resolution,cycles=trace_cycles,progress=progress)
            row['trace']=traced
            if traced['residues']<=0:
                # Nothing was traced, so there is no Fcalc to feed into the next
                # density-modification cycle.  Stop the feedback here and keep
                # this cycle's DM map: crashing on an empty REFMAC input would
                # throw away the map the user already has.
                row['feedback_stopped']='trace built no residues, so no Fcalc could be fed back; this cycle\'s DM map is kept'
                progress('feedback cycle %d: no residues traced, stopping the Fcalc feedback (keeping this cycle\'s DM map)'%(cycle+1))
                history.append(row)
                break
            model=traced['model']
            progress('feedback cycle %d: %d residues in %d fragments'%(cycle+1,traced['residues'],traced['fragments']))
        history.append(row)
    final=out/'map_phases.mtz'
    merge_phase_map(merged,dm_path,final)
    report={'engine':'parrot-trace-feedback','history':history,'final_map':str(final),
        'work_only':True,'free_reflections_in_map':0,'reference_coordinates_used':False,
        'shelxe_phases_used':False,'rfree_used_for_selection':False,
        'inputs':{k:fingerprint(v) for k,v in [('work',work_mtz),('prior',prior_mtz),('sequence',sequence)]},
        'parameters':{'solvent_fraction':solvent_fraction,'dm_cycles':cycles,
            'feedback_cycles':feedback_cycles,'trace_resolution':trace_resolution,
            'mask_filter_radius':mask_radius,'anisotropy_correction':anisotropy}}
    (out/'feedback_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def phase_input(work_mtz, dm_mtz, destination):
    """Work rows with FP,SIGFP,FreeR_flag plus DM PHIB/FOM for the tracer."""
    work=gemmi.read_mtz_file(str(work_mtz));work.ensure_asu()
    dm=gemmi.read_mtz_file(str(dm_mtz));dm.ensure_asu()
    wl=work.column_labels();dl=dm.column_labels()
    fp_label,sig_label=_amp_labels(wl)
    quartet=['HL'+k for k in 'ABCD']
    if not all(c in dl for c in quartet):
        for prefix in ('DM','PIRATE','PARROT'):
            if all(prefix+'.ABCD.'+k in dl for k in 'ABCD'):
                quartet=[prefix+'.ABCD.'+k for k in 'ABCD'];break
        else:raise ValueError('No HL quartet in the density-modified MTZ')
    lookup={tuple(int(v) for v in row[:3]):row for row in np.array(dm)}
    rows=[]
    for row in np.array(work):
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        values=[match[dl.index(c)] for c in quartet]
        if not np.isfinite(values).all():continue
        kappa=float(np.hypot(values[0],values[1]))
        if kappa<=0:fom=0.0
        elif kappa>30:fom=1.0
        else:
            from scipy.special import i0e,i1e
            fom=float(i1e(kappa)/max(i0e(kappa),1e-300))
        phase=float(np.degrees(np.arctan2(values[1],values[0]))%360.)
        rows.append([*row[:3],row[wl.index(fp_label)],row[wl.index(sig_label)],
                     row[wl.index('FreeR_flag')],phase,fom])
    m=gemmi.Mtz(with_base=True);m.spacegroup=work.spacegroup;m.set_cell_for_all(work.cell);m.add_dataset('trace')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I'),('PHIB','P'),('FOM','W')]:m.add_column(label,kind)
    m.set_data(np.asarray(rows,np.float32));m.write_to_file(str(destination))
    return {'rows':len(rows),'hl_columns':quartet}


def _amp_labels(labels):
    """Accept both the staged F/SIGF schema and the FP/SIGFP schema."""
    for pair in (('FP','SIGFP'),('F','SIGF')):
        if pair[0] in labels and pair[1] in labels:return pair
    raise ValueError('No amplitude/sigma pair (FP/SIGFP or F/SIGF) in the input MTZ')


def merge_work(work_mtz, destination):
    """Restrict the observed amplitudes to work reflections, keeping FP/SIGFP/FreeR_flag."""
    m=gemmi.read_mtz_file(str(work_mtz));m.ensure_asu();labels=m.column_labels()
    data=np.array(m)
    fp_label,sig_label=_amp_labels(labels)
    keep=data[:,labels.index('FreeR_flag')]!=0
    out=gemmi.Mtz(with_base=True);out.spacegroup=m.spacegroup;out.set_cell_for_all(m.cell);out.add_dataset('work')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I')]:out.add_column(label,kind)
    rows=np.column_stack([data[keep][:,:3],data[keep][:,labels.index(fp_label)],
        data[keep][:,labels.index(sig_label)],data[keep][:,labels.index('FreeR_flag')]])
    out.set_data(rows.astype(np.float32));out.write_to_file(str(destination))
    return {'work_rows':int(keep.sum()),'excluded':int((~keep).sum())}


def work_from_experimental(experimental,sad,free_label,free_value,destination):
    """Write FP/SIGFP/FreeR_flag from the phase-stage map plus the SAD free set.

    A phases-only run has no ``01_inputs/native.mtz``; Parrot still needs the
    observed amplitudes and the test set, both of which the SAD file carries.
    """
    exp=gemmi.read_mtz_file(str(experimental));exp.ensure_asu()
    obs=gemmi.read_mtz_file(str(sad));obs.ensure_asu()
    labels=obs.column_labels();data=np.array(obs)
    exp_labels=exp.column_labels();exp_data=np.array(exp)
    fp,sig=_amp_labels(exp_labels)
    flags={}
    for row in data:
        flags[(int(row[0]),int(row[1]),int(row[2]))]=row[labels.index(free_label)]
    rows=[]
    for row in exp_data:
        key=(int(row[0]),int(row[1]),int(row[2]))
        rows.append([row[0],row[1],row[2],row[exp_labels.index(fp)],row[exp_labels.index(sig)],
                     flags.get(key,free_value)])
    out=gemmi.Mtz(with_base=True);out.spacegroup=exp.spacegroup;out.set_cell_for_all(exp.cell)
    out.add_dataset('work')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I')]:out.add_column(label,kind)
    out.set_data(np.array(rows,dtype=np.float32));out.write_to_file(str(destination))
    return str(destination)
