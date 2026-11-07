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
          cycles=5, jobs=1, timeout=3600):
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
    lookup={tuple(int(v) for v in row[:3]):row for row in np.array(dm)}
    quartet=['HL'+k for k in 'ABCD']
    rows=[]
    for row in np.array(work):
        if row[work_labels.index('FreeR_flag')]==free_value:continue
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        f=row[work_labels.index('FP')];sig=row[work_labels.index('SIGFP')]
        values=[match[dm_labels.index(c)] for c in quartet]
        if not np.isfinite(values).all():continue
        kappa=float(np.hypot(values[0],values[1]));fom=min(kappa,1.0)
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
        executable=None, trace_cycles=5, progress=print):
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
            runtime=runtime,scratch=folder/'scratch')
        row={'cycle':cycle+1,'dm':str(dm_path),'dm_info':info}
        if cycle<feedback_cycles:
            traced=trace(merged,dm_path,sequence,runtime,folder/'trace',
                         resolution=trace_resolution,cycles=trace_cycles,progress=progress)
            row['trace']=traced;model=traced['model']
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


def merge_work(work_mtz, destination):
    """Restrict the observed amplitudes to work reflections, keeping FP/SIGFP/FreeR_flag."""
    m=gemmi.read_mtz_file(str(work_mtz));m.ensure_asu();labels=m.column_labels()
    data=np.array(m)
    keep=data[:,labels.index('FreeR_flag')]!=0
    out=gemmi.Mtz(with_base=True);out.spacegroup=m.spacegroup;out.set_cell_for_all(m.cell);out.add_dataset('work')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I')]:out.add_column(label,kind)
    rows=np.column_stack([data[keep][:,:3],data[keep][:,labels.index('FP')],
        data[keep][:,labels.index('SIGFP')],data[keep][:,labels.index('FreeR_flag')]])
    out.set_data(rows.astype(np.float32));out.write_to_file(str(destination))
    return {'work_rows':int(keep.sum()),'excluded':int((~keep).sum())}
