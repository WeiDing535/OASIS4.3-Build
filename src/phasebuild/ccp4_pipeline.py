"""CCP4 buccaneer_pipeline building recipe (same engine as the reference route).

The builder iterates Buccaneer construction with REFMAC refinement and its own
density modification. The free set is preserved verbatim and is never used to
choose between cycles; the fixed final cycle is reported, optionally guarded by
the existing work-R rule.
"""
from pathlib import Path
import json
import os
import shutil
import subprocess
import gemmi
import numpy as np

from .data import fingerprint


def _environment(executable, scratch):
    binary=Path(executable).resolve();base=binary.parent.parent
    env=dict(os.environ)
    env.update({'PATH':str(binary.parent)+':'+env.get('PATH','/usr/bin:/bin'),'CCP4':str(base),
        'CCP4_SCR':str(scratch),'CLIBD':str(base/'lib/data'),
        'CLIBD_MON':str(base/'lib/data/monomers')+'/','SYMINFO':str(base/'lib/data/syminfo.lib'),
        'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
        'VECLIB_MAXIMUM_THREADS':'1'})
    Path(scratch).mkdir(parents=True,exist_ok=True)
    return env


def _find(name, runtime):
    for candidate in (Path(runtime)/'bin'/name,):
        if candidate.is_file():return str(candidate)
    found=shutil.which(name)
    if found:return found
    for base in ('/Applications/ccp4-9','/opt/ccp4','/usr/local/ccp4'):
        candidate=Path(base)/'bin'/name
        if candidate.is_file():return str(candidate)
    raise ValueError(name+' was not found')


def input_mtz(native, phase_map, destination, free_value=0):
    """Merge native observations with the work-only DM phase/FOM columns.

    Every native row is kept so REFMAC can evaluate its own free set; phase and
    FOM are left undefined on free reflections so no held-out observation enters
    construction or density modification.
    """
    n=gemmi.read_mtz_file(str(native));n.ensure_asu();na=np.array(n);nl=n.column_labels()
    fp='FP' if 'FP' in nl else 'F'
    sig='SIGFP' if 'SIGFP' in nl else 'SIGF'
    d=gemmi.read_mtz_file(str(phase_map));d.ensure_asu();da=np.array(d);dl=d.column_labels()
    lookup={tuple(int(v) for v in row[:3]):row for row in da}
    rows=[]
    for row in na:
        key=tuple(int(v) for v in row[:3]);match=lookup.get(key)
        f=row[nl.index(fp)];s=row[nl.index(sig)];flag=row[nl.index('FreeR_flag')]
        if not np.isfinite([f,s,flag]).all() or min(f,s)<=0:continue
        if flag==free_value or match is None:
            rows.append([*key,f,s,flag,float('nan'),float('nan')]);continue
        ph=match[dl.index('PHIB')];fom=match[dl.index('FOM')]
        if not np.isfinite([ph,fom]).all():continue
        rows.append([*key,f,s,flag,ph,fom])
    m=gemmi.Mtz(with_base=True);m.spacegroup=n.spacegroup;m.set_cell_for_all(n.cell);m.add_dataset('build')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I'),('PHIB','P'),('FOM','W')]:
        m.add_column(label,kind)
    m.set_data(np.asarray(rows,np.float32));m.write_to_file(str(destination))
    return {'rows':len(rows),'free_value':free_value}


def run(native, phase_map, sequence, runtime, output, *, cycles=5, free_value=0,
        resolve_cycles=2.0, jobs=4, timeout=21600):
    """Run buccaneer_pipeline and expose the endpoint in the pipeline schema."""
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    work=out/'build_input.mtz'
    staged=input_mtz(native,phase_map,work,free_value)
    executable=_find('buccaneer_pipeline',runtime)
    reference=Path(runtime)/'lib/data/reference_structures/reference-1tqw'
    command=[executable,'-mtzin',str(work),'-seqin',str(sequence),
        '-colin-fo','/*/*/[FP,SIGFP]','-colin-phifom','/*/*/[PHIB,FOM]',
        '-colin-free','/*/*/[FreeR_flag]',
        '-mtzin-ref',str(reference.with_suffix('.mtz')),'-pdbin-ref',str(reference.with_suffix('.pdb')),
        '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
        '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
        '-cycles',str(int(cycles)),'-jobs',str(int(jobs)),'-buccaneer-fast',
        '-buccaneer-anisotropy-correction','-buccaneer-new-residue-name','UNK',
        '-buccaneer-resolution','%.3f'%resolve_cycles,'-refmac-mlhl','1',
        '-pdbout',str(out/'model.pdb')]
    (out/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    log_path=out/'pipeline.log'
    with log_path.open('w') as log:
        completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,env=_environment(executable,out/'scratch'),timeout=timeout)
    if completed.returncode:
        raise RuntimeError('buccaneer_pipeline failed; see '+str(log_path))
    model=out/'model.pdb'
    if not model.is_file():raise RuntimeError('buccaneer_pipeline produced no model')
    structure=gemmi.read_structure(str(model))
    structure.make_mmcif_document().write_file(str(out/'primary_model.cif'))
    text=log_path.read_text(errors='ignore')
    rwork=rfree=None;residues=0
    for line in text.splitlines():
        if 'Overall R factor' in line:rwork=float(line.split('=')[1])
        if 'Free R factor' in line:rfree=float(line.split('=')[1])
        if 'residues were built' in line:
            try:residues=int(line.split()[0])
            except ValueError:pass
    if rwork is None or rfree is None:
        raise RuntimeError('Could not read the final R factors from the builder log')
    report={'engine':'ccp4-buccaneer-pipeline','status':'completed','parameters':dict(cycles=cycles,
            free_value=free_value,resolution=resolve_cycles),
        'inputs':{k:fingerprint(v) for k,v in [('native',native),('phase_map',phase_map),('sequence',sequence)]},
        'staged':staged,'primary_model':fingerprint(out/'primary_model.cif'),
        'primary_metric':{'cycle':int(cycles),'r_work':float(rwork),'r_free':float(rfree),
                          'bond_rmsd_angstrom':None,'angle_rmsd_degrees':None},
        'residues_built':residues,
        'free_flags_preserved':True,'selection':'fixed final build/refine cycle; no Rfree used',
        'scope':'Same builder as the reference route; only the input phases differ'}
    (out/'hybrid_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
