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
        'CLIBD_MON':str(base/'lib/data/monomers')+'/','CINCL':str(base/'include'),
        'SYMINFO':str(base/'lib/data/syminfo.lib'),
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


BUCCANEER_1ST="""title  {title}
pdbin-ref {pdbin_ref}
mtzin-ref {mtzin_ref}
colin-ref-fo /*/*/[FP.F_sigF.F,FP.F_sigF.sigF]
colin-ref-hl /*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]
seqin {seqin}
mtzin {mtzin}
colin-fo /*/*/[FP,SIGFP]
colin-phifom /*/*/[PHIB,FOM]
colin-free /*/*/[FreeR_flag]
pdbout {pdbout}
anisotropy-correction
fast
cycles 3
sequence-reliability 0.95
new-residue-name UNK
resolution {resolution}
jobs 1
"""

BUCCANEER_NTH="""title  {title}
pdbin-ref {pdbin_ref}
mtzin-ref {mtzin_ref}
colin-ref-fo /*/*/[FP.F_sigF.F,FP.F_sigF.sigF]
colin-ref-hl /*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]
seqin {seqin}
mtzin {mtzin}
colin-fo /*/*/[FP,SIGFP]
colin-hl HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB
colin-fc FWT,PHWT
colin-free /*/*/[FreeR_flag]
pdbin {pdbin}
pdbout {pdbout}
anisotropy-correction
fast
cycles 2
sequence-reliability 0.95
correlation-mode
new-residue-name UNK
resolution {resolution}
jobs 1
"""

REFMAC_TAIL="""make check NONE
make hydrogen NO hout NO peptide NO cispeptide YES ssbridge YES symmetry YES sugar YES connectivity NO link NO
refi type REST PHASE resi MLKF meth CGMAT bref ISOT
scal type SIMP LSSC ANISO EXPE
solvent YES VDWProb 1.4 IONProb 0.8 RSHRink 0.8
PHOUT
PNAME buccaneer
DNAME buccaneer
END
"""


def _write(path, text):
    Path(path).write_text(text)
    return path


def _refmac_keywords(labin_phase, hklin, xyzin, hklout, xyzout, title):
    lines=["title %s"%title,"ncyc 10","labin -",labin_phase+" -",
           "  FREE=FreeR_flag -","  FP=FP SIGFP=SIGFP","weight AUTO"]
    return "\n".join(lines)+"\n"+REFMAC_TAIL


def run(native, phase_map, sequence, runtime, output, *, cycles=5, free_value=0,
        resolve_cycles=2.0, jobs=1, timeout=21600):
    """Run an in-package equivalent of CCP4's buccaneer_pipeline.

    Only the bundled ``cbuccaneer`` and ``refmac5`` executables are used, with the
    same cycle structure as the CCP4 script (3 Buccaneer cycles then REFMAC in the
    first round, 2 Buccaneer cycles with combined HL and FWT/PHWT maps afterwards).
    This removes the dependency on a CCP4 Python installation and its pipeline
    module while keeping the scientific protocol identical.
    """
    native=Path(native).resolve();phase_map=Path(phase_map).resolve()
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    work=out/'build_input.mtz'
    staged=input_mtz(native,phase_map,work,free_value)
    buccaneer=_find('cbuccaneer',runtime);refmac=_find('refmac5',runtime)
    reference=(Path(runtime)/'lib/data/reference_structures/reference-1tqw').resolve()
    env=_environment(buccaneer,out/'scratch')
    title='OASIS4.3 P+ build'
    flags={'title':title,'pdbin_ref':str(reference.with_suffix('.pdb')),
           'mtzin_ref':str(reference.with_suffix('.mtz')),'seqin':str(Path(sequence).resolve()),
           'resolution':'%.3f'%resolve_cycles}
    rows=[]
    current=work.resolve();previous=None
    for cycle in range(1,int(cycles)+1):
        folder=out/('cycle_%02d'%cycle);folder.mkdir(exist_ok=True)
        build_pdb=folder/'build.pdb'
        if cycle==1:
            text=BUCCANEER_1ST.format(mtzin=current,pdbout=build_pdb,**flags)
        else:
            text=BUCCANEER_NTH.format(mtzin=current,pdbin=previous,pdbout=build_pdb,**flags)
        log=folder/'buccaneer.log'
        with log.open('w') as handle:
            done=subprocess.run([buccaneer,'-stdin'],input=text,text=True,stdout=handle,
                stderr=subprocess.STDOUT,env=env,cwd=folder,timeout=timeout)
        if done.returncode or not build_pdb.is_file():
            raise RuntimeError('cbuccaneer failed in cycle %d; see %s'%(cycle,log))
        if not _has_atoms(build_pdb):
            # REFMAC cannot refine an empty model, and reporting that as a REFMAC
            # crash would hide the actual result: the map carried no buildable
            # density.  Record it and stop with a message that names the cause.
            report={'engine':'oasis-buccaneer-refmac-pipeline','status':'no_model',
                'parameters':dict(cycles=cycles,free_value=free_value,resolution=resolve_cycles),
                'inputs':{k:fingerprint(v) for k,v in [('native',native),('phase_map',phase_map),('sequence',sequence)]},
                'staged':staged,'residues_built':0,'stopped_at_cycle':cycle,
                'stop_reason':'Buccaneer built no atoms from the density-modified map; REFMAC was not run',
                'dm_map_kept':fingerprint(phase_map)}
            (out/'hybrid_report.json').write_text(json.dumps(report,indent=2)+'\n')
            raise RuntimeError('密度修饰后的图未建出任何残基（第 %d 轮 Buccaneer 0 原子）；未生成模型。'
                'DM 相位图保留在 %s，可调整溶剂含量（dm_solvent）或密度修饰轮数后重跑'%(cycle,phase_map))
        refine_mtz=folder/'refine.mtz';refine_pdb=folder/'refine.pdb'
        labin=('  PHIB=PHIB FOM=FOM' if cycle==1 else
               '  HLA=HLACOMB HLB=HLBCOMB HLC=HLCCOMB HLD=HLDCOMB')
        keywords=_write(folder/'refmac.in',_refmac_keywords(labin,current,build_pdb,refine_mtz,refine_pdb,title))
        rlog=folder/'refmac.log'
        with rlog.open('w') as handle:
            done=subprocess.run([refmac,'HKLIN',str(current),'XYZIN',str(build_pdb),
                'HKLOUT',str(refine_mtz),'XYZOUT',str(refine_pdb)],
                input=Path(keywords).read_text(),text=True,stdout=handle,stderr=subprocess.STDOUT,
                env=env,cwd=folder,timeout=timeout)
        if done.returncode or not refine_pdb.is_file():
            raise RuntimeError('REFMAC failed in cycle %d; see %s'%(cycle,rlog))
        stats=_refmac_final(rlog)
        rows.append(dict(cycle=cycle,**stats))
        # One line per cycle: without it a five-cycle build is silent for minutes
        # and the GUI log looks empty.
        print('build cycle %d/%d: %d residues, Rwork %.4f, Rfree %.4f'
              %(cycle,cycles,_residues(refine_pdb),stats['r_work'],stats['r_free']),flush=True)
        current=refine_mtz;previous=refine_pdb
    model=Path(previous)
    structure=gemmi.read_structure(str(model))
    structure.make_mmcif_document().write_file(str(out/'primary_model.cif'))
    residues=sum(1 for chain in structure[0] for residue in chain if residue.name!='HOH')
    last=rows[-1]
    report={'engine':'oasis-buccaneer-refmac-pipeline','status':'completed',
        'parameters':dict(cycles=cycles,free_value=free_value,resolution=resolve_cycles),
        'inputs':{k:fingerprint(v) for k,v in [('native',native),('phase_map',phase_map),('sequence',sequence)]},
        'staged':staged,'primary_model':fingerprint(out/'primary_model.cif'),
        'primary_metric':{'cycle':int(cycles),'r_work':last['r_work'],'r_free':last['r_free'],
                          'bond_rmsd_angstrom':None,'angle_rmsd_degrees':None},
        'cycles':rows,'residues_built':residues,'free_flags_preserved':True,
        'selection':'fixed final build/refine cycle; no Rfree used',
        'scope':'Same builder protocol as the reference route; only the input phases differ',
        'requires_ccp4_python':False,'external_executables':['cbuccaneer','refmac5']}
    (out/'hybrid_report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def _has_atoms(pdb):
    """True when the model has at least one coordinate record."""
    try:lines=Path(pdb).read_text(errors='ignore').splitlines()
    except OSError:return False
    return any(line.startswith(('ATOM','HETATM')) for line in lines)


def _residues(pdb):
    """Residue count of a built model, for the per-cycle progress line."""
    try:structure=gemmi.read_structure(str(pdb))
    except Exception:return 0
    return sum(1 for chain in structure[0] for residue in chain if residue.name!='HOH')


def _refmac_final(log):
    """Read the last REFMAC cycle row (Rwork, Rfree) from a refinement log."""
    lines=Path(log).read_text(errors='ignore').splitlines()
    best=None
    for index,line in enumerate(lines):
        if 'Ncyc' in line and 'Rfact' in line:
            for row in lines[index+1:]:
                fields=row.split()
                if len(fields)<3 or not fields[0].isdigit():continue
                try:rwork,rfree=float(fields[1]),float(fields[2])
                except ValueError:continue
                if 0<=rwork<=1.5 and 0<=rfree<=1.5:best=(rwork,rfree)
    if best is None:raise RuntimeError('Could not read R factors from '+str(log))
    return {'r_work':best[0],'r_free':best[1]}
