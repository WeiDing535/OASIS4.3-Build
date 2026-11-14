"""Structure-factor calculation for the feedback loop (no refinement).

The density-modification feedback needs Fcalc/PHIC for a work-only partial
model without moving a single atom.  ``attach_fc`` writes those two columns onto
the staged MTZ so the next handdm cycle can build the combined likelihood:
a bundled REFMAC pass with zero refinement cycles (the default), or Phenix
``fmodel`` when an explicit binary directory is supplied.
"""
from pathlib import Path
import os
import subprocess
import gemmi
import numpy as np


def attach_fc(staged, model_pdb, destination, phenix_bin=None, runtime=None, fcalc='fc_all'):
    """Attach bulk-solvent-corrected FC/PHIC to the staged MTZ.

    Both paths compute structure factors without moving any atom, so the traced
    partial model enters the next density-modification cycle unchanged. Phenix
    ``fmodel`` is used only when an explicit binary directory is supplied;
    otherwise the bundled REFMAC runs with zero refinement cycles, which yields
    the same FC/PHIC (and the same FWT/PHWT) columns.
    """
    if phenix_bin:
        return _attach_fc_phenix(staged,model_pdb,destination,phenix_bin)
    return _attach_fc_refmac0(staged,model_pdb,destination,runtime,fcalc)


def _refmac_fc_keywords(cycles=0):
    return ("title fcalc\n"
            "ncyc %d\n"%int(cycles) +
            "labin -\n  FP=FP SIGFP=SIGFP -\n  FREE=FreeR_flag\n"
            "weight AUTO\n"
            "make check NONE\n"
            "make hydrogen NO hout NO peptide NO cispeptide YES ssbridge YES symmetry YES sugar YES "
            "connectivity NO link NO\n"
            "refi type REST PHASE resi MLKF meth CGMAT bref ISOT\n"
            "scal type SIMP LSSC ANISO EXPE\n"
            "solvent YES VDWProb 1.4 IONProb 0.8 RSHRink 0.8\n"
            "PHOUT\nPNAME buccaneer\nDNAME buccaneer\nEND\n")


def _attach_fc_refmac0(staged, model_pdb, destination, runtime=None, fcalc='fc_all'):
    if fcalc not in ('fc','fc_all','refined'):raise ValueError('fcalc must be fc, fc_all or refined')
    """FC/PHIC from a zero-cycle REFMAC pass with the bundled runtime."""
    if runtime is None:raise ValueError('A runtime directory is required for the bundled REFMAC Fcalc pass')
    executable=Path(runtime)/'bin'/'refmac5'
    if not executable.is_file():raise ValueError('Bundled refmac5 is missing')
    source=gemmi.read_mtz_file(str(staged));labels=source.column_labels()
    fp_label='FP' if 'FP' in labels else 'F'
    sig_label='SIGFP' if 'SIGFP' in labels else 'SIGF'
    base=np.array(source,copy=True)
    destination=Path(destination).resolve()
    work=Path(destination).parent/('refmac0_'+Path(destination).stem)
    if work.exists():raise ValueError('REFMAC Fcalc workspace already exists: '+str(work))
    work.mkdir(parents=True)
    native=(work/'native.mtz').resolve()
    m=gemmi.read_mtz_file(str(staged));copy=gemmi.Mtz(with_base=True)
    copy.spacegroup=m.spacegroup;copy.set_cell_for_all(m.cell);copy.add_dataset('native')
    for label,kind in [('FP','F'),('SIGFP','Q'),('FreeR_flag','I')]:copy.add_column(label,kind)
    data=np.column_stack([base[:,:3],base[:,labels.index(fp_label)],base[:,labels.index(sig_label)],
                          (np.arange(base.shape[0])%20!=0).astype(float)])
    copy.set_data(data.astype(np.float32));copy.write_to_file(str(native))
    fc_mtz=(work/'fc.mtz').resolve();fc_pdb=(work/'fc.pdb').resolve()
    env={'CCP4':str(Path(runtime).resolve())}
    environment=dict(os.environ)
    environment.update({'PATH':str(executable.parent)+':'+environment.get('PATH','/usr/bin:/bin'),
        'CCP4':str(Path(runtime).resolve()),'CCP4_SCR':str(work/'scratch'),
        'CLIBD':str(Path(runtime)/'lib/data'),'CLIBD_MON':str(Path(runtime)/'lib/data/monomers')+'/',
        'CINCL':str(Path(runtime)/'include'),
        'SYMINFO':str(Path(runtime)/'lib/data/syminfo.lib')})
    (work/'scratch').mkdir(exist_ok=True)
    with (work/'refmac.log').open('w') as log:
        completed=subprocess.run([str(executable),'HKLIN',str(native),'XYZIN',str(Path(model_pdb).resolve()),
            'HKLOUT',str(fc_mtz),'XYZOUT',str(fc_pdb)],
            input=_refmac_fc_keywords(1 if fcalc=='refined' else 0),text=True,
            stdout=log,stderr=subprocess.STDOUT,env=environment,cwd=str(work))
    if completed.returncode or not fc_mtz.is_file():
        raise RuntimeError('REFMAC Fcalc pass failed for '+str(model_pdb))
    model_mtz=gemmi.read_mtz_file(str(fc_mtz));ml=model_mtz.column_labels();ma=np.array(model_mtz)
    lookup={tuple(int(v) for v in row[:3]):row for row in ma}
    fcs=np.zeros(source.nreflections);phis=np.zeros(source.nreflections)
    for i,row in enumerate(base):
        match=lookup.get(tuple(int(v) for v in row[:3]))
        if match is None:continue
        columns=('FC_ALL','PHIC_ALL') if fcalc=='fc_all' else ('FC','PHIC')
        fcs[i]=match[ml.index(columns[0])];phis[i]=match[ml.index(columns[1])]
    source.add_column('FC','F');source.add_column('PHIC','P')
    source.set_data(np.column_stack([base,fcs,phis]).astype(np.float32))
    source.write_to_file(str(destination))
    return {'engine':'bundled-refmac5-zero-cycle','binary':str(executable)}


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
