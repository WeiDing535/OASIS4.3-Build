"""Established full Buccaneer control, with explicit inputs and provenance."""
import json
import os
import subprocess
import time
from pathlib import Path
import gemmi
from .data import audit,fingerprint
from .sequencing import read_sequences
from .runtime import resolve_runtime, runtime_environment


def buccaneer_full(mtz,sequence,output,cycles=3,ccp4=None,timeout=600,correlation=False,
                   starting_model=None):
    if cycles<1:
        raise ValueError('Positive Buccaneer cycle count required')
    records=read_sequences(sequence)
    info=audit(mtz)
    labels={c['label'] for c in info['columns']}
    if not {'F','SIGF','FWT','PHWT','FOM','FreeR_flag'}<=labels:
        raise ValueError('Control requires matching observations, map coefficients and free flags in one MTZ')
    root=resolve_runtime(ccp4,required=('cbuccaneer',));binary=root/'bin/cbuccaneer'
    reference=root/'lib/data/reference_structures/reference-1tqw'
    required=[binary,reference.with_suffix('.mtz'),reference.with_suffix('.pdb')]
    if not all(p.is_file() for p in required):
        raise ValueError('CCP4 executable or independent training reference unavailable')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    # Normalise only the file format, without inferring biological chain breaks.
    (out/'sequence.fasta').write_text(''.join('>%s\n%s\n'%(name,seq) for name,seq in records))
    command=[str(binary),'-mtzin',str(Path(mtz).resolve()),
             '-mtzin-ref',str(reference.with_suffix('.mtz')),'-pdbin-ref',str(reference.with_suffix('.pdb')),
             '-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
             '-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
             '-colin-fo','/*/*/[F,SIGF]','-colin-phifom','/*/*/[PHWT,FOM]',
             '-colin-fc','/*/*/[FWT,PHWT]','-colin-free','/*/*/[FreeR_flag]',
             '-seqin',str(out/'sequence.fasta'),'-pdbout',str(out/'model.pdb'),
             '-xmlout',str(out/'buccaneer.xml')]
    steps=['find','grow','join','link','sequence','correct','filter','ncsbuild','prune','rebuild']
    command+=['-'+step for step in steps]+['-cycles',str(cycles),'-jobs','1']
    if correlation:
        command.append('-correlation-mode')
    if starting_model is not None:
        command+=['-pdbin',str(Path(starting_model).resolve())]
    report={'engine':'external-Buccaneer-full','input':info['input'],'sequence':fingerprint(sequence),
            'training_reference':[fingerprint(p) for p in required[1:]],'binary':fingerprint(binary),
            'source':fingerprint(__file__),'command':command,'cycles':cycles,'steps':steps,
            'correlation_mode':bool(correlation),
            'starting_model':fingerprint(starting_model) if starting_model is not None else None,
            'scope':'Full sequence/sidechain building, no target reference. '
                    'Initial FWT/PHWT and sequence match the PhaseBuild inputs; internal weighting differs. '
                    'Uses phase-file F/SIGF on its map scale; subsequent common refinement uses original native F/SIGF. '
                    'This is a fixed single build plus refinement control, not the complete CCP4i2 iterative pipeline.'}
    environment=runtime_environment(root)
    started=time.perf_counter()
    try:
        with (out/'buccaneer.log').open('w') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,cwd=str(out),
                                  env=environment,timeout=timeout)
        report['returncode']=result.returncode
        if result.returncode or not (out/'model.pdb').is_file():
            raise RuntimeError('Buccaneer failed; see '+str(out/'buccaneer.log'))
        model=gemmi.read_structure(str(out/'model.pdb'))
        model.make_mmcif_document().write_file(str(out/'model.cif'))
        report.update(status='success',output=fingerprint(out/'model.cif'),
                      residues=sum(len(c) for c in model[0]),
                      assigned_residues=sum(r.name!='UNK' for c in model[0] for r in c),
                      heavy_atoms=sum(a.element.name not in ['H','D'] for c in model[0] for r in c for a in r))
    except Exception as error:
        report.update(status='failed',error=str(error))
        raise
    finally:
        report['seconds']=time.perf_counter()-started
        (out/'build_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
