"""Explicit external Buccaneer postprocessing of a PhaseBuild-built model.

The adapter never searches/grows new backbone. Its output is a hybrid result,
not a native PhaseBuild result. Local CCP4 supplies the executable/library.
"""
import json,os,subprocess,time
from pathlib import Path
import gemmi
from .data import fingerprint
from .runtime import resolve_runtime, runtime_environment
STAGES=('join','link','sequence','correct','filter','ncsbuild','prune','rebuild')

def postprocess(mtz,model,sequence,output,stages=('sequence','rebuild'),cycles=1,ccp4=None,build_semet=False,timeout=3600):
 if not stages or any(s not in STAGES for s in stages):raise ValueError('Only explicit backbone postprocessing stages are allowed')
 if not isinstance(cycles,int) or not 1<=cycles<=20:raise ValueError('Cycles must be 1..20')
 root=resolve_runtime(ccp4,required=('cbuccaneer',));library=root/'lib/data/reference_structures/reference-1tqw';out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
 binary=root/'bin/cbuccaneer';st=gemmi.read_structure(str(model));st.write_pdb(str(out/'input.pdb'))
 command=[str(binary),'-mtzin',str(Path(mtz).resolve()),'-mtzin-ref',str(library.with_suffix('.mtz')),'-pdbin-ref',str(library.with_suffix('.pdb')),'-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]','-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]','-colin-fo','/*/*/[F,SIGF]','-colin-phifom','/*/*/[PHWT,FOM]','-colin-fc','/*/*/[FWT,PHWT]','-seqin',str(Path(sequence).resolve()),'-pdbin',str(out/'input.pdb'),'-pdbout',str(out/'model.pdb'),'-xmlout',str(out/'buccaneer.xml')]
 command+=['-'+s for s in STAGES if s in stages]+['-cycles',str(cycles),'-jobs','1','-no-correlation-mode']
 if build_semet:command+=['-build-semet']
 report=dict(engine='hybrid-phasebuild-with-external-buccaneer-postprocessing',command=command,input=fingerprint(mtz),model=fingerprint(model),sequence=fingerprint(sequence),binary=fingerprint(binary),training_library=[fingerprint(library.with_suffix(s)) for s in ['.mtz','.pdb']],stages=list(stages),cycles=cycles,scope='No find/grow stages requested. External join/link/sequence/correction/prune/rebuild contribution is explicit; not a native-engine result.')
 start=time.time()
 try:
  env=runtime_environment(root);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  with (out/'buccaneer.log').open('w') as log:result=subprocess.run(command,cwd=out,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
  if result.returncode:raise RuntimeError('Buccaneer postprocessing failed; inspect log')
  text=(out/'buccaneer.log').read_text()
  # Stage counters in the log provide an independent check against accidental full building.
  if 'C-alphas after finding:' in text or 'C-alphas after growing:' in text:raise RuntimeError('Unexpected find/grow stage in postprocessing')
  st=gemmi.read_structure(str(out/'model.pdb'));st.make_mmcif_document().write_file(str(out/'model.cif'));report.update(status='success',residues=sum(1 for c in st[0] for r in c if r.find_atom('CA','*')),output_model=fingerprint(out/'model.cif'))
 except BaseException as error:report.update(status='failed',error=str(error));raise
 finally:report['seconds']=time.time()-start;(out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
 return report
