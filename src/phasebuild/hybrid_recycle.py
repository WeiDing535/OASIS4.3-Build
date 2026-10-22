"""Explicit PhaseBuild-seeded external Buccaneer/REFMAC recycling.

This is a composite method, never reported as an autonomous/native tracer.
No target-reference argument exists. The fixed endpoint is the final round.
"""
import json,os,shutil,subprocess,time,signal,threading
import xml.etree.ElementTree as ET
from pathlib import Path
import gemmi
from .data import fingerprint,audit
from .refinement import refmac_refine
from .runtime import resolve_runtime, runtime_environment, runtime_info

def reuse_first_build(source, destination, expected_command, report):
 """Reuse only a proven completed first construction, never a refined model."""
 source=Path(source).resolve();destination=Path(destination)
 old=json.loads((source.parent/'hybrid_report.json').read_text())
 if old['parameters']!=report['parameters'] or old['inputs']!=report['inputs']:
  raise ValueError('Cached construction inputs/protocol differ')
 if old['binary']['sha256']!=report['binary']['sha256'] or [r['sha256'] for r in old['training_library']]!=[r['sha256'] for r in report['training_library']]:
  raise ValueError('Cached construction binary/training resources differ')
 actual=json.loads((source/'command.json').read_text())
 def normalized(command):
  command=list(command)
  for flag in ['-pdbout','-xmlout']:
   if command.count(flag)!=1:raise ValueError('Invalid cached construction command')
   command[command.index(flag)+1]='<output>'
  return command
 if normalized(actual)!=normalized(expected_command):raise ValueError('Cached construction scientific command differs')
 cycles=ET.parse(str(source/'buccaneer.xml')).getroot().findall('./Cycles/Cycle')
 if [int(c.findtext('Number')) for c in cycles]!=[1,2,3] or 'cbuccaneer: Normal termination' not in (source/'buccaneer.log').read_text():
  raise ValueError('Cached first construction did not complete three cycles')
 files=[{'relative_path':str(p.relative_to(source)),'original':fingerprint(p)} for p in sorted(source.rglob('*')) if p.is_file()]
 for name in ['model.cif','model.pdb']:
  if not (source/name).is_file():raise ValueError('Cached construction model missing: '+name)
 shutil.copytree(source,destination)
 for item in files:
  item['copied']=fingerprint(destination/item['relative_path'])
  if item['copied']['sha256']!=item['original']['sha256']:raise RuntimeError('Cached construction changed during copy')
 return {'source':str(source),'source_report':fingerprint(source.parent/'hybrid_report.json'),
         'original_model':fingerprint(source/'model.cif'),'copied_model':fingerprint(destination/'model.cif'),
         'files':files,'construction_executed_again':False,
         'execution':'Completed first3cycles copied byte-for-byte; resume at failed first refinement'}

def recycle(mtz,native,model,sequence,output,rounds=5,refine_cycles=10,ccp4=None,build_semet=False,progress=print,existing_first_build=None):
 if not isinstance(rounds,int) or not 1<=rounds<=20:raise ValueError('Rounds must be 1..20')
 if not isinstance(refine_cycles,int) or not 1<=refine_cycles<=200:raise ValueError('REFMAC cycles must be 1..200')
 info=audit(mtz);required={'F':'F','SIGF':'Q','FWT':'F','PHWT':'P','FOM':'W'}
 for label,kind in required.items():
  found=[c for c in info['columns'] if c['label']==label]
  if len(found)!=1 or found[0]['type']!=kind:raise ValueError('Missing/ambiguous required map column '+label)
 root=resolve_runtime(ccp4,required=('cbuccaneer','refmac5'));library=root/'lib/data/reference_structures/reference-1tqw';binary=root/'bin/cbuccaneer'
 for suffix in ['.mtz','.pdb']:
  if not library.with_suffix(suffix).is_file():raise ValueError('Runtime training resource unavailable: '+str(library.with_suffix(suffix)))
 out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);previous=Path(model).resolve() if model is not None else None;current=Path(mtz).resolve();started=time.time()
 report=dict(engine='hybrid-phasebuild-seeded-external-buccaneer-refmac',status='running',pid=os.getpid(),parameters=dict(rounds=rounds,refine_cycles=refine_cycles,build_semet=build_semet,model_filter=True),inputs={k:fingerprint(v) for k,v in [('initial_map',mtz),('native',native),('phasebuild_seed_model',model),('sequence',sequence)] if v is not None},binary=fingerprint(binary),training_library=[fingerprint(library.with_suffix(s)) for s in ['.mtz','.pdb']],stages=[],scope='PhaseBuild supplies starting coordinates; external Buccaneer performs all subsequent construction. This is a composite method, not native PhaseBuild. Fixed final round, no Rfree-based model selection.')
 if model is None:
  report.update(engine='phasebuild-integrated-external-buccaneer-refmac',
   starting_model_source='None; first Buccaneer construction starts from experimental phases',
   scope='OASIS experimental phases initialize external Buccaneer construction without a model seed. PhaseBuild orchestrates external Buccaneer/REFMAC. Fixed final round; no Rfree-based model selection.')
 report['runtime']=runtime_info(root)
 sources=[fingerprint(f) for f in sorted((Path(__file__).parent).glob('*.py'))]
 def save():
  temp=out/'hybrid_report.tmp';temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(out/'hybrid_report.json')
 save();env=runtime_environment(root,scratch=out/'scratch');env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 proc=None;old_term=None
 if threading.current_thread() is threading.main_thread():
  def stop(signum,frame):raise InterruptedError('Hybrid construction cancelled')
  old_term=signal.signal(signal.SIGTERM,stop)
 try:
  for number in range(1,rounds+1):
   folder=out/('round_%02d_build'%number)
   # CIF input retains labels for fragments beyond the PDB single-chain alphabet.
   command=[str(binary),'-mtzin',str(current),'-mtzin-ref',str(library.with_suffix('.mtz')),'-pdbin-ref',str(library.with_suffix('.pdb')),'-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]','-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]','-colin-fo','/*/*/[F,SIGF]','-colin-phifom','/*/*/[PHWT,FOM]','-colin-fc','/*/*/[FWT,PHWT]','-pdbin',str(previous),'-model-filter','-seqin',str(Path(sequence).resolve()),'-pdbout',str(folder/'model.pdb'),'-xmlout',str(folder/'buccaneer.xml')]
   if previous is None:
    index=command.index('-pdbin');del command[index:index+2]
   command+=['-'+s for s in ['find','grow','join','link','sequence','correct','filter','ncsbuild','prune','rebuild']]+['-cycles','3' if number==1 else '2','-jobs','1','-no-correlation-mode' if number==1 else '-correlation-mode']
   if build_semet:command+=['-build-semet']
   if number==1 and existing_first_build is not None:
    report['resumed_first_construction']=reuse_first_build(existing_first_build,folder,command,report)
    report.update(current='reused first construction',current_round=number);save()
    progress('Reused completed first3constructioncycles; resume at original-HL refinement1')
   else:
    folder.mkdir();progress('Hybrid construction round %d'%number)
    (folder/'command.json').write_text(json.dumps(command,indent=2)+'\n');report.update(current='buccaneer',current_round=number);save()
    with (folder/'buccaneer.log').open('w') as log:
     proc=subprocess.Popen(command,cwd=folder,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);report['child_pid']=proc.pid;save()
     try:code=proc.wait(timeout=21600)
     except subprocess.TimeoutExpired:
      os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=15);raise
    if code:raise RuntimeError('External Buccaneer failed in round '+str(number))
    structure=gemmi.read_structure(str(folder/'model.pdb'));structure.make_mmcif_document().write_file(str(folder/'model.cif'));report.pop('child_pid',None)
   report['current']='refmac';save();progress('Hybrid REFMAC round %d'%number)
   refinement=out/('round_%02d_refine'%number);ref=refmac_refine(native,folder/'model.cif',refinement,refine_cycles,ccp4,timeout=3600)
   report['stages'].append(dict(round=number,construction='external-Buccaneer',buccaneer_cycles=3 if number==1 else 2,refmac_cycles=refine_cycles,statistics=ref['statistics']['final'],free_flags_preserved=ref['free_flags_preserved']));save();current=refinement/'refined.mtz';previous=refinement/'model.cif'
  shutil.copyfile(previous,out/'primary_model.cif');report.update(status='completed',primary_model=fingerprint(out/'primary_model.cif'),primary_metric=report['stages'][-1]['statistics'])
 except BaseException as error:report.update(status='failed',error=str(error));raise
 finally:
  # GUI cancellation targets the CLI process group. Buccaneer has its own
  # session for timeout control, so forward cancellation to that child session.
  if proc is not None and proc.poll() is None:
   try:
    os.killpg(proc.pid,signal.SIGTERM)
    # Finish child cleanup before the GUI's 3-second hard-stop deadline.
    try:proc.wait(timeout=1)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait(timeout=1)
   except ProcessLookupError:pass
  if old_term is not None:signal.signal(signal.SIGTERM,old_term)
  report.update(seconds=time.time()-started,sources=sources,source_changed_during_run=any(fingerprint(f['path'])!=f for f in sources));save()
 return report
