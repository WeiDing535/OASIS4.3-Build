"""Experimental weak-phase PhaseBuild pipeline (not yet a validated default).

Inputs are map coefficients FWT/PHWT, sequence and optional native observations.
No target reference or Buccaneer model is accepted. With native data the primary
output is the last planned round, not the model selected by minimum Rfree.
"""
import argparse,json,math,os,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))

def plan(args):
 import gemmi
 from phasebuild.gui_workflow import defaults,make_plan
 from phasebuild.sequencing import read_sequences
 default_bank=ROOT/'runs/sad_four_optimization/models'
 bank=Path(args.profile_bank or default_bank).resolve()
 resolution=gemmi.read_mtz_file(args.map).resolution_high()
 if resolution>=2.35:backbone,residue=bank/'backbone_2p8.npz',bank/'residue_2p8.npz'
 else:backbone,residue=ROOT/'models/1tqw_1p9_symmetric.npz',ROOT/'models/residue_1tqw_1p9.npz'
 n=sum(len(seq) for _,seq in read_sequences(args.sequence));capacity=min(20000,int(math.ceil(n*1.25)))
 cfg=defaults(ROOT);cfg.update(mode='fresh',map=args.map,sequence=args.sequence,output=args.output,profile=str(backbone),residue_profile=str(residue),refine=bool(args.native),rounds=args.rounds if args.native else 1,cycles=args.cycles,map_mode='rebuild',max_residues=capacity)
 if args.native:cfg['native']=args.native
 cfg,steps=make_plan(cfg,ROOT,python=sys.executable)
 out=Path(cfg['output']);seeds=out/'00_seeds';trace=out/'01_trace';sequence=out/'02_sequence'
 first=[dict(name='Contextual tripeptide search',kind='search',output=str(seeds),command=[sys.executable,'-u',str(ROOT/'run_context_seeds.py'),cfg['map'],'--profile',str(backbone),'--count',str(args.seeds),'--output',str(seeds)]),dict(name='Bounded weak-density trace',kind='trace',output=str(trace),command=[sys.executable,'-u',str(ROOT/'run_weak_trace.py'),cfg['map'],'--profile',str(backbone),'--candidates',str(seeds/'candidates.npz'),'--capacity',str(capacity),'--output',str(trace)]),dict(name='Adaptive sequence placement',kind='build',output=str(sequence),command=[sys.executable,'-u',str(ROOT/'run_weak_sequence.py'),cfg['map'],'--model',str(trace/'backbone.cif'),'--sequence',cfg['sequence'],'--residue-profile',str(residue),'--output',str(sequence)])]
 return cfg,first+steps[2:]

def stage_command(step):
 """Choose profiles from the actual map, preserving stage argument boundaries."""
 import gemmi
 command=list(step['command']);script=next(x for x in command if x.endswith('.py'))
 index=command.index(script);name=Path(script).stem;map_path=command[index+1]
 if '--profile' in command and gemmi.read_mtz_file(map_path).resolution_high()<2.35:
  command[command.index('--profile')+1]=str(ROOT/'models/1tqw_1p9_symmetric.npz')
  if '--residue-profile' in command:
   command[command.index('--residue-profile')+1]=str(ROOT/'models/residue_1tqw_1p9.npz')
 if name in ['run_improve','run_complete','run_rebuild']:
  command[index:index+1]=[str(ROOT/'run_adaptive_stage.py'),name[4:]]
 return command

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--map',required=True);p.add_argument('--sequence',required=True);p.add_argument('--native');p.add_argument('--output',required=True);p.add_argument('--profile-bank');p.add_argument('--seeds',type=int,default=512);p.add_argument('--rounds',type=int,default=5);p.add_argument('--cycles',type=int,default=10);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
 if not 1<=a.seeds<=2048:p.error('--seeds must be 1..2048')
 cfg,steps=plan(a)
 if a.dry_run:print(json.dumps(dict(config=cfg,steps=steps),indent=2));return
 from phasebuild.data import fingerprint
 from phasebuild.gui_workflow import collect_results
 import gemmi
 out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=False);state=dict(config=cfg,steps=steps,status='running',pid=os.getpid(),completed=0,started=time.time(),experimental=True)
 sources=[fingerprint(f) for f in sorted((ROOT/'src/phasebuild').glob('*.py'))]+[fingerprint(f) for f in sorted(ROOT.glob('run_*.py'))]
 def save():
  temp=out/'gui_job.tmp';temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(out/'gui_job.json')
 save();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 try:
  with (out/'pipeline.log').open('w') as log:
   for step in steps:
    command=stage_command(step)
    state.update(current=step['name'],command=command);save();print(step['name'],flush=True);log.write(json.dumps(command)+'\n');log.flush()
    proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT);state['child_pid']=proc.pid;save();code=proc.wait()
    if code:raise RuntimeError('Stage failed: '+step['name'])
    state['completed']+=1;save()
  state.update(status='success',finished=time.time(),sources=sources,source_changed_during_run=any(fingerprint(f['path'])!=f for f in sources));state.pop('child_pid',None)
  results=collect_results(out);(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
  import shutil
  model=Path(steps[-1]['output'])/'model.cif';shutil.copyfile(model,out/'primary_model.cif')
 except BaseException as error:
  state.update(status='failed',error=str(error),traceback=traceback.format_exc());raise
 finally:save()
if __name__=='__main__':main()
