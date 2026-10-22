"""Opt-in weak-map PhaseBuild search/growth using only map and unrelated profiles."""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
from phasebuild.data import Density,fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.trace import find_seeds,grow,grow_backward,write_backbone
from phasebuild.weak_trace import grow_weak,supported_intervals
from phasebuild.selection import independent_set

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mtz');p.add_argument('--profile',required=True);p.add_argument('--output',required=True);p.add_argument('--candidates');p.add_argument('--refine-count',type=int,default=512);p.add_argument('--max-length',type=int,default=30);p.add_argument('--capacity',type=int,required=True);p.add_argument('--seed-limit',type=int);p.add_argument('--legacy-growth',action='store_true');p.add_argument('--legacy-selection',action='store_true');args=p.parse_args()
 out=Path(args.output);out.mkdir(parents=True,exist_ok=False);started=time.time();d=Density(args.mtz);d.profile=DensityProfile.load(args.profile)
 sources=[fingerprint(p) for p in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
 if args.candidates:
  with np.load(args.candidates,allow_pickle=False) as z:
   assert json.loads(str(z['metadata']))['input']['sha256']==d.source['sha256']
   seeds=list(zip(z['seed_scores'],z['seeds'],z['seed_psi']))
 else:seeds,_=find_seeds(d,refine_count=args.refine_count,progress=lambda m:print(m,flush=True))
 if args.seed_limit:seeds=seeds[:args.seed_limit]
 # Re-rank geometries by the selected resolution-matched map target.
 seeds=sorted(seeds,key=lambda s:float(d.score(s[1]).mean()),reverse=True)
 grown=[]
 for i,(score,xyz,psi) in enumerate(seeds):
  if np.mean(d.score(xyz))<2.5 or np.min(d.score(xyz))<1.5:continue
  used=[ca for c in grown for ca in c[:,1]]
  if used and any(d.distance(ca,old)<2.5 for ca in xyz[:,1] for old in used):continue
  if args.legacy_growth:
   chain=grow(d,xyz,float(psi),args.max_length,4,2.);chain=grow_backward(d,chain,args.max_length,4,2.,used)
  else:
   chain=grow_weak(d,xyz,args.max_length,4,used=used);chain=grow_weak(d,chain,args.max_length,4,backward=True,used=used)
  grown.append(chain)
  print('Seed %d/%d: grew %d residues'%(i+1,len(seeds),len(chain)),flush=True)
 pieces=[]
 for chain in grown:
  scores=d.score(chain)
  if args.legacy_selection:
   strong=scores>=2.5;starts=np.flatnonzero(strong&~np.r_[False,strong[:-1]]);ends=np.flatnonzero(strong&~np.r_[strong[1:],False])+1;intervals=[(a,b) for a,b in zip(starts,ends) if b-a>=3]
  else:intervals=supported_intervals(scores)
  for a,b in intervals:pieces.append((float(np.sum(scores[a:b]-2.5)),chain[a:b]))
 pieces.sort(key=lambda x:x[0],reverse=True);pieces=pieces[:max(160,min(512,args.capacity))];conflicts=[0]*len(pieces)
 for i,(_,a) in enumerate(pieces):
  for j,(_,b) in enumerate(pieces[:i]):
   if any(d.distance(x,y)<2.5 for x in a[:,1] for y in b[:,1]):conflicts[i]|=1<<j;conflicts[j]|=1<<i
 selected=independent_set([x[0] for x in pieces],[len(x[1]) for x in pieces],conflicts,args.capacity)
 chains=[pieces[i][1] for i in selected];write_backbone(out/'backbone.cif',chains,d)
 np.savez_compressed(out/'candidates.npz',seeds=np.array([x[1] for x in seeds]),seed_scores=np.array([x[0] for x in seeds]),seed_psi=np.array([x[2] for x in seeds]),grown=np.concatenate(grown) if grown else np.empty((0,4,3)),grown_lengths=np.array([len(c) for c in grown]),metadata=json.dumps({'input':d.source,'profile_training':d.profile.metadata}))
 report=dict(engine='phasebuild-bounded-weak-map-trace',parameters=vars(args),residues=sum(map(len,chains)),fragments=len(chains),fragment_lengths=list(map(len,chains)),before_selection=sum(map(len,grown)),sources=sources,source_changed_during_run=any(fingerprint(s['path'])!=s for s in sources),seconds=time.time()-started,scope='Map and unrelated training target only; no reference/Buccaneer input. Development path, not a blind validation.')
 (out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
