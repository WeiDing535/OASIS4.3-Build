"""Frozen Buccaneer/REFMAC comparison on the same 100 OASIS SAD datasets."""
import argparse,concurrent.futures,csv,hashlib,importlib.util,json,os,shutil,signal,subprocess,sys,time,traceback,threading
from pathlib import Path
PROJECT=Path(__file__).resolve().parent
PREVIOUS=Path('/Volumes/Fast_4T/GPT_WORK/phasebuild/runs/sad100_20261019')
PYTHON='/Applications/ccp4-9/bin/ccp4-python'
CCP4=Path('/Applications/ccp4-9')
ENV={k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.%s.%s.tmp'%(os.getpid(),threading.get_ident()));t.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n');t.replace(p)

def prepare(root):
 import gemmi,numpy as np
 if (root/'manifest.json').exists():raise ValueError('Already frozen')
 original=load(PREVIOUS/'manifest.json');inventory={c['id']:c['case'] for c in load(PREVIOUS/'inventory.json')}
 shutil.copytree(PREVIOUS/'engine',root/'engine',ignore=shutil.ignore_patterns('__pycache__'))
 ref=root/'engine/src/phasebuild/refinement.py';s=ref.read_text()
 assert 'not all(0.<=v<=1. for v in values[:2])' in s
 s=s.replace('not all(0.<=v<=1. for v in values[:2])','not all(v>=0. for v in values[:2])').replace('free_value=0,timeout=300):','free_value=0,timeout=3600):');ref.write_text(s)
 shutil.copyfile(__file__,root/'runner_snapshot.py')
 cases=[]
 for old in original['cases']:
  c={k:old[k] for k in ['id','sample_class','supported','preparation_status']};cid=c['id']
  if not c['supported']:cases.append(c);continue
  for rel,h in old['prepared_sha256'].items():assert sha(PREVIOUS/rel)==h
  assert sha(old['phase_source'])==old['phase_sha256']
  seq=Path(old['sequence']);native=Path(old['native']) if old['has_free'] else None
  c.update(sequence=str(seq),sequence_sha256=sha(seq),sequence_residues=old['sequence_residues'],native=str(native) if native else '',native_sha256=sha(native) if native else None,has_free=old['has_free'],phase_source=old['phase_source'],phase_sha256=old['phase_sha256'],build_semet='Se' in inventory[cid]['heavy_elements'])
  mtz=gemmi.read_mtz_file(old['phase_source']);data=np.array(mtz,copy=False);labels=mtz.column_labels()
  f=data[:,labels.index('FP')];sig=data[:,labels.index('SIGFP')];phi=data[:,labels.index('PHIB')];fom=data[:,labels.index('FOM')]
  out=gemmi.Mtz(with_base=True);out.spacegroup=mtz.spacegroup;out.set_cell_for_all(mtz.cell);out.add_dataset('OASIS_SAD_ONLY')
  for label,kind in [('F','F'),('SIGF','Q'),('PHWT','P'),('FOM','W'),('FWT','F')]:out.add_column(label,kind)
  out.set_data(np.column_stack([data[:,:3],f,sig,phi,fom,f*fom]).astype('float32'));out.update_reso()
  folder=root/'inputs'/cid;folder.mkdir(parents=True);out.write_to_file(str(folder/'initial.mtz'))
  c.update(initial_mtz=str(folder/'initial.mtz'),initial_sha256=sha(folder/'initial.mtz'));cases.append(c)
 rules=dict(original['rules']);rules.update(builder='CCP4 cbuccaneer 1.6.12',buccaneer_cycles=[3,2,2,2,2],correlation_mode=[False,True,True,True,True],phase_schedule='Initial same OASIS FP*FOM/PHIB; subsequent maps and models are previous Buccaneer/REFMAC outputs',building_free_column='Not passed to builder, as both builders use supplied all-reflection maps; native refinement always preserves existing free flags',semet='Enable -build-semet when experiment metadata contains Se; no target-coordinate decision',adapter_changes='Same REFMAC keywords/cycles/data; accept finite nonnegative R>1 instead of rejecting; timeout 3600s instead of prior 300s. Prior frozen PhaseBuild results unchanged.',missing_free_policy='caoy: one initial 3-cycle Buccaneer build, no refinement/no Rfree; excluded from validated successes',building_timeout_seconds=21600)
 reference=CCP4/'lib/data/reference_structures'
 provenance=[CCP4/'bin/cbuccaneer',CCP4/'bin/refmac5',reference/'reference-1tqw.pdb',reference/'reference-1tqw.mtz',CCP4/'html/cbuccaneer.html']
 save(root/'evaluation_manifest.json',load(PREVIOUS/'evaluation_manifest_v2.json'))
 save(root/'manifest.json',dict(created=time.time(),baseline=str(PREVIOUS),rules=rules,cases=cases,engine_sha256={str(p.relative_to(root/'engine')):sha(p) for p in (root/'engine').rglob('*') if p.is_file()},executables_and_reference={str(p):sha(p) for p in provenance},runner_sha256=sha(root/'runner_snapshot.py')))
 summarize(root)
 print('Prepared 100; supported',sum(c['supported'] for c in cases),flush=True)

def bucc_command(c,mtz,folder,number,previous=None):
 ref=CCP4/'lib/data/reference_structures/reference-1tqw'
 cmd=[str(CCP4/'bin/cbuccaneer'),'-mtzin',str(mtz),'-mtzin-ref',str(ref.with_suffix('.mtz')),'-pdbin-ref',str(ref.with_suffix('.pdb')),'-colin-ref-fo','/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]','-colin-ref-hl','/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]','-colin-fo','/*/*/[F,SIGF]','-colin-phifom','/*/*/[PHWT,FOM]','-colin-fc','/*/*/[FWT,PHWT]','-seqin',c['sequence'],'-pdbout',str(folder/'model.pdb'),'-xmlout',str(folder/'buccaneer.xml')]
 cmd+=['-'+s for s in ['find','grow','join','link','sequence','correct','filter','ncsbuild','prune','rebuild']]
 cmd+=['-cycles','3' if number==1 else '2','-jobs','1','-no-correlation-mode' if number==1 else '-correlation-mode']
 if previous:cmd+=['-pdbin',str(previous)]
 if c['build_semet']:cmd+=['-build-semet']
 return cmd

def execute(cmd,folder,logname,timeout,state,statefile):
 folder.mkdir(parents=True,exist_ok=False);save(folder/'command.json',cmd)
 env=dict(os.environ,**ENV,CCP4=str(CCP4),CLIBD=str(CCP4/'lib/data'),CLIBD_MON=str(CCP4/'lib/data/monomers')+'/')
 with (folder/logname).open('w') as log:
  proc=subprocess.Popen(cmd,cwd=str(folder),env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  state['child_pid']=proc.pid;save(statefile,state)
  try:code=proc.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGTERM)
   try:proc.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
   raise RuntimeError('Step exceeded time limit '+str(timeout))
  if code:raise RuntimeError('Program exit '+str(code)+': '+str(folder/logname))
 state.pop('child_pid',None)

def run_case(c,root):
 import gemmi
 path=root/'states'/(c['id']+'.json')
 if path.exists() and load(path)['status'] not in ('running','queued'):return load(path)
 state=dict(id=c['id'],sample_class=c['sample_class'],status='running',started=time.time(),current='preflight',completed_rounds=0,stages=[]);save(path,state)
 if not c['supported']:
  state.update(status='unsupported_nonprotein',current='',seconds=0);save(path,state);return state
 try:
  for name,digest in [('phase_source',c['phase_sha256']),('sequence',c['sequence_sha256']),('initial_mtz',c['initial_sha256'])]:assert sha(c[name])==digest,name
  if c['has_free']:assert sha(c['native'])==c['native_sha256']
  case=root/'cases'/c['id'];case.mkdir(parents=True,exist_ok=False)
  mtz=Path(c['initial_mtz']);previous=None;rounds=5 if c['has_free'] else 1
  for number in range(1,rounds+1):
   build=case/('round_%02d_build'%number);state['current']='Buccaneer round %d'%number;save(path,state)
   cmd=bucc_command(c,mtz,build,number,previous);execute(cmd,build,'buccaneer.log',21600,state,path)
   st=gemmi.read_structure(str(build/'model.pdb'));st.make_mmcif_document().write_file(str(build/'model.cif'))
   residues=sum(1 for ch in st[0] for r in ch if r.find_atom('CA','*')) if len(st) else 0
   save(build/'build_report.json',dict(status='success',engine='external-Buccaneer',command=cmd,residues=residues,cycles=3 if number==1 else 2))
   if residues==0:state.update(status='no_trace',current='');break
   if not c['has_free']:
    state.update(status='completed_no_free',completed_rounds=1,primary_model=str(build/'model.cif'),primary_metric=None);break
   refine=case/('round_%02d_refine'%number);state['current']='REFMAC round %d'%number;save(path,state)
   # run_refine creates its own output folder, so the outer log has a separate path.
   logfolder=case/('round_%02d_refine_runner'%number)
   cmd=[PYTHON,'-u',str(root/'engine/run_refine.py'),c['native'],'--model',str(build/'model.cif'),'--output',str(refine),'--cycles','10']
   execute(cmd,logfolder,'runner.log',3900,state,path)
   report=load(refine/'refinement_report.json');assert report['free_flags_preserved'] and report['statistics']['final']['cycle']==10
   metric=dict(stage=refine.name,r_work=report['statistics']['final']['r_work'],r_free=report['statistics']['final']['r_free'],free_preserved=True)
   state['stages'].append(metric);state['completed_rounds']=number;save(path,state)
   mtz=refine/'refined.mtz';previous=refine/'refined.pdb'
   if number==5:
    shutil.copyfile(refine/'model.cif',case/'primary_model.cif');state.update(status='completed',primary_model=str(case/'primary_model.cif'),primary_metric=metric)
 except Exception as e:
  state.update(status='failed',error=str(e),traceback=traceback.format_exc())
 state.update(current='',finished=time.time(),seconds=time.time()-state['started']);state.pop('child_pid',None);save(path,state);return state

def evaluate_case(c,state,root,refs):
 if not state.get('primary_model'):return state
 import gemmi
 from phasebuild.evaluate import evaluate
 item=refs[c['id']];assert sha(item['reference'])==item['sha256']
 out=root/'evaluation'/c['id'];out.mkdir(parents=True,exist_ok=True)
 previous=PREVIOUS/'evaluation'/c['id']/'protein_reference_v2.cif'
 if previous.exists():shutil.copyfile(previous,out/'protein_reference.cif')
 else:
  st=gemmi.read_structure(item['reference']);source_cell=st.cell;target=gemmi.read_mtz_file(c['initial_mtz'])
  assert st.find_spacegroup().number==target.spacegroup.number and source_cell.is_similar(target.cell,.01,1.)
  for ch in st[0]:
   for i in range(len(ch)-1,-1,-1):
    res=ch[i]
    if not gemmi.find_tabulated_residue(res.name).is_amino_acid() or res.find_atom('CA','*') is None:del ch[i];continue
    best={}
    for a in res:
     if a.name not in best or a.occ>best[a.name].occ:best[a.name]=a.clone()
    for j in range(len(res)-1,-1,-1):del res[j]
    for a in best.values():a.altloc='\x00';a.pos=target.cell.orthogonalize(source_cell.fractionalize(a.pos));res.add_atom(a)
  st.cell=target.cell;st.spacegroup_hm=target.spacegroup.hm;st.make_mmcif_document().write_file(str(out/'protein_reference.cif'))
 model=state['primary_model'];shift=item['fractional_origin_shift']
 if any(shift):
  st=gemmi.read_structure(model);t=st.cell.orthogonalize(gemmi.Fractional(*shift))
  for ch in st[0]:
   for res in ch:
    for a in res:a.pos=gemmi.Position(a.pos.x+t.x,a.pos.y+t.y,a.pos.z+t.z)
  st.make_mmcif_document().write_file(str(out/'evaluation_only_shift.cif'));model=str(out/'evaluation_only_shift.cif')
 a=evaluate(model,out/'protein_reference.cif');a['evaluation_only_origin_shift']=shift;save(out/'accuracy.json',a)
 n=a['reference_residues'];joint=a['sequence_evaluation']['correct_type_and_backbone'];metric=state.get('primary_metric');t=load(root/'manifest.json')['rules']['thresholds']
 passed=state['status']=='completed' and metric['r_free']<=t['r_free_max'] and a['n_ca_c_completeness']>=t['backbone_completeness_min'] and a['n_ca_c_precision']>=t['backbone_precision_min'] and joint/n>=t['joint_sequence_backbone_completeness_min']
 state['evaluation']=dict(built=a['built_residues'],reference_protein_residues=n,strict_completeness=a['n_ca_c_completeness'],strict_precision=a['n_ca_c_precision'],joint_completeness=joint/n,protein_success=bool(passed),full_structure_success=bool(passed and c['sample_class']=='protein'),revision='same-residue-aware-reference-as-PhaseBuild')
 save(root/'states'/(c['id']+'.json'),state);return state

def summarize(root):
 from collections import Counter
 m=load(root/'manifest.json');rows=[];baseline={d['id']:d for d in load(PREVIOUS/'summary.json')['cases']}
 for c in m['cases']:
  p=root/'states'/(c['id']+'.json');s=load(p) if p.exists() else dict(status='pending' if c['supported'] else 'unsupported_nonprotein')
  a=s.get('evaluation') or {};metric=s.get('primary_metric') or {};b=baseline[c['id']]
  rows.append(dict(id=c['id'],sample_class=c['sample_class'],status=s['status'],current=s.get('current',''),completed_rounds=s.get('completed_rounds',0),r_work=metric.get('r_work'),r_free=metric.get('r_free'),built=a.get('built'),strict_completeness=a.get('strict_completeness'),strict_precision=a.get('strict_precision'),joint_completeness=a.get('joint_completeness'),success=a.get('protein_success') if s['status']=='completed' else None,phasebuild_status=b['status'],phasebuild_r_free=b.get('r_free'),phasebuild_strict_completeness=b.get('strict_completeness'),error=s.get('error',''),evaluation_error=s.get('evaluation_error','')))
 summary=dict(updated=time.time(),registered=100,status_counts=dict(Counter(r['status'] for r in rows)),pure_protein_success=sum(r['success'] is True and r['sample_class']=='protein' for r in rows),complex_protein_component_success=sum(r['success'] is True and r['sample_class']=='protein_nucleic_complex' for r in rows),rules=m['rules'],cases=rows)
 save(root/'summary.json',summary)
 tmp=root/('results.csv.%s.tmp'%os.getpid())
 with tmp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 tmp.replace(root/'results.csv')
 md=['# 同一100套OASIS SAD相位：Buccaneer 对照','',
 '状态：'+json.dumps(summary['status_counts'],ensure_ascii=False)+'。未全部结束前为累计结果，不能当作最终成功率。','',
 '纯蛋白质达标：%d /84；复合物蛋白部分达标：%d /6。10套非蛋白不支持，caoy缺原自由集单列。'%(summary['pure_protein_success'],summary['complex_protein_component_success']),'',
 '起始OASIS相位、序列、native观测与自由集复用上一批文件与哈希，未重新定相或读取目标模型作建模输入。首次Buccaneer默认似然目标3循环，后4轮各2循环相关目标；每轮10个REFMAC周期，共50周期，固定第5轮作主终点。每次后续建模只用自身上一轮模型与更新图。依据本机CCP4文档选择模式；Se实验按元数据启用-build-semet。','',
 '成功标准与PhaseBuild相同：Rfree≤35%、严格N/CA/C覆盖率≥80%、精确率≥90%、序列与严格主链联合覆盖率≥80%，1 Å阈值。复合物仅评估蛋白部分，不能声称完整复合物解析。','',
 '适配器差别公开记录：仅放宽有限非负R>1的错误拦截，并将REFMAC墙钟超时从300秒增至3600秒；精修关键词、反射、自由集和周期数不变。两建模器内部循环与目标函数不同，不等同CPU预算完全相同。旧PhaseBuild失败记录没有被覆盖。','',
 '已知重原子、历史上游组成/溶剂元数据、上游定相包含自由反射等限制仍适用；Rfree为诊断，不是严格独立盲测。','',
 '| 案例 | 类型 | 状态 | Buccaneer Rfree | PhaseBuild Rfree | 严格覆盖率 | 联合覆盖率 | 达标 |','|---|---|---|---:|---:|---:|---:|---|']
 pct=lambda v:'—' if v is None else '{:.2%}'.format(v)
 for r in rows:md.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(r['id'],r['sample_class'],r['status'],pct(r['r_free']),pct(r['phasebuild_r_free']),pct(r['strict_completeness']),pct(r['joint_completeness']),'—' if r['success'] is None else '是' if r['success'] else '否'))
 (root/'REPORT.md').write_text('\n'.join(md)+'\n');return summary

def run(root,workers):
 import fcntl
 lock=(root/'batch.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 m=load(root/'manifest.json');assert sha(root/'runner_snapshot.py')==m['runner_sha256']
 for rel,h in m['engine_sha256'].items():assert sha(root/'engine'/rel)==h
 for p,h in m['executables_and_reference'].items():assert sha(p)==h
 os.environ.update(ENV);sys.path.insert(0,str(root/'engine/src'));refs={c['id']:c for c in load(root/'evaluation_manifest.json')}
 save(root/'batch_process.json',dict(pid=os.getpid(),status='running',started=time.time(),workers=workers))
 with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
  futures={pool.submit(run_case,c,root):c for c in sorted(m['cases'],key=lambda c:c.get('sequence_residues',0))}
  while futures:
   done,_=concurrent.futures.wait(futures,timeout=15,return_when=concurrent.futures.FIRST_COMPLETED)
   for f in done:
    c=futures.pop(f)
    try:
     s=evaluate_case(c,f.result(),root,refs);print(c['id'],s['status'],s.get('primary_metric'),s.get('evaluation'),flush=True)
    except Exception as e:
     p=root/'states'/(c['id']+'.json');s=load(p);s.update(evaluation_error=str(e),evaluation_traceback=traceback.format_exc());save(p,s);print(c['id'],'evaluation_error',str(e),flush=True)
   summarize(root)
 summarize(root);save(root/'batch_process.json',dict(pid=os.getpid(),status='finished',finished=time.time(),workers=workers));print('FINISHED',flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','status']);p.add_argument('--root',type=Path,default=PROJECT/'runs/buccaneer_sad100_20261019');p.add_argument('--workers',type=int,default=8);a=p.parse_args();r=a.root.resolve();r.mkdir(parents=True,exist_ok=True)
 if a.action=='prepare':prepare(r)
 elif a.action=='run':run(r,a.workers)
 else:print(json.dumps({k:v for k,v in summarize(r).items() if k not in ['cases','rules']},indent=2))
if __name__=='__main__':main()
