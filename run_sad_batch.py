"""Frozen OASIS SAD -> PhaseBuild CLI benchmark with resumable bounded parallelism."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
import threading

PROJECT=Path(__file__).resolve().parent
OASIS=Path('/Volumes/Fast_4T/PAPER/Writing/4d_stem/gpt/oasis_modern')
PYTHON='/Applications/ccp4-9/bin/ccp4-python'
THREAD_ENV={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
RULES=dict(primary_endpoint='fixed_final_round_5_50_REFMAC_cycles; no reference-based or Rfree-based round selection',
    phase_route='corrected continuous SAD 5 internal cycles + native MLHL best-feedback DM 5 cycles; kernel_best.mtz',
    phase_dmin=2.8,map_amplitude='FP*FOM (applied once)',map_phase='PHIB',
    native='original prepared experimental mean amplitude, full measured resolution; no model coefficients',
    rounds=5,cycles_per_round=10,rotations=512,refine_count='min(512,max(96,ceil(sequence_ASU_residues/4)))',
    capacity='ceil(1.25*sequence_ASU_residues), bounded by 20000',
    thresholds=dict(r_free_max=.35,backbone_completeness_min=.8,backbone_precision_min=.9,joint_sequence_backbone_completeness_min=.8,matching_cutoff_angstrom=1.),
    denominator='All 100 registered; 84 pure proteins primary; 6 complexes protein-component exploratory; 10 unsupported nonprotein',
    missing_free_policy='No new free set generated; one-pass build only, not an Rfree-validated solution',
    reference_use='Post-build evaluation only, documented caoy b/2 gauge correction evaluation only',
    limitations='Retrospective corpus with supplied heavy sites and model-derived upstream composition/solvent metadata. Upstream OASIS used all observed reflections including free reflections. Rfree is therefore diagnostic, not a certified independent blind-test statistic. No hand or phase-route search.')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.{}.{}.tmp'.format(os.getpid(),threading.get_ident()))
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temp.replace(path)


def parse_sequence(path):
    records=[];header=None;parts=[];skip=False
    for line in Path(path).read_text(errors='strict').splitlines():
        if line.startswith('>'):
            if header is not None:records.append((header,''.join(parts)))
            header=line[1:].strip();parts=[];skip=header[:3].upper() in ('P1;','F1;')
        elif skip:skip=False
        elif line.strip():parts.append(''.join(line.split()).rstrip('*').upper())
    if header is not None:records.append((header,''.join(parts)))
    protein=[];excluded=[];normalizations=[]
    for header,seq in records:
        for source,target in [('(MSE)','M'),('(SEC)','X'),('(PYL)','X'),('(UNK)','X')]:
            if source in seq:normalizations.append(dict(header=header,source=source,target=target));seq=seq.replace(source,target)
        if set(seq)<=set('ACGTUX') or not seq:
            excluded.append(dict(header=header,reason='nucleic sequence alphabet'));continue
        if set(seq)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO'):
            raise ValueError('Unsupported sequence symbols: '+header+' '+str(set(seq)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO')))
        for ambiguous in 'BZJUO':
            if ambiguous in seq:normalizations.append(dict(header=header,source=ambiguous,target='X'));seq=seq.replace(ambiguous,'X')
        fields=header.split('|');chain_field=next((v for v in fields if re.match(r'^Chains?\s',v.strip())),None)
        chains=[str(i+1) for i in range(1)]
        if chain_field:
            chains=[v.strip() for v in re.sub(r'^Chains?\s+','',chain_field.strip()).split(',')]
        for chain in chains:protein.append((header+' declared_copy='+chain,seq))
    return protein,dict(excluded_records=excluded,normalizations=normalizations,copy_source='FASTA declared chain list; no target coordinate chain count used')


def prepare(root):
    import numpy as np
    import gemmi
    if (root/'manifest.json').exists():raise ValueError('Manifest exists; use run/status to resume')
    inventory=json.loads((root/'inventory.json').read_text())
    engine=root/'engine';engine.mkdir(exist_ok=False)
    shutil.copytree(PROJECT/'src',engine/'src',ignore=shutil.ignore_patterns('__pycache__'))
    for p in PROJECT.glob('run_*.py'):
        if p.name!='run_sad_batch.py':shutil.copyfile(p,engine/p.name)
    (engine/'models').mkdir()
    for name in ('1tqw_1p9_symmetric.npz','residue_1tqw_1p9.npz','rotamers_1tqw_v08.json'):
        shutil.copyfile(PROJECT/'models'/name,engine/'models'/name)
    shutil.copyfile(__file__,root/'runner_snapshot.py')
    cases=[];evaluation=[]
    for item in inventory:
        c=item['case'];cid=c['id'];folder=root/'inputs'/cid;folder.mkdir(parents=True)
        row=dict(id=cid,sample_class=c['sample_class'],supported=c['protein_residues']>0)
        # References belong exclusively to a separate post-build evaluator manifest.
        evaluation.append(dict(id=cid,reference=c['reference'],reference_sha256=digest(c['reference']),
                               fractional_origin_shift=[0,.5,0] if cid=='caoy_6e9c' else [0,0,0]))
        phase=Path(item['phase']);row.update(phase_source=str(phase),phase_sha256=digest(phase),sequence_source=item['sequence_files'][0],sequence_sha256=digest(item['sequence_files'][0]))
        if not row['supported']:
            row.update(preparation_status='unsupported_nonprotein',reason='PhaseBuild currently constructs protein only')
            cases.append(row);continue
        try:
            records,seq_audit=parse_sequence(row['sequence_source'])
            if not records:raise ValueError('No protein sequence')
            seq=folder/'sequence.fasta';seq.write_text(''.join('>'+h+'\n'+s+'\n' for h,s in records))
            row.update(sequence=str(seq),sequence_audit=seq_audit,sequence_records=len(records),sequence_residues=sum(len(s) for _,s in records))
            m=gemmi.read_mtz_file(str(phase));a=np.array(m,copy=False);labels=m.column_labels()
            f=a[:,labels.index('FP')];fom=a[:,labels.index('FOM')];phi=a[:,labels.index('PHIB')]
            if not np.all(np.isfinite(f)&np.isfinite(fom)&np.isfinite(phi)&(fom>=0)&(fom<=1)&(f>0)):
                raise ValueError('Invalid OASIS map coefficient input')
            output=gemmi.Mtz(with_base=True);output.spacegroup=m.spacegroup;output.set_cell_for_all(m.cell);output.add_dataset('OASIS_SAD_MAP')
            output.add_column('FWT','F');output.add_column('PHWT','P')
            output.set_data(np.column_stack([a[:,:3],f*fom,phi]).astype('float32'));output.update_reso()
            output.history=['Only OASIS FP*FOM and PHIB; no reference/model map coefficients']
            map_path=folder/'oasis_map.mtz';output.write_to_file(str(map_path));row['map']=str(map_path)
            native=gemmi.read_mtz_file(c['mtz']);nlabels=native.column_labels();arr=np.array(native,copy=False)
            flags=[col for col in native.columns if col.type=='I' and any(s in col.label.lower() for s in ('free','flag'))]
            has_free=bool(flags);row['has_free']=has_free
            if has_free:
                fp,sp=next(((f,s) for f,s in [('FP','SIGFP'),('F','SIGF'),('FMEAN','SIGFMEAN')] if f in nlabels and s in nlabels),(None,None))
                fl=flags[0].label;flag=arr[:,nlabels.index(fl)]
                if fp is not None:
                    if native.column_with_label(fp).type!='F' or native.column_with_label(sp).type!='Q':raise ValueError('Invalid native amplitude/sigma column types')
                    f=arr[:,nlabels.index(fp)];sig=arr[:,nlabels.index(sp)]
                elif all(label in nlabels for label in ['F(+)','SIGF(+)','F(-)','SIGF(-)']):
                    plus,ps,minus,ms=[arr[:,nlabels.index(label)] for label in ['F(+)','SIGF(+)','F(-)','SIGF(-)']]
                    vp=np.isfinite(plus)&np.isfinite(ps)&(plus>=0)&(ps>0)
                    vm=np.isfinite(minus)&np.isfinite(ms)&(minus>=0)&(ms>0)
                    f=np.where(vp&vm,(plus+minus)/2,np.where(vp,plus,np.where(vm,minus,np.nan)))
                    sig=np.where(vp&vm,np.sqrt(ps**2+ms**2)/2,np.where(vp,ps,np.where(vm,ms,np.nan)))
                    fp='arithmetic mean F(+),F(-), single mates retained';sp='propagated sigma of mean, single mates retained'
                else:raise ValueError('No experimental amplitude/standard deviation columns')
                finite_obs=np.isfinite(f)&np.isfinite(sig)&(f>=0)&(sig>0)
                valid_flag=np.isfinite(flag)&(flag>=0)&(flag==np.rint(flag))
                mask=finite_obs&valid_flag
                if not np.any(mask&(flag==0)) or not np.any(mask&(flag!=0)):raise ValueError('No existing free/work partition')
                output=gemmi.Mtz(with_base=True);output.spacegroup=native.spacegroup;output.set_cell_for_all(native.cell);output.add_dataset('EXPERIMENTAL_ONLY')
                for name,kind in [('F','F'),('SIGF','Q'),('FreeR_flag','I')]:output.add_column(name,kind)
                output.set_data(np.column_stack([arr[mask,:3],f[mask],sig[mask],flag[mask]]).astype('float32'));output.update_reso()
                native_path=folder/'native.mtz';output.write_to_file(str(native_path))
                row.update(native=str(native_path),native_source=c['mtz'],native_source_sha256=digest(c['mtz']),native_audit=dict(amplitude_source=fp,sigma_source=sp,flag_source=fl,test_value=0,flags_regenerated=False,input_rows=len(arr),retained_rows=int(mask.sum()),missing_observation_rows=int((~finite_obs).sum()),observed_rows_missing_valid_flag=int((finite_obs&~valid_flag).sum()),usable_free=int((mask&(flag==0)).sum()),usable_work=int((mask&(flag!=0)).sum())))
            else:row.update(native='',native_source=c['mtz'],native_source_sha256=digest(c['mtz']),no_free_reason='No supplied free-set flags; no new split generated')
            row['preparation_status']='ready'
            row['prepared_sha256']={str(p.relative_to(root)):digest(p) for p in folder.iterdir() if p.is_file()}
        except Exception as e:row.update(preparation_status='input_error',error=str(e),traceback=traceback.format_exc())
        cases.append(row)
    engine_hashes={str(p.relative_to(engine)):digest(p) for p in engine.rglob('*') if p.is_file()}
    save(root/'evaluation_manifest.json',evaluation)
    save(root/'manifest.json',dict(created=time.time(),rules=RULES,python=PYTHON,workers=8,cases=cases,engine_sha256=engine_hashes,runner_sha256=digest(root/'runner_snapshot.py')))
    summarize(root)
    print('Prepared',len(cases),'cases',flush=True)


def config_for(c,root):
    sys.path.insert(0,str(root/'engine/src'))
    from phasebuild.gui_workflow import defaults
    cfg=defaults(root/'engine');n=c['sequence_residues']
    cfg.update(mode='fresh',map=c['map'],native=c['native'],sequence=c['sequence'],output=str(root/'cases'/c['id']),
               refine=c['has_free'],rounds=5 if c['has_free'] else 1,cycles=10,map_mode='rebuild',
               refine_count=min(512,max(96,int(math.ceil(n/4.)))),max_residues=min(20000,int(math.ceil(n*1.25))))
    return cfg


def run_case(c,root,stage_timeout):
    state_path=root/'states'/(c['id']+'.json')
    if state_path.exists():
        state=json.loads(state_path.read_text())
        if state['status'] in ('completed','completed_no_free','no_trace','unsupported_nonprotein','input_error'):return state
        if state['status'] in ('failed','timeout'):return state
    start=time.time();state=dict(id=c['id'],sample_class=c['sample_class'],status='running',started=start,completed_steps=0,current='preflight')
    save(state_path,state)
    if c['preparation_status']!='ready':
        state.update(status=c['preparation_status'],error=c.get('error',c.get('reason')),seconds=0)
        save(state_path,state);return state
    try:
        sys.path.insert(0,str(root/'engine/src'))
        from phasebuild.gui_workflow import make_plan,save_selected_result,collect_results
        for relative,expected in c['prepared_sha256'].items():
            if digest(root/relative)!=expected:raise ValueError('Prepared input changed: '+relative)
        cfg=config_for(c,root);cfg,steps=make_plan(cfg,root/'engine',python=PYTHON)
        output=Path(cfg['output']);output.mkdir(parents=True,exist_ok=False)
        job=dict(config=cfg,steps=steps,status='running',completed=0,output=str(output))
        save(output/'gui_job.json',job)
        with (output/'pipeline.log').open('a',buffering=1) as log:
            for step in steps:
                state['current']=step['name'];save(state_path,state)
                log.write('\nSTEP '+step['name']+'\n'+json.dumps(step['command'],ensure_ascii=False)+'\n');log.flush()
                proc=subprocess.Popen(step['command'],cwd=str(root/'engine'),env=dict(os.environ,**THREAD_ENV),stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                state['child_pid']=proc.pid;save(state_path,state)
                try:code=proc.wait(timeout=stage_timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                    state['status']='timeout';raise RuntimeError('Stage time budget exceeded: '+step['name'])
                state.pop('child_pid',None)
                if code:raise RuntimeError(step['name']+' failed, exit '+str(code))
                state['completed_steps']+=1;job['completed']+=1;save(output/'gui_job.json',job);save(state_path,state)
                if step['kind']=='trace':
                    report=json.loads((Path(step['output'])/'build_report.json').read_text())
                    if report.get('residues',0)<3:
                        state['status']='no_trace';break
        if state['status']=='no_trace':job['status']='no_trace'
        else:
            results=save_selected_result(output)
            # Benchmark primary endpoint is the last planned round, not min Rfree.
            primary=output/results['last_model'];shutil.copyfile(primary,output/'primary_model.cif')
            state.update(primary_model=str(output/'primary_model.cif'),primary_metric=results['metrics'][-1] if results['metrics'] else None,
                         exploratory_best_metric=results['best_metric'],status='completed' if c['has_free'] else 'completed_no_free')
            job['status']='success'
        save(output/'gui_job.json',job)
    except Exception as e:
        state.update(status='timeout' if state.get('status')=='timeout' else 'failed',error=str(e),traceback=traceback.format_exc())
        if 'output' in locals() and (output/'gui_job.json').exists():
            job.update(status=state['status'],error=str(e));save(output/'gui_job.json',job)
    state.update(current='',seconds=time.time()-start,finished=time.time());state.pop('child_pid',None)
    save(state_path,state);return state


def evaluate_completed(c,state,root,evaluation):
    if 'primary_model' not in state:return state
    audit=root/'evaluation'/c['id'];audit.mkdir(parents=True,exist_ok=True)
    report_path=audit/'accuracy.json'
    if report_path.exists():accuracy=json.loads(report_path.read_text())
    else:
        import gemmi
        from phasebuild.evaluate import evaluate
        item=evaluation[c['id']]
        if digest(item['reference'])!=item['reference_sha256']:raise ValueError('Reference changed')
        reference=gemmi.read_structure(item['reference'])
        # Only protein CA residues enter the denominator for complexes.
        for chain in reference[0]:
            for i in range(len(chain)-1,-1,-1):
                if not gemmi.find_tabulated_residue(chain[i].name).is_amino_acid() or chain[i].find_atom('CA','*') is None:del chain[i]
        ref_path=audit/'protein_reference.cif';reference.make_mmcif_document().write_file(str(ref_path))
        model=state['primary_model'];shift=item['fractional_origin_shift']
        if any(shift):
            st=gemmi.read_structure(model);translation=st.cell.orthogonalize(gemmi.Fractional(*shift))
            for chain in st[0]:
                for residue in chain:
                    for atom in residue:atom.pos=gemmi.Position(atom.pos.x+translation.x,atom.pos.y+translation.y,atom.pos.z+translation.z)
            shifted=audit/'evaluation_origin_only.cif';st.make_mmcif_document().write_file(str(shifted));model=str(shifted)
        accuracy=evaluate(model,ref_path);accuracy['evaluation_only_origin_shift']=shift;save(report_path,accuracy)
    n=accuracy['reference_residues'];strict=accuracy['n_ca_c_correct_residues'];joint=accuracy['sequence_evaluation']['correct_type_and_backbone']
    metric=state.get('primary_metric');rfree=metric['r_free'] if metric else None;t=RULES['thresholds']
    solved=state['status']=='completed' and rfree is not None and rfree<=t['r_free_max'] and strict/n>=t['backbone_completeness_min'] and accuracy['n_ca_c_precision']>=t['backbone_precision_min'] and joint/n>=t['joint_sequence_backbone_completeness_min']
    state.update(evaluation=dict(built=accuracy['built_residues'],reference_protein_residues=n,ca_matched=accuracy['matched_residues'],strict_backbone_correct=strict,strict_completeness=strict/n,strict_precision=accuracy['n_ca_c_precision'],joint_sequence_backbone_correct=joint,joint_completeness=joint/n,protein_success_by_frozen_rule=solved,full_structure_success=solved and c['sample_class']=='protein'))
    save(root/'states'/(c['id']+'.json'),state);return state


def summarize(root):
    manifest=json.loads((root/'manifest.json').read_text());rows=[]
    for c in manifest['cases']:
        path=root/'states'/(c['id']+'.json')
        state=json.loads(path.read_text()) if path.exists() else dict(id=c['id'],sample_class=c['sample_class'],status=c['preparation_status'] if c['preparation_status']!='ready' else 'pending')
        m=state.get('primary_metric') or {};a=state.get('evaluation') or {}
        rows.append(dict(id=c['id'],sample_class=c['sample_class'],status=state['status'],completed_steps=state.get('completed_steps',0),current=state.get('current',''),r_work=m.get('r_work'),r_free=m.get('r_free'),built=a.get('built'),reference_protein_residues=a.get('reference_protein_residues'),strict_completeness=a.get('strict_completeness'),strict_precision=a.get('strict_precision'),joint_completeness=a.get('joint_completeness'),protein_success=a.get('protein_success_by_frozen_rule',False),full_structure_success=a.get('full_structure_success',False),seconds=state.get('seconds'),error=state.get('error','')))
    from collections import Counter
    status=dict(Counter(r['status'] for r in rows))
    report=dict(updated=time.time(),registered=100,status_counts=status,pure_protein_registered=84,complex_registered=6,unsupported_nonprotein_registered=10,pure_protein_success=sum(r['full_structure_success'] for r in rows),complex_protein_component_success=sum(r['protein_success'] and r['sample_class']=='protein_nucleic_complex' for r in rows),all100_full_structure_success=sum(r['full_structure_success'] for r in rows),rules=manifest['rules'],cases=rows)
    save(root/'summary.json',report)
    tmp=root/('results.csv.'+str(os.getpid())+'.tmp')
    with tmp.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    tmp.replace(root/'results.csv')
    md=['# OASIS SAD → PhaseBuild：100 套批量重构','',
        '后台状态：'+json.dumps(status,ensure_ascii=False)+'。未完成时以下成功数是累计值，不是最终成功率。','',
        '纯蛋白质成功：{} /84；复合物蛋白部分成功：{} /6；全100条中的完整蛋白结构成功：{} /100。10套非蛋白目标不支持，保留记录。'.format(report['pure_protein_success'],report['complex_protein_component_success'],report['all100_full_structure_success']),'',
        '主终点为固定第5轮、累计50个REFMAC周期，不按答案或最低Rfree挑轮次。缺少原始自由集时只构建一轮，单列且不算Rfree验证成功。判据：Rfree≤35%、严格N/CA/C覆盖率≥80%、精确率≥90%、序列与主链联合正确覆盖率≥80%，匹配阈值1 Å。','',
        '相位统一取既有修复后连续SAD五次内部迭代＋原生MLHL best反馈密度修饰五次的输出，使用FP×FOM/PHIB生成地图。只读取实验均值振幅及原自由集进行精修，没有读取源MTZ中的FC/PHIC作建模输入。','',
        '这是回顾性测试：上游重原子位点已提供，组成及溶剂含量沿用既有基准元数据；上游定相未排除自由反射，因此Rfree只作诊断，不能称为严格独立盲测。6套复合物仅评估蛋白部分，不算完整复合物解析。相位图为2.8 Å；精修可使用原始数据全部分辨率。','',
        '程序、模型库及输入哈希固定在manifest.json；构建与后验参考评估清单分开保存。进程数8，每个数值库1线程。每步预算6小时，失败与超时单列，已有输出保留。结果由后台程序自动更新。','',
        '| 案例 | 类型 | 状态 | Rwork | Rfree | 严格覆盖率 | 联合覆盖率 | 成功 |','|---|---|---|---:|---:|---:|---:|---|']
    pct=lambda v:'—' if v is None else '{:.2%}'.format(v)
    for r in rows:md.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(r['id'],r['sample_class'],r['status'],pct(r['r_work']),pct(r['r_free']),pct(r['strict_completeness']),pct(r['joint_completeness']),'是' if r['protein_success'] else '否'))
    (root/'REPORT.md').write_text('\n'.join(md)+'\n')
    return report


def run(root,workers,stage_timeout):
    import fcntl
    lock=(root/'batch.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=json.loads((root/'manifest.json').read_text())
    if digest(root/'runner_snapshot.py')!=manifest['runner_sha256']:raise ValueError('Frozen runner changed')
    for relative,expected in manifest['engine_sha256'].items():
        if digest(root/'engine'/relative)!=expected:raise ValueError('Frozen engine changed: '+relative)
    evaluation={c['id']:c for c in json.loads((root/'evaluation_manifest.json').read_text())}
    os.environ.update(THREAD_ENV);sys.path.insert(0,str(root/'engine/src'))
    save(root/'batch_process.json',dict(pid=os.getpid(),started=time.time(),workers=workers,status='running'))
    # Small targets start first; this is a size rule, not reference-phase quality.
    cases=sorted(manifest['cases'],key=lambda c:c.get('sequence_residues',0))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(run_case,c,root,stage_timeout):c for c in cases}
        while futures:
            done,_=concurrent.futures.wait(futures,timeout=15,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:
                c=futures.pop(f)
                try:
                    state=f.result();state=evaluate_completed(c,state,root,evaluation)
                    print(c['id'],state['status'],state.get('primary_metric'),state.get('evaluation'),flush=True)
                except Exception as e:
                    path=root/'states'/(c['id']+'.json');state=json.loads(path.read_text()) if path.exists() else dict(id=c['id'],sample_class=c['sample_class'])
                    state.update(evaluation_error=str(e),evaluation_traceback=traceback.format_exc());save(path,state)
                    print(c['id'],'evaluation_error',str(e),flush=True)
            summarize(root)
    report=summarize(root)
    save(root/'batch_process.json',dict(pid=os.getpid(),finished=time.time(),workers=workers,status='finished'))
    print('BATCH FINISHED',report['status_counts'],'pure protein successes',report['pure_protein_success'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','run','status']);parser.add_argument('--root',type=Path,default=PROJECT/'runs/sad100_20261019');parser.add_argument('--workers',type=int,default=8);parser.add_argument('--stage-timeout',type=int,default=21600)
    args=parser.parse_args();root=args.root.resolve();root.mkdir(parents=True,exist_ok=True)
    if not 1<=args.workers<=16:parser.error('workers must be 1..16')
    if args.action=='prepare':prepare(root)
    elif args.action=='run':run(root,args.workers,args.stage_timeout)
    else:
        report=summarize(root);print(json.dumps({k:v for k,v in report.items() if k not in ('cases','rules')},indent=2))


if __name__=='__main__':main()
