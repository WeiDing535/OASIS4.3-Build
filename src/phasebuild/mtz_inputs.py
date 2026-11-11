"""Explicit MTZ column selection and auditable copies for the bundled SAD route."""
import hashlib
import json
from contextlib import contextmanager
from pathlib import Path
import tempfile
import gemmi
import numpy as np
from .data import audit, fingerprint

SAD_KEYS=('fplus','sigfplus','fminus','sigfminus')
SAD_NAMES=('F(+)','SIGF(+)','F(-)','SIGF(-)')
SAD_TYPES=('G','L','G','L')


def mapping(text):
    result={}
    for item in text.split(',') if text else []:
        if '=' not in item:raise ValueError('列映射须为 key=列名')
        key,value=(s.strip() for s in item.split('=',1))
        if not value or key in result:raise ValueError('空或重复列映射：'+key)
        result[key]=value
    return result


def identifier(m,c):
    # Indices disambiguate duplicate labels and are validated again on every run.
    return '@'+str(list(m.columns).index(c))


def resolve(m,value,types):
    if value.startswith('@') and value[1:].isdigit():
        index=int(value[1:]);col=m.columns[index] if index<len(m.columns) else None
    else:
        matches=[c for c in m.columns if c.label==value]
        col=matches[0] if len(matches)==1 else None
    if col is None or col.type not in types:raise ValueError('缺少唯一且类型正确的 MTZ 列：'+value+'（应为 '+types+'）')
    return col


def candidate_groups(m,role):
    groups=[]
    for ds in m.datasets:
        cols=[(i,c) for i,c in enumerate(m.columns) if c.dataset_id==ds.id]
        if role=='sad':
            # Pair by matching labels, sigma spelling and dataset; never pair unrelated arrays.
            for i,c in cols:
                if c.type not in 'FG' or '(+)' not in c.label:continue
                minus=c.label.replace('(+)','(-)')
                def find(name,types):return next(((j,x) for j,x in cols if x.label==name and x.type in types),None)
                fm=find(minus,'FG');sp=find('SIG'+c.label,'QL');sm=find('SIG'+minus,'QL')
                if fm and sp and sm:
                    indices=[i,sp[0],fm[0],sm[0]]
                    groups.append(dict(mapping=dict(zip(SAD_KEYS,['@'+str(k) for k in indices])),
                        title=ds.dataset_name+' / '+', '.join(m.columns[k].label for k in indices)))
        else:
            for i,c in cols:
                if c.type!='F' or c.label.upper() in ('FWT','DELFWT','FC','FC_ALL') or '(+)' in c.label or '(-)' in c.label:continue
                sig=next(((j,x) for j,x in cols if x.label=='SIG'+c.label and x.type=='Q'),None)
                if sig:groups.append(dict(mapping={'f':'@'+str(i),'sigf':'@'+str(sig[0])},title=ds.dataset_name+' / '+c.label+', '+sig[1].label))
    return groups


def choices(path):
    info=audit(path);m=gemmi.read_mtz_file(str(path))
    for i,col in enumerate(info['columns']):col['id']='@'+str(i)
    info['sad_candidates']=candidate_groups(m,'sad');info['native_candidates']=candidate_groups(m,'native')
    info['free_candidates']=[dict(id='@'+str(i),label=c.label,dataset=m.dataset(c.dataset_id).dataset_name)
        for i,c in enumerate(m.columns) if c.type=='I' and ('FREE' in c.label.upper() or c.label.upper() in ('RFREE','R-FREE-FLAGS'))]
    return info


def select(m,text,role):
    selected=mapping(text);allowed=set(SAD_KEYS) if role=='sad' else {'f','sigf','free'}
    if set(selected)-allowed:raise ValueError('未知列映射：'+','.join(sorted(set(selected)-allowed)))
    required=SAD_KEYS if role=='sad' else ('f','sigf')
    if role=='native' and selected.get('f')=='__derive__' and selected.get('sigf')=='__derive__':
        return {k:v for k,v in selected.items() if k=='free'}
    if role=='native' and not any(k in selected for k in required):
        candidates=candidate_groups(m,role)
        if len(candidates)>1:raise ValueError('原始观测有多组振幅，请在输入界面选择 F / SIGF 列')
        if candidates:selected.update(candidates[0]['mapping'])
    if role=='sad' and not selected:
        candidates=candidate_groups(m,role)
        if len(candidates)!=1:raise ValueError('SAD 列未能唯一识别，请选择 F(+)、SIGF(+)、F(-)、SIGF(-)；当前只接受振幅，不把强度 I 当作 F')
        selected=candidates[0]['mapping']
    if role=='native' and not any(k in selected for k in required):return selected
    if any(k not in selected for k in required):raise ValueError('请完整选择 '+', '.join(required))
    return selected


def canonical(source,names,arrays,types,dataset_id=None):
    from .sad_bridge import keys
    keys(source)
    result=gemmi.Mtz(with_base=True);result.spacegroup=source.spacegroup
    result.set_cell_for_all(source.cell)
    original=source.dataset(dataset_id) if dataset_id is not None else next((d for d in source.datasets if d.id!=0),source.datasets[0])
    ds=result.add_dataset(original.dataset_name)
    ds.crystal_name=original.crystal_name;ds.project_name=original.project_name;ds.wavelength=original.wavelength;ds.cell=source.cell
    result.title=source.title
    # Preserve MTZ dataset metadata; this does NOT fill the mandatory SAD wavelength input.
    for name,kind in zip(names,types):result.add_column(name,kind,ds.id)
    result.set_data(np.column_stack([source.make_miller_array(),*arrays]).astype(np.float32));result.update_reso()
    return result


def build(cfg):
    """Pure in-memory selection; original files and existing flags never changed."""
    from .sad_bridge import keys,crystal
    sad=gemmi.read_mtz_file(cfg['sad']);native_path=cfg.get('native') or cfg['sad']
    native=gemmi.read_mtz_file(native_path);crystal(sad,native)
    # No silent reindexing, including for newly generated free flags.
    asu=gemmi.ReciprocalAsu(sad.spacegroup);ops=sad.spacegroup.operations()
    for m in (sad,native):
        if any(tuple(asu.to_asu(list(h),ops)[0])!=h for h in keys(m)):
            raise ValueError('需要已合并且位于同一 reciprocal ASU 的 HKL；请显式重索引')
    ss=select(sad,cfg.get('labels',''),'sad');ns=select(native,cfg.get('native_labels',''),'native')
    sad_cols=[resolve(sad,ss[k],t) for k,t in zip(SAD_KEYS,('FG','QL','FG','QL'))]
    if len({list(sad.columns).index(c) for c in sad_cols})!=4:raise ValueError('SAD 的正负振幅与误差必须选择四个不同的列')
    if len({c.dataset_id for c in sad_cols})!=1:raise ValueError('SAD 四列必须来自同一个 MTZ dataset')
    sa=[np.array(c,copy=True) for c in sad_cols]
    if 'f' in ns:
        ncols=[resolve(native,ns[k],t) for k,t in [('f','F'),('sigf','Q')]]
        if len({c.dataset_id for c in ncols})!=1:raise ValueError('平均 F / SIGF 必须来自同一个 dataset')
        f,sigma=[np.array(c,copy=True) for c in ncols];mean_source='selected original mean amplitude columns'
    else:
        # Derive means only from the SAME explicitly selected SAD pair, not a guessed array.
        if Path(native_path).resolve()!=Path(cfg['sad']).resolve():raise ValueError('原始观测文件缺少平均振幅 F / SIGF；请选列，或留空以使用 SAD 文件')
        fp,sp,fm,sm=sa
        plus=np.isfinite(fp)&np.isfinite(sp)&(fp>0)&(sp>0);minus=np.isfinite(fm)&np.isfinite(sm)&(fm>0)&(sm>0)
        both=plus&minus;f=np.full(len(fp),np.nan);sigma=f.copy()
        f[both]=(fp[both]+fm[both])/2;sigma[both]=np.hypot(sp[both],sm[both])/2
        for mask,amp,sig in [(plus&~minus,fp,sp),(minus&~plus,fm,sm)]:f[mask]=amp[mask];sigma[mask]=sig[mask]
        ncols=[];mean_source='derived from selected Bijvoet amplitudes: arithmetic mean; sigma=hypot(sigma+,sigma-)/2; measured singletons retained'
    free=ns.get('free');free_candidates=[(i,c) for i,c in enumerate(native.columns) if c.type=='I' and ('FREE' in c.label.upper() or c.label.upper() in ('RFREE','R-FREE-FLAGS'))]
    if free is None:
        if len(free_candidates)>1:raise ValueError('有多列自由集标记，请明确选择原自由集列')
        if free_candidates:free='@'+str(free_candidates[0][0])
    free_source=native
    if free is None:
        existing_sad=[(i,c) for i,c in enumerate(sad.columns) if c.type=='I' and ('FREE' in c.label.upper() or c.label.upper() in ('RFREE','R-FREE-FLAGS'))]
        if len(existing_sad)>1:raise ValueError('SAD 中有多列自由集，请选择 SAD 文件为原始观测并明确指定自由集列')
        if len(existing_sad)==1:free='@'+str(existing_sad[0][0]);free_source=sad
    generated=free is None
    if free=='__generate__' and free_candidates:raise ValueError('已有自由集，禁止重新生成；请选择已有标记')
    generated=generated or free=='__generate__'
    if generated:
        hkl=keys(native);valid=np.isfinite(f)&np.isfinite(sigma)&(f>=0)&(sigma>0)
        eligible=np.flatnonzero(valid)
        if len(eligible)<2:raise ValueError('没有足够观测生成工作集和自由集')
        # Stable across input row order; canonical ASU HKLs keep symmetry mates together.
        ordered=sorted(eligible,key=lambda i:hashlib.sha256((str(cfg['seed'])+':'+','.join(map(str,hkl[i]))).encode()).digest())
        flags=np.full(len(f),1 if cfg['free_value']==0 else 0,dtype=np.float32)
        count=min(len(ordered)-1,max(1,round(.05*len(ordered))));flags[ordered[:count]]=cfg['free_value'];free_col=None
    else:
        free_col=resolve(free_source,free,'I')
        if free_source is native:flags=np.array(free_col,copy=True)
        else:
            by_hkl=dict(zip(keys(sad),np.asarray(free_col)))
            flags=np.array([by_hkl.get(h,np.nan) for h in keys(native)],dtype=np.float32)
    # An existing SAD flag is checked by HKL against the chosen original flag, even when names differ.
    sad_flags=[c for c in sad.columns if c.type=='I' and ('FREE' in c.label.upper() or c.label.upper() in ('RFREE','R-FREE-FLAGS'))]
    if not generated and len(sad_flags)==1:
        lookup=dict(zip(keys(native),flags));sf=np.asarray(sad_flags[0])
        for h,old in zip(keys(sad),sf):
            value=lookup.get(h,np.nan)
            if np.isfinite(old) and np.isfinite(value) and old!=value:raise ValueError('SAD 与原始观测文件的自由集按 HKL 不一致')
    def descriptor(m,c):return dict(index=list(m.columns).index(c),label=c.label,type=c.type,dataset_id=c.dataset_id,dataset=m.dataset(c.dataset_id).dataset_name)
    report=dict(sad=fingerprint(cfg['sad']),native=fingerprint(native_path),native_omitted=not bool(cfg.get('native')),
        sad_columns={k:descriptor(sad,c) for k,c in zip(SAD_KEYS,sad_cols)},
        native_columns={k:descriptor(native,c) for k,c in zip(('f','sigf'),ncols)},mean_source=mean_source,
        free_source=fingerprint(cfg['sad'] if free_source is sad else native_path),
        free_column=descriptor(free_source,free_col) if free_col is not None else None,
        flags_generated=generated,free_value=cfg['free_value'],seed=cfg['seed'] if generated else None,
        free_policy='new deterministic 5% HKL partition; not the original experimental free set; Rfree cannot be compared directly to prior runs' if generated else 'existing original flags copied unchanged by HKL',
        intensities_converted=False,input_files_modified=False)
    return canonical(sad,SAD_NAMES,sa,SAD_TYPES,sad_cols[0].dataset_id),canonical(native,('F','SIGF','FreeR_flag'),(f,sigma,flags),('F','Q','I'),ncols[0].dataset_id if ncols else sad_cols[0].dataset_id),report


@contextmanager
def checked_config(cfg):
    """Temporary canonical views for legacy readers; nothing written to user inputs/output."""
    sad,native,report=build(cfg)
    with tempfile.TemporaryDirectory(prefix='oasis-mtz-check-') as temp:
        root=Path(temp);sad.write_to_file(str(root/'sad.mtz'));native.write_to_file(str(root/'native.mtz'))
        yield {**cfg,'sad':str(root/'sad.mtz'),'native':str(root/'native.mtz'),'labels':'','native_labels':'','input_columns':False},report


def prepare(cfg,folder):
    from .sad_bridge import prepare_sad
    folder=Path(folder);sad,native,report=build(cfg)
    raw=folder/'selected_sad.mtz';obs=folder/'native.mtz'
    sad.write_to_file(str(raw));native.write_to_file(str(obs))
    result=prepare_sad(raw,obs,folder/'sad.mtz',cfg['free_value']);result['column_selection']=report
    result['flags_regenerated']=report['flags_generated']
    (folder/'sad.json').write_text(json.dumps(result,indent=2)+'\n')
    (folder/'columns.json').write_text(json.dumps(report,indent=2)+'\n')
    return result
