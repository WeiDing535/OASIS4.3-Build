"""Disjoint-relation stability shrinkage for genuine KL-bounded P+ seeds.

A research robustness heuristic, not independent observation cross-validation.
Original likelihoods and conditional shapes stay fixed; no Rfree is accepted.
"""
import hashlib
import numpy as np
from oasis_py.relations import Relations
from oasis_py.clipper_likelihood import reflection_classes
from .pplus_conditional import ConditionalPplusResult,prepare_conditional_branches
from .pplus_trust import iterate_trusted_pplus,bernoulli_kl


def relation_groups(data,relations):
    groups=[]
    for i,j,k,sj,sk,shift in zip(relations.target,relations.left,relations.right,relations.left_sign,relations.right_sign,relations.shift):
        terms=sorted([(tuple(map(int,data.hkl[i])),1),(tuple(map(int,data.hkl[j])),-int(sj)),(tuple(map(int,data.hkl[k])),-int(sk))])
        sign=terms[0][1];key=tuple((h,c*sign) for h,c in terms)
        # Crystallographic symmetry translations have denominator24 or a divisor.
        phase=int(np.rint((-float(shift)*sign)/(2*np.pi)*24))%24
        digest=hashlib.sha256(repr((key,phase)).encode()).digest();groups.append((digest,int(digest[0]&1)))
    return groups


def stability_shrink(base,full,left,right,count_left,count_right,min_groups=4):
    arrays=[np.asarray(x,float) for x in [base,full,left,right]]
    if not arrays or any(x.shape!=arrays[0].shape or not np.isfinite(x).all() or np.any((x<0)|(x>1)) for x in arrays):raise ValueError('Invalid branch masses')
    b,f,l,r=arrays;dl=l-b;dr=r-b;df=f-b
    valid=(dl*dr>0)&(dl*df>0)&(np.asarray(count_left)>=min_groups)&(np.asarray(count_right)>=min_groups)
    denominator=abs(dl)+abs(dr);alpha=np.divide(2*np.minimum(abs(dl),abs(dr)),denominator,out=np.zeros_like(b),where=denominator>0);alpha=np.where(valid,alpha,0)
    return b+alpha*df,alpha


def phase_stability_pplus(data,observation_hl,reference,relations,*,samples=1440,kl_budget=.0005,cycles=100,damping=.25,min_groups=4):
    branches=prepare_conditional_branches(data,observation_hl,reference,samples=samples)
    group=relation_groups(data,relations);halves=[];counts=[]
    for half in [0,1]:
        mask=np.array([g[1]==half for g in group]);ids=np.flatnonzero(mask)
        fields=[np.asarray(getattr(relations,n))[mask] for n in ['target','left','right','left_sign','right_sign','shift']]
        graph=Relations(*fields,relations.kappa[mask]*2,relations.nreflections,{'half':half,'scaling':2,'relations':int(mask.sum())})
        groups=[set() for _ in data.f]
        for index in ids:groups[relations.target[index]].add(group[index][0])
        counts.append(np.array([len(x) for x in groups]));halves.append(iterate_trusted_pplus(data,branches,graph,kl_budget=kl_budget,cycles=cycles,damping=damping))
    full=iterate_trusted_pplus(data,branches,relations,kl_budget=kl_budget,cycles=cycles,damping=damping)
    base=np.exp(branches.log_mass[:,0]);p,alpha=stability_shrink(base,full.probability,halves[0].probability,halves[1].probability,*counts,min_groups=min_groups)
    first=p*branches.first[:,0]+(1-p)*branches.first[:,1];second=p*branches.second[:,0]+(1-p)*branches.second[:,1];phase=np.angle(first)%(2*np.pi);fom=np.minimum(abs(first),1)
    _,centric,allowed=reflection_classes(data);signed=np.real(first[centric]/np.exp(1j*allowed[centric]));phase[centric]=(allowed[centric]+np.where(signed<0,np.pi,0))%(2*np.pi)
    kl=bernoulli_kl(p,branches.log_mass)
    if np.max(kl)>kl_budget+1e-12:raise RuntimeError('Stability projection exceeded original KL budget')
    metadata={**full.metadata,'method':'disjoint-canonical-triplet-stability-shrunk-Pplus','approximation':'Robustness heuristic with fixed conditional branch shapes; not independent data cross-validation','minimum_distinct_groups_each_half':min_groups,'stable_update_rows':int((alpha>0).sum()),'mean_stability_weight':float(alpha.mean()),'mean_abs_probability_change':float(abs(p-base).mean()),'rows_probability_change_over_0p01':int((abs(p-base)>.01).sum()),'full_mean_abs_probability_change':full.metadata['mean_abs_probability_change'],'max_KL_from_original':float(kl.max()),'original_observation_HL_modified':False,'posterior_HL_exported':False,'half_graph_iterations':[r.metadata['iterations'] for r in halves],'group_assignment':'Canonical phase constraint modulo global sign and symmetry translation; HKL-stable SHA256 split','selection_uses_Rfree':False}
    return ConditionalPplusResult(phase,fom,p,first,second,full.history,metadata)
