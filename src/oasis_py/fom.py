# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Case-balanced empirical FOM calibration for FINAL MAP WEIGHTS only.

Learned from other solved structures; not a physical likelihood or a new phase
solver. Apply requires an explicit target id and refuses training-case reuse.
The fitted probability shape is unknown, so calibrated output omits HL/PPLUS.
"""
from dataclasses import dataclass
import json
import numpy as np
from .phasing import PhaseResult


def isotonic(values,weights):
    blocks=[]
    for i,(v,w) in enumerate(zip(values,weights)):
        blocks.append([i,i,float(v*w),float(w)])
        while len(blocks)>1 and blocks[-2][2]/blocks[-2][3]>blocks[-1][2]/blocks[-1][3]:
            b=blocks.pop();a=blocks.pop();blocks.append([a[0],b[1],a[2]+b[2],a[3]+b[3]])
    out=np.empty(len(values))
    for lo,hi,total,weight in blocks:out[lo:hi+1]=total/weight
    return np.clip(out,0,1)


@dataclass
class FomCalibration:
    model: dict

    @classmethod
    def fit(cls,records,method,initial_scale_fit_version='legacy_bounded'):
        """records: case id, raw FOM, observed cos(error), centric mask.

        Target correct phases are training labels, so fitting must NEVER use
        the evaluated target or its sequence family. Each case has equal mass.
        """
        if not records or len({r['id'] for r in records})!=len(records):raise ValueError('Calibration requires unique, nonempty training cases')
        for r in records:
            f=np.asarray(r['fom']);y=np.asarray(r['cosine']);c=np.asarray(r['centric'])
            if not len(f) or f.shape!=y.shape or f.shape!=c.shape or not np.isfinite(f).all() or not np.isfinite(y).all() or np.any((f<0)|(f>1)) or np.any(abs(y)>1+1e-10):raise ValueError('Invalid calibration training arrays')
        tables={}
        for name in ['all','centric','acentric']:
            sums=np.zeros(100);xs=np.zeros(100);weights=np.zeros(100)
            for r in records:
                f=np.asarray(r['fom']);y=np.asarray(r['cosine']);c=np.asarray(r['centric'],bool)
                keep=np.ones(len(f),bool) if name=='all' else c if name=='centric' else ~c
                ix=np.minimum((f[keep]*100).astype(int),99);w=np.full(len(ix),1/len(f))
                weights+=np.bincount(ix,weights=w,minlength=100)
                sums+=np.bincount(ix,weights=w*y[keep],minlength=100)
                xs+=np.bincount(ix,weights=w*f[keep],minlength=100)
            good=weights>0
            if good.sum()<3:continue
            x=xs[good]/weights[good];y=isotonic(sums[good]/weights[good],weights[good])
            tables[name]={'x':x.tolist(),'y':y.tolist()}
        if 'all' not in tables:raise ValueError('Insufficient FOM calibration data')
        return cls({'schema':1,'method':method,'initial_scale_fit_version':initial_scale_fit_version,'training_cases':[r['id'] for r in records],'tables':tables,'usage':'final map weights only; no calibrated HL or branch probabilities'})

    def predict(self,fom,centric):
        fom=np.asarray(fom,dtype=float);centric=np.asarray(centric,dtype=bool)
        if fom.shape!=centric.shape or not np.isfinite(fom).all() or np.any((fom<0)|(fom>1)):raise ValueError('Invalid FOM calibration inputs')
        result=np.empty_like(fom)
        for name,keep in [('centric',centric),('acentric',~centric)]:
            table=self.model['tables'].get(name,self.model['tables']['all'])
            result[keep]=np.interp(fom[keep],table['x'],table['y'])
        # Uniform initial phase cannot acquire phase information from calibration.
        result[fom<=1e-8]=0
        return np.clip(result,0,1)

    def apply(self,data,result,target_id):
        if not target_id:raise ValueError('An explicit target id is required for calibration provenance')
        if target_id.lower() in {s.lower() for s in self.model['training_cases']}:raise ValueError('Target occurs in calibration training set; use an excluded-target model')
        if result.metadata.get('method')!=self.model['method']:raise ValueError('Calibration model does not match the phase method')
        if result.metadata.get('initial_scale_fit_version','legacy_bounded')!=self.model.get('initial_scale_fit_version','legacy_bounded'):raise ValueError('FOM model and phase scale-fit versions differ; revalidation required')
        if self.model['method']=='native_classical_dm_experimental':
            if result.metadata.get('initial_method') not in {None,'sad_likelihood'}:raise ValueError('Existing FOM model requires the continuous SAD initial method')
            if result.metadata.get('bias_correction') or result.metadata.get('anisotropy',{}).get('status','disabled')!='disabled' or result.metadata.get('ncs_operator_count',0):
                raise ValueError('This calibration model was fitted without unsupported corrections, anisotropy or NCS; configuration mismatch')
            if result.metadata.get('cycles') and len(result.metadata['cycles'])!=5:raise ValueError('Calibration requires the five-cycle native DM configuration')
        if self.model.get('ensemble_stage') and result.metadata.get('ensemble_stage')!=self.model['ensemble_stage']:raise ValueError('Calibration ensemble stage mismatch')
        if self.model.get('ensemble_stage'):
            for component in result.metadata.get('final_component_metadata',[]):
                if component.get('anisotropy',{}).get('status','disabled')!='disabled' or component.get('ncs_operator_count',0) or len(component.get('cycles',[]))!=5:raise ValueError('Unsupported ensemble calibration configuration')
        calibrated=self.predict(result.fom,data.centric)
        return PhaseResult(result.phase.copy(),calibrated,None,None,{**result.metadata,'fom_calibration':{'training_cases':self.model['training_cases'],'target_id':target_id,'initial_scale_fit_version':self.model.get('initial_scale_fit_version','legacy_bounded'),'purpose':'final_map_weights_only','HL_omitted':True,'PPLUS_omitted':True,'phase_unchanged':True,'homology_exclusion_must_be_checked_by_caller':True}})

    def save(self,path):
        with open(path,'w') as f:json.dump(self.model,f,indent=2,allow_nan=False)

    @classmethod
    def read(cls,path):
        with open(path) as f:model=json.load(f)
        if model.get('schema')!=1 or not model.get('tables'):raise ValueError('Unsupported calibration model')
        return cls(model)
