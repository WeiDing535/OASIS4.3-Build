"""Unit tests for the work-only Parrot + trace-feedback density-modification layer.

These tests use synthetic MTZs and never execute Parrot, Buccaneer or REFMAC.
"""
from pathlib import Path
import json
import numpy as np
import gemmi
import pytest

from phasebuild import parrot_dm, trace_feedback, ccp4_pipeline


def _mtz(path, rows, cell=(40.0, 40.0, 40.0, 90.0, 90.0, 90.0), sg='P 1', extra=()):
    m=gemmi.Mtz(with_base=True)
    m.spacegroup=gemmi.SpaceGroup(sg)
    m.set_cell_for_all(gemmi.UnitCell(*cell))
    m.add_dataset('test')
    for label,kind in extra:
        m.add_column(label,kind)
    m.set_data(np.asarray(rows,np.float32))
    m.write_to_file(str(path))
    return path


def _hl(kappa, phase_deg):
    r=np.deg2rad(phase_deg)
    return [kappa*np.cos(r),kappa*np.sin(r),0.0,0.0]


def test_to_phifom_recovers_phase_and_monotone_fom(tmp_path):
    rows=[]
    for i,(kappa,phase) in enumerate([(0.0,0.0),(0.5,90.0),(5.0,180.0),(60.0,270.0)]):
        rows.append([i+1,0,1,10.0,1.0]+_hl(kappa,phase))
    src=_mtz(tmp_path/'hl.mtz',rows,extra=[('FP','F'),('SIGFP','Q'),
        ('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')])
    dst=tmp_path/'phifom.mtz'
    info=parrot_dm.to_phifom(src,dst)
    m=gemmi.read_mtz_file(str(dst));labels=m.column_labels()
    for label in ('F','SIGF','FWT','PHWT','PHIB','FOM'):
        assert label in labels
    phases=np.asarray(m.column_with_label('PHIB'))
    foms=np.asarray(m.column_with_label('FOM'))
    assert np.allclose(phases%360,[0.0,90.0,180.0,270.0],atol=1e-4)
    assert np.all(np.diff(foms)>0),'FOM must increase with the HL magnitude'
    assert 0.0<=float(foms.min()) and float(foms.max())<=1.0

    # FWTs are the FOM-weighted amplitudes actually used downstream.
    fwt=np.asarray(m.column_with_label('FWT'))
    assert np.allclose(fwt,10.0*foms,atol=1e-4)
    assert info['mean_fom']==pytest.approx(float(foms.mean()))


def test_work_only_input_excludes_free_rows(tmp_path):
    native_rows=[[1,0,0,10.0,1.0,0],[2,0,0,11.0,1.0,1],[3,0,3,12.0,1.0,1]]
    native=_mtz(tmp_path/'native.mtz',native_rows,extra=[('F','F'),('SIGF','Q'),('FreeR_flag','I')])
    staged=tmp_path/'work.mtz'
    info=trace_feedback.merge_work(native,staged)
    m=gemmi.read_mtz_file(str(staged))
    flags=np.asarray(m.column_with_label('FreeR_flag'))
    assert info['work_rows']==2 and info['excluded']==1
    assert 0 not in flags.tolist(),'free reflections must never enter density modification'


def test_empty_trace_stops_the_feedback_and_keeps_the_dm_map(tmp_path,monkeypatch):
    """A marginal map can trace nothing; that must not discard the DM map."""
    work=tmp_path/'work.mtz';work.write_text('x');prior=tmp_path/'prior.mtz';prior.write_text('x')
    sequence=tmp_path/'seq.fasta';sequence.write_text('>a\nAAAA\n')
    dm_models=[]
    monkeypatch.setattr(trace_feedback,'merge_work',lambda src,dest:Path(dest).write_text('x'))
    monkeypatch.setattr(trace_feedback,'phase_input',lambda w,d,dest:Path(dest).write_text('x'))
    monkeypatch.setattr(trace_feedback,'merge_phase_map',lambda w,d,dest:Path(dest).write_text('x') or {'work_rows':1})
    def fake_dm(work_mtz,prior_mtz,destination,solvent_fraction,**kw):
        Path(destination).write_text('x');dm_models.append(kw.get('model_pdb'))
        return {'dm_cycles':kw.get('cycles')}
    def fake_trace(work_mtz,phase_mtz,seq,runtime,output_dir,**kw):
        model=Path(output_dir)/'model.pdb';model.parent.mkdir(parents=True,exist_ok=True);model.write_text('END\n')
        return {'model':str(model),'residues':0,'fragments':0}
    monkeypatch.setattr(parrot_dm,'run',fake_dm)
    monkeypatch.setattr(trace_feedback,'trace',fake_trace)
    report=trace_feedback.run(work,prior,sequence,'runtime',tmp_path/'out',solvent_fraction=.5,feedback_cycles=3)
    assert len(report['history'])==1                 # no second cycle was attempted
    assert 'no residues' in report['history'][0]['feedback_stopped']
    assert dm_models==[None]                          # the empty model was never fed back
    assert Path(report['final_map']).is_file()
    assert (tmp_path/'out/feedback_report.json').is_file()


def test_build_cycle_with_no_atoms_stops_before_refmac(tmp_path,monkeypatch):
    """An empty Buccaneer model must be reported as "no model", not as a REFMAC crash."""
    native=_mtz(tmp_path/'native.mtz',[[1,0,0,10.0,1.0,0]],extra=[('FP','F'),('SIGFP','Q'),('FreeR_flag','I')])
    dm=_mtz(tmp_path/'dm.mtz',[[1,0,0,10.0,1.0,45.0,0.8]],extra=[('FP','F'),('SIGFP','Q'),
        ('PHIB','P'),('FOM','W')])
    sequence=tmp_path/'seq.fasta';sequence.write_text('>a\nAAAA\n')
    runtime=tmp_path/'runtime';runtime.mkdir()
    calls=[]
    class Done:returncode=0
    def fake_run(command,**kw):
        calls.append(command[0])
        if 'cbuccaneer' in str(command[0]):
            out=Path(kw['input'].split('pdbout ')[1].splitlines()[0].strip())
            out.write_text('REMARK nothing was built\nEND\n')
        return Done()
    monkeypatch.setattr(ccp4_pipeline,'_find',lambda name,rt:str(tmp_path/name))
    monkeypatch.setattr(ccp4_pipeline.subprocess,'run',fake_run)
    with pytest.raises(RuntimeError) as caught:
        ccp4_pipeline.run(native,dm,sequence,runtime,tmp_path/'out',cycles=3)
    assert '未建出任何残基' in str(caught.value)
    assert not any('refmac' in str(c) for c in calls)   # REFMAC never ran on an empty model
    report=json.loads((tmp_path/'out/hybrid_report.json').read_text())
    assert report['status']=='no_model' and report['residues_built']==0


def test_merge_phase_map_writes_pipeline_schema_and_skips_free(tmp_path):
    native_rows=[[1,0,0,10.0,1.0,0],[2,0,0,11.0,1.0,1]]
    native=_mtz(tmp_path/'native.mtz',native_rows,extra=[('FP','F'),('SIGFP','Q'),('FreeR_flag','I')])
    dm_rows=[[2,0,0,11.0,1.0]+_hl(3.0,45.0)]
    dm=_mtz(tmp_path/'dm.mtz',dm_rows,extra=[('FP','F'),('SIGFP','Q'),
        ('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')])
    out=tmp_path/'map_phases.mtz'
    info=trace_feedback.merge_phase_map(native,dm,out)
    m=gemmi.read_mtz_file(str(out));labels=m.column_labels()
    for label in ('F','SIGF','FWT','PHWT','PHIB','FOM','HLA','HLB','HLC','HLD'):
        assert label in labels
    assert info['work_rows']==1 and m.nreflections==1
    assert float(np.asarray(m.column_with_label('PHIB'))[0])==pytest.approx(45.0,abs=1e-3)


def test_phase_input_never_injects_free_reflections(tmp_path):
    """The tracer input is the intersection with the work-only map; a free row is
    dropped rather than given a fabricated experimental phase."""
    native_rows=[[1,0,0,10.0,1.0,0],[2,0,0,11.0,1.0,1]]
    native=_mtz(tmp_path/'native.mtz',native_rows,extra=[('F','F'),('SIGF','Q'),('FreeR_flag','I')])
    dm_rows=[[2,0,0]+_hl(2.5,120.0)]
    dm=_mtz(tmp_path/'dm.mtz',dm_rows,extra=[('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')])
    out=tmp_path/'trace_input.mtz'
    info=trace_feedback.phase_input(native,dm,out)
    m=gemmi.read_mtz_file(str(out))
    flags=np.asarray(m.column_with_label('FreeR_flag')).astype(int)
    phases=np.asarray(m.column_with_label('PHIB'))
    assert info['rows']==1 and m.nreflections==1
    assert flags.tolist()==[1],'only work reflections may reach the tracer'
    assert float(phases[0])==pytest.approx(120.0,abs=1e-3)


def test_ccp4_pipeline_input_masks_free_phases(tmp_path):
    native_rows=[[1,0,0,10.0,1.0,0],[2,0,0,11.0,1.0,1]]
    native=_mtz(tmp_path/'native.mtz',native_rows,extra=[('F','F'),('SIGF','Q'),('FreeR_flag','I')])
    dm_rows=[[1,0,0,45.0,0.4],[2,0,0,90.0,0.6]]
    dm=_mtz(tmp_path/'dm.mtz',dm_rows,extra=[('PHIB','P'),('FOM','W')])
    out=tmp_path/'build_input.mtz'
    ccp4_pipeline.input_mtz(native,dm,out,free_value=0)
    m=gemmi.read_mtz_file(str(out));labels=m.column_labels()
    flags=np.asarray(m.column_with_label('FreeR_flag')).astype(int)
    ph=np.asarray(m.column_with_label('PHIB'))
    fom=np.asarray(m.column_with_label('FOM'))
    free=flags==0
    assert m.nreflections==2 and free.sum()==1
    assert not np.isfinite(ph[free]).all() and not np.isfinite(fom[free]).all()
    assert np.isfinite(ph[~free]).all() and np.isfinite(fom[~free]).all()
    assert np.allclose(np.asarray(m.column_with_label('FP'))[[0,1]],[10.0,11.0])


def test_work_from_experimental_supplies_the_free_set_a_phases_only_run_lacks(tmp_path):
    h=gemmi.make_miller_array(gemmi.UnitCell(40,40,40,90,90,90),gemmi.SpaceGroup('P 1'),3.0)
    flags=np.where(np.arange(len(h))%10==0,7,3)
    sad=_mtz(tmp_path/'sad.mtz',np.column_stack([h,np.full(len(h),30.),np.ones(len(h)),
                                                 np.zeros(len(h)),np.ones(len(h)),flags]),
             extra=[('FP','F'),('SIGFP','Q'),('DANO','D'),('SIGDANO','Q'),('PIPELINE_FREE','I')])
    experimental=_mtz(tmp_path/'experimental.mtz',
                      np.column_stack([h,np.full(len(h),30.),np.ones(len(h)),np.full(len(h),60.),
                                       np.full(len(h),.5),np.zeros((len(h),4))]),
                      extra=[('FP','F'),('SIGFP','Q'),('PHIB','P'),('FOM','W'),
                             ('HLA','A'),('HLB','A'),('HLC','A'),('HLD','A')])
    out=trace_feedback.work_from_experimental(experimental,sad,'PIPELINE_FREE',7,tmp_path/'work.mtz')
    merged=gemmi.read_mtz_file(out)
    assert merged.column_labels()==['H','K','L','FP','SIGFP','FreeR_flag']
    kept=np.array(merged)[:,5]!=7
    assert kept.sum()==int((flags!=7).sum())
