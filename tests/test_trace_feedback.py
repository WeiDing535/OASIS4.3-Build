"""Unit tests for the work-only Parrot + trace-feedback density-modification layer.

These tests use synthetic MTZs and never execute Parrot, Buccaneer or REFMAC.
"""
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
