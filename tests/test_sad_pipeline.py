"""Synthetic data handoff, provenance, restart and actual OASIS stage contracts."""
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace
import gemmi
import numpy as np
import pytest
from phasebuild import sad_bridge as bridge
from phasebuild import sad_pipeline as pipeline
from phasebuild.data import fingerprint
from phasebuild.geometry import template
from phasebuild.trace import write_backbone


def mtz(path, names, rows):
    m=gemmi.Mtz(with_base=True);m.spacegroup=gemmi.SpaceGroup('P 1');m.set_cell_for_all(gemmi.UnitCell(16,16,16,90,90,90));m.add_dataset('synthetic')
    for name,kind in names:m.add_column(name,kind)
    m.set_data(np.array(rows,dtype=np.float32));m.update_reso();m.write_to_file(str(path));return m


@pytest.fixture
def inputs(tmp_path):
    cell=gemmi.UnitCell(16,16,16,90,90,90);sg=gemmi.SpaceGroup('P 1')
    h=gemmi.make_miller_array(cell,sg,3.)
    f=40+5*np.cos(np.arange(len(h)));sig=np.ones(len(h));flags=np.where(np.arange(len(h))%10==0,7,3)
    native=tmp_path/'native.mtz';sad=tmp_path/'sad.mtz'
    mtz(native,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],np.column_stack([h,f,sig,flags]))
    mtz(sad,[('FP','F'),('SIGFP','Q'),('DANO','D'),('SIGDANO','Q'),('FreeR_flag','I')],np.column_stack([h,f,sig,np.sin(np.arange(len(h))),sig,flags]))
    seq=tmp_path/'sequence.fasta';seq.write_text('>chainA\n'+'A'*12+'\n')
    heavy=tmp_path/'heavy.pdb';heavy.write_text('CRYST1   16.000   16.000   16.000  90.00  90.00  90.00 P 1           1\nHETATM    1  S   SUL A   1       1.560   3.240   4.680  1.00 10.00           S  \nEND\n')
    ref=tmp_path/'reference.mtz';mtz(ref,[('FC_BEST.F_phi.F','F'),('FC_BEST.F_phi.phi','P')],np.column_stack([h,f,(h@np.array([.13,.27,.39]))*360]))
    runtime=tmp_path/'runtime';(runtime/'bin').mkdir(parents=True)
    for name in ['cbuccaneer','refmac5']:
        p=runtime/'bin'/name;p.write_text('#!/bin/sh\nexit 0\n');p.chmod(0o755)
    lib=runtime/'lib/data/reference_structures';lib.mkdir(parents=True)
    shutil.copyfile(ref,lib/'reference-1tqw.mtz');shutil.copyfile(heavy,lib/'reference-1tqw.pdb')
    return dict(sad=str(sad),native=str(native),sequence=str(seq),heavy=str(heavy),reference=str(ref),
                runtime=str(runtime),fpp='S=0.56',output=str(tmp_path/'output'),free_value=7,copies=2,
                phase_cycles=1,dm_cycles=1,rounds=1,cycles=1)


def write_phases(path, native_path, free=7, dm=False):
    native=gemmi.read_mtz_file(native_path);a=np.array(native,copy=False);work=a[:,5]!=free;h=a[work,:3]
    rows=np.column_stack([h,a[work,3],a[work,4],np.full(len(h),60),np.full(len(h),.5),
                          np.full(len(h),99 if dm else 2),np.zeros((len(h),3))])
    mtz(path,[('FP','F'),('SIGFP','Q'),('PHIB','P'),('FOM','W')]+[(n,'A') for n in ['HLA','HLB','HLC','HLD']],rows)
    path.with_suffix('.json').write_text(json.dumps({'synthetic':True,'post_dm':dm}))


class FakeStages:
    def __init__(self, cfg, fail_dm=False, real_phase=False):self.cfg=cfg;self.calls=[];self.fail_dm=fail_dm;self.real_phase=real_phase
    def __call__(self,command,log,env):
        kind='build' if command[2].endswith('run_build_from_phases.py') else command[3]
        self.calls.append(kind)
        if kind=='dm' and self.fail_dm:
            self.fail_dm=False;raise RuntimeError('synthetic interrupted DM')
        if self.real_phase and kind!='build':return pipeline.run_command(command,log,env)
        out=Path(command[command.index('--output')+1])
        if kind in ['phase','dm']:write_phases(out,self.cfg['native'],dm=kind=='dm');return
        out.mkdir();model=out/'primary_model.cif'
        write_backbone(model,[template(-63,-42,length=3)+[5,5,5]],SimpleNamespace(cell=gemmi.UnitCell(16,16,16,90,90,90),sg=gemmi.SpaceGroup('P 1')))
        (out/'hybrid_report.json').write_text(json.dumps({'engine':'phasebuild-integrated-external-buccaneer-refmac',
            'status':'completed','primary_model':fingerprint(model),'stages':[{'round':1}],
            'primary_metric':{'r_work':.2,'r_free':.3}}))


def test_bridge_preserves_native_and_uses_fom_once(tmp_path,inputs):
    phase=tmp_path/'raw.mtz';dm=tmp_path/'dm.mtz';write_phases(phase,inputs['native']);write_phases(dm,inputs['native'],dm=True)
    before=Path(inputs['native']).read_bytes();target=tmp_path/'building.mtz'
    report=bridge.adapt_maps(inputs['native'],phase,dm,target,7)
    result=gemmi.read_mtz_file(str(target));original=gemmi.read_mtz_file(inputs['native'])
    np.testing.assert_array_equal(np.array(result)[:,:6],np.array(original))
    work=np.asarray(result.column_with_label('FreeR_flag'))!=7
    np.testing.assert_allclose(np.asarray(result.column_with_label('FWT'))[work],np.asarray(original.column_with_label('F'))[work]*.5)
    assert np.isnan(np.asarray(result.column_with_label('FWT'))[~work]).all()
    assert Path(inputs['native']).read_bytes()==before and report['experimental_prior']['sha256']==fingerprint(phase)['sha256']


def test_conflicting_original_flags_fail_before_output(tmp_path,inputs):
    p=Path(inputs['sad']);m=gemmi.read_mtz_file(str(p));a=np.array(m,copy=True);a[0,-1]=3;m.set_data(a);m.write_to_file(str(p))
    target=tmp_path/'prepared.mtz'
    with pytest.raises(ValueError,match='flags disagree'):bridge.prepare_sad(p,inputs['native'],target,7)
    assert not target.exists()


def test_malformed_dm_fom_and_wrong_cell_are_rejected(tmp_path,inputs):
    phase=tmp_path/'raw.mtz';dm=tmp_path/'dm.mtz';write_phases(phase,inputs['native']);write_phases(dm,inputs['native'],dm=True)
    m=gemmi.read_mtz_file(str(dm));a=np.array(m,copy=True);a[0,m.column_labels().index('FOM')]=1.2;m.set_data(a);m.write_to_file(str(dm))
    with pytest.raises(ValueError,match='FOM'):bridge.adapt_maps(inputs['native'],phase,dm,tmp_path/'none.mtz',7)
    m.set_cell_for_all(gemmi.UnitCell(13,12,12,90,90,90));m.write_to_file(str(dm))
    with pytest.raises(ValueError,match='Crystal'):bridge.adapt_maps(inputs['native'],phase,dm,tmp_path/'none.mtz',7)


@pytest.mark.skipif(sys.version_info<(3,12),reason='Unified engine uses Python 3.12')
def test_pipeline_restarts_from_verified_checkpoint_without_rephasing(inputs):
    runner=FakeStages(inputs,fail_dm=True)
    with pytest.raises(RuntimeError,match='interrupted'):pipeline.run_pipeline(inputs,runner=runner,progress=lambda x:None)
    root=Path(inputs['output']);raw=(root/'02_sad/experimental.mtz').read_bytes()
    report=pipeline.run_pipeline(inputs,resume=True,runner=runner,progress=lambda x:None)
    assert report['status']=='completed' and runner.calls==['phase','dm','dm','build']
    assert (root/'02_sad/experimental.mtz').read_bytes()==raw
    assert list(root.glob('03_dm.interrupted-*'))
    assert len((root/'01_inputs/building.fasta').read_text().split('>'))-1==2
    assert report['experimental_prior']['sha256']!=report['dm_phases']['sha256']


@pytest.mark.skipif(sys.version_info<(3,12),reason='Unified engine uses Python 3.12')
def test_checkpoint_tampering_and_protocol_changes_refuse_resume(inputs):
    runner=FakeStages(inputs);pipeline.run_pipeline(inputs,runner=runner,progress=lambda x:None)
    with pytest.raises(ValueError,match='Resume refused'):pipeline.run_pipeline({**inputs,'phase_cycles':2},resume=True,runner=runner)
    p=Path(inputs['output'])/'02_sad/experimental.mtz';p.write_bytes(p.read_bytes()+b'tamper')
    with pytest.raises(ValueError,match='artifact changed'):pipeline.run_pipeline(inputs,resume=True,runner=runner)


def test_output_lock_rejects_a_second_runner(tmp_path):
    with pipeline.output_lock(tmp_path):
        with pytest.raises(ValueError,match='already running'):
            with pipeline.output_lock(tmp_path):pass


@pytest.mark.skipif(sys.version_info<(3,12),reason='Unified engine uses Python 3.12')
def test_gui_full_plan_is_read_only_and_ignores_hidden_map_and_seed(inputs):
    from phasebuild.gui_workflow import make_plan
    root=Path(__file__).resolve().parents[1]
    cfg,steps=make_plan(dict(mode='sad_to_model',sad_input=inputs['sad'],native=inputs['native'],sequence=inputs['sequence'],
        heavy=inputs['heavy'],fpp=inputs['fpp'],dm_reference=inputs['reference'],ccp4=inputs['runtime'],copies=2,
        output=inputs['output'],free_value=7,map='/missing/hidden-map',model='/missing/hidden-model',phase_prior='/missing/hidden-prior'),root)
    assert not Path(inputs['output']).exists()
    command=steps[0]['command'];assert 'run_sad_to_model.py' in command[2]
    assert '--model' not in command and command[command.index('--free-value')+1]=='7'
    assert command[command.index('--copies')+1]=='2' and cfg['model']==''


@pytest.mark.skipif(sys.version_info<(3,12),reason='Unified engine uses Python 3.12')
def test_real_oasis_sad_and_dm_feed_a_stub_construction_stage(inputs):
    runner=FakeStages(inputs,real_phase=True)
    report=pipeline.run_pipeline(inputs,runner=runner,progress=lambda x:None)
    assert report['status']=='completed' and runner.calls==['phase','dm','build']
    root=Path(inputs['output']);prepared=gemmi.read_mtz_file(str(root/'01_inputs/sad.mtz'))
    np.testing.assert_array_equal(prepared.column_with_label(bridge.FLAG).array,
                                 gemmi.read_mtz_file(inputs['native']).column_with_label('FreeR_flag').array)
    raw=gemmi.read_mtz_file(str(root/'02_sad/experimental.mtz'));final=gemmi.read_mtz_file(str(root/'03_dm/map_phases.mtz'))
    assert raw.nreflections==final.nreflections and raw.nreflections<prepared.nreflections


@pytest.mark.skipif(sys.version_info<(3,12),reason='Unified engine uses Python 3.12')
def test_resume_rejects_runtime_resource_changes(inputs):
    runner=FakeStages(inputs);pipeline.run_pipeline(inputs,runner=runner,progress=lambda x:None)
    reference=Path(inputs['runtime'])/'lib/data/reference_structures/reference-1tqw.pdb'
    reference.write_text(reference.read_text()+'REMARK changed runtime\n')
    with pytest.raises(ValueError,match='Resume refused'):
        pipeline.run_pipeline(inputs,resume=True,runner=runner)


def test_missing_hkl_stays_nan_and_existing_output_is_preserved(tmp_path,inputs):
    phase=tmp_path/'raw.mtz';dm=tmp_path/'dm.mtz';write_phases(phase,inputs['native']);write_phases(dm,inputs['native'],dm=True)
    m=gemmi.read_mtz_file(str(dm));a=np.array(m,copy=True);missing=tuple(a[0,:3].astype(int));m.set_data(a[1:]);m.write_to_file(str(dm))
    out=tmp_path/'map.mtz';bridge.adapt_maps(inputs['native'],phase,dm,out,7)
    result=gemmi.read_mtz_file(str(out));values=np.array(result,copy=False)
    row=next(row for row in values if tuple(row[:3].astype(int))==missing)
    assert np.isnan(row[result.column_labels().index('FWT')])
    before=out.read_bytes()
    with pytest.raises(FileExistsError):bridge.adapt_maps(inputs['native'],phase,dm,out,7)
    assert out.read_bytes()==before


def test_origin_choice_changes_are_not_silently_accepted():
    first=SimpleNamespace(cell=gemmi.UnitCell(20,20,20,90,90,90),spacegroup=gemmi.SpaceGroup('P n n n:1'))
    second=SimpleNamespace(cell=first.cell,spacegroup=gemmi.SpaceGroup('P n n n:2'))
    assert first.spacegroup.hm==second.spacegroup.hm
    with pytest.raises(ValueError,match='setting differs'):bridge.crystal(first,second)


def test_command_cancellation_cleans_up_its_independent_process_group(tmp_path):
    import os,signal,subprocess,time
    if os.name!='posix':pytest.skip('POSIX process groups')
    root=Path(__file__).resolve().parents[1];pidfile=tmp_path/'child.pid';driver=tmp_path/'driver.py'
    child="import os,time; from pathlib import Path; Path("+repr(str(pidfile))+").write_text(str(os.getpid())); time.sleep(120)"
    driver.write_text('import sys,os,signal\nfrom pathlib import Path\nsys.path.insert(0,'+repr(str(root/'src'))+')\nfrom phasebuild.sad_pipeline import run_command\ndef stop(s,f):raise InterruptedError()\nsignal.signal(signal.SIGTERM,stop)\ntry:\n with open('+repr(str(tmp_path/'child.log'))+',"w") as log:run_command([sys.executable,"-c",'+repr(child)+'],log,dict(os.environ))\nexcept InterruptedError:pass\n')
    process=subprocess.Popen([sys.executable,str(driver)])
    try:
        deadline=time.time()+15
        while not pidfile.exists() and time.time()<deadline:
            if process.poll() is not None:pytest.fail('Driver exited before child start')
            time.sleep(.03)
        assert pidfile.exists()
        pid=int(pidfile.read_text());process.terminate();process.wait(timeout=8)
        with pytest.raises(ProcessLookupError):os.kill(pid,0)
    finally:
        if process.poll() is None:process.kill();process.wait(timeout=3)
