"""Native-engine boundary contracts, especially HKL joins and honest priors."""
from pathlib import Path
import gemmi,numpy as np,pytest
from test_sad_pipeline import inputs,mtz
from test_trust_workflow import trusted
from phasebuild import bundled_pplus as b,sad_pipeline as pipeline
from phasebuild.sad_bridge import prepare_sad,FLAG
from phasebuild.sad_input import read_sad_pairs

@pytest.fixture
def bundled(trusted,tmp_path):
 cfg={**trusted,'method':'pplus_bundled','observation_prior':'','wavelength':.9792}
 runtime=Path(cfg['runtime'])
 # P+ needs exactly these three CCP4 programs; density modification is in-tree.
 for name in ['bin/phaser','bin/cbuccaneer','bin/refmac5']:
  p=runtime/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('#!/bin/sh\nexit 0\n');p.chmod(0o755)
 return cfg

def test_bundle_plan_and_explicit_experimental_wavelength(bundled):
 # The bundled profile picks a cycle budget (12 for oasis-dm v2), but an explicit
 # request must survive (a silently ignored CLI flag is a defect).
 cfg=pipeline.validate({k:v for k,v in bundled.items() if k!='dm_cycles'});plan=pipeline.stage_plan(cfg)
 assert cfg['phenix_bin']=='' and cfg['dm_cycles']==12
 assert pipeline.validate({**bundled,'dm_cycles':20})['dm_cycles']==20
 assert 'frontend.json' in '\n'.join(plan[1]['outputs'])
 assert plan[4]['command'][plan[4]['command'].index('--phase-prior')+1].endswith('02_sad/observations.mtz')
 for wave in [None,True,0,float('nan'),11]:
  with pytest.raises(ValueError,match='波长'):pipeline.validate({**bundled,'wavelength':wave})
 assert not Path(cfg['output']).exists()

def test_phase_only_keeps_bundled_runtime_and_fails_closed(bundled):
 cfg=pipeline.validate({**bundled,'build_model':False})
 assert cfg['runtime']==bundled['runtime'] and len(pipeline.stage_plan(cfg))==3
 (Path(cfg['runtime'])/'bin/phaser').unlink()
 with pytest.raises(ValueError,match='引擎缺失'):pipeline.validate(bundled)

def test_work_native_amplitudes_join_by_key_and_free_excluded(bundled,tmp_path):
 cfg=pipeline.validate(bundled);prepared=tmp_path/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
 native=gemmi.read_mtz_file(cfg['native']);original=np.array(native);native.set_data(original[::-1]);native.write_to_file(cfg['native'])
 out=tmp_path/'work.mtz';b.stage_work(cfg,out);m=gemmi.read_mtz_file(str(out));lookup={tuple(a[:3]):a for a in original}
 for a in np.array(m):
  expected=lookup[tuple(a[:3])];assert a[-1]!=cfg['free_value'];np.testing.assert_array_equal(a[-3:],expected[3:6])
 assert len(m.array)==sum(original[:,5]!=cfg['free_value'])
 sad=gemmi.read_mtz_file(str(prepared));sa=np.array(sad);fl=sad.column_labels().index(FLAG)
 np.testing.assert_array_equal(m.make_miller_array(),sa[sa[:,fl]!=cfg['free_value'],:3])

def test_prior_missing_hl_uniform_and_moments_consistent(bundled,tmp_path):
 cfg=pipeline.validate(bundled);prepared=tmp_path/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value'])
 data,_,_=read_sad_pairs(prepared,free_label=FLAG,free_value=cfg['free_value']);h=data.hkl
 hl=np.tile([2.,-.5,.3,.1],(len(h),1));hl[0]=np.nan
 raw=tmp_path/'raw.mtz';mtz(raw,[(x,'A') for x in ['HLA','HLB','HLC','HLD']],np.column_stack([h,hl])[::-1])
 out=tmp_path/'prior.mtz';meta=b.extract_prior(raw,data,out);assert meta['uninformative_missing_HL_rows']==1
 m=gemmi.read_mtz_file(str(out));lookup={tuple(a[:3]):a for a in np.array(m)}
 row=lookup[tuple(h[0])];np.testing.assert_array_equal(row[-4:],0)
 from oasis_py.sad import posterior_moments
 phase,fom=posterior_moments(data,np.nan_to_num(hl))
 for i,key in enumerate(h):
  row=lookup[tuple(key)];assert row[6]==pytest.approx(fom[i],abs=1e-6)


def test_phase_confidence_correction_reaches_build_command_only(bundled):
 cfg=pipeline.validate({**bundled,'refinement_prior':'dm_phifom','phase_blur_scale':.7,'phase_blur_b':20.})
 stages=pipeline.stage_plan(cfg);build=next(row['command'] for row in stages if row['id']=='build')
 assert build[build.index('--phase-blur-scale')+1]=='0.7'
 assert build[build.index('--phase-blur-b')+1]=='20.0'
 for row in stages:
  if row['id']!='build':assert '--phase-blur-scale' not in (row.get('command') or [])
 assert not Path(cfg['output']).exists()

@pytest.mark.parametrize('update',[{'phase_blur_b':None},{'phase_blur_scale':None},{'phase_blur_scale':False},{'refinement_prior':'original_sad'},{'method':'likelihood'},{'build_model':False}])
def test_incompatible_or_incomplete_phase_confidence_rejected(bundled,update):
 with pytest.raises(ValueError,match='[Pp]hase'):
  pipeline.validate({**bundled,'refinement_prior':'dm_phifom','phase_blur_scale':.7,'phase_blur_b':20.,**update})
 assert not Path(bundled['output']).exists()


def test_per_round_cleanup_reaches_actual_build_command(bundled):
 cfg=pipeline.validate({**bundled,'copies':1,'refinement_prior':'dm_phifom','per_round_cleanup':True})
 stages=pipeline.stage_plan(cfg);build=next(r['command'] for r in stages if r['id']=='build')
 assert '--per-round-cleanup' in build
 assert build[build.index('--cleanup-map')+1].endswith('/03_dm/map_phases.mtz')
 assert build[build.index('--asu-copies')+1]=='1'
 assert not Path(cfg['output']).exists()

@pytest.mark.parametrize('update',[{'copies':2},{'build_model':False},{'refinement_prior':'original_sad'},{'method':'pplus_internal','refinement_prior':'original_sad'},{'per_round_cleanup':'true'}])
def test_incompatible_per_round_cleanup_fails_before_output(bundled,update):
 with pytest.raises(ValueError,match='[Pp]er-round'):
  pipeline.validate({**bundled,'copies':1,'refinement_prior':'dm_phifom','per_round_cleanup':True,**update})
 assert not Path(bundled['output']).exists()


def test_per_round_cleanup_rejects_multiple_sequences_at_preflight(bundled):
 Path(bundled['sequence']).write_text('>one\nAAAAAAAAAAAA\n>two\nCCCCCCCCCCCC\n')
 with pytest.raises(ValueError,match='one sequence record'):
  pipeline.validate({**bundled,'copies':1,'refinement_prior':'dm_phifom','per_round_cleanup':True})
 assert not Path(bundled['output']).exists()
