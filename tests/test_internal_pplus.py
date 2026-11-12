"""No-suite scientific route and explicit portability boundaries."""
import json
from pathlib import Path
import subprocess
import sys
import gemmi
import numpy as np
import pytest
from test_sad_pipeline import inputs,mtz
from test_trust_workflow import trusted
from phasebuild import sad_pipeline as pipeline,internal_pplus as internal,oasis_workflow
from phasebuild.sad_bridge import prepare_sad,FLAG
from phasebuild.sad_input import read_sad_pairs,read_phase_prior


@pytest.fixture
def own(trusted):
    cfg={**trusted,'method':'pplus_internal','observation_prior':'','phenix_bin':'/not-installed','reference':str(Path(trusted['sad']).parent/'reference.mtz')}
    m=gemmi.read_mtz_file(cfg['sad']);a=np.array(m);d=12*np.sin(np.arange(len(a)));f=(a[:,3]+a[:,5])/2;a[:,3]=f+d/2;a[:,5]=f-d/2;m.set_data(a);m.write_to_file(cfg['sad'])
    return cfg


def test_internal_preflight_needs_no_phenix_or_external_phase_file(own,monkeypatch):
    monkeypatch.setenv('OASIS_INSTALLATION_FREE','1')
    cfg=pipeline.validate(own)
    assert cfg['phenix_bin']=='' and cfg['observation_prior']=='' and cfg['phase_cycles']==100
    stages=pipeline.stage_plan(cfg)
    assert len(stages)==6 and '内部' in stages[1]['name'] and '内部' in stages[2]['name']
    checked,steps=oasis_workflow.make_plan(own,Path(__file__).resolve().parents[1]);assert checked['method']=='pplus_internal'
    assert not Path(own['output']).exists()
    with pytest.raises(ValueError,match='禁止调用 Phenix'):pipeline.validate({**own,'method':'pplus_trust'})
    with pytest.raises(ValueError,match='不接受外部'):pipeline.validate({**own,'observation_prior':'from_phaser.mtz'})


def test_installation_free_trusted_route_is_allowed_when_parrot_replaces_phenix(trusted,monkeypatch):
    """The ban exists for the RESOLVE/CCTBX backend; the Parrot route never calls Phenix."""
    monkeypatch.setenv('OASIS_INSTALLATION_FREE','1')
    cfg=pipeline.validate({**trusted,'dm_engine':'parrot','phenix_bin':'/not-installed'})
    assert cfg['phenix_bin']=='' and cfg['dm_engine']=='parrot'
    assert cfg['dm_cycles']==40 and cfg['phase_cycles']==100
    with pytest.raises(ValueError,match='禁止调用 Phenix'):
        pipeline.validate({**trusted,'dm_engine':None,'phenix_bin':'/not-installed'})


def test_math_phase_dm_execute_no_external_program_and_keep_original_prior(own,tmp_path,monkeypatch):
    cfg=pipeline.validate({**own,'dm_cycles':1});root=tmp_path/'math';root.mkdir();phase=root/'02_sad';dm=root/'03_dm';phase.mkdir();dm.mkdir()
    prepared=root/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
    def prohibited(*args,**kwargs):raise AssertionError('Phasing and DM must execute no external program')
    monkeypatch.setattr(subprocess,'Popen',prohibited)
    internal.phase(cfg,phase);original=(phase/'observations.mtz').read_bytes();internal.dm(cfg,dm)
    assert original==(phase/'observations.mtz').read_bytes()
    meta=json.loads((phase/'experimental.json').read_text());assert meta['Phenix_executed'] is False and meta['Phaser_input_used'] is False and meta['mean_abs_probability_change']>0
    result=gemmi.read_mtz_file(str(dm/'map_phases.mtz'))
    assert {'F','SIGF','FWT','PHWT','PHIB','FOM','HLA','HLB','HLC','HLD'}<=set(result.column_labels())
    assert np.isfinite(np.array(result)).all()
    flags={tuple(h):v for h,v in zip(gemmi.read_mtz_file(cfg['native']).make_miller_array(),gemmi.read_mtz_file(cfg['native']).column_with_label('FreeR_flag').array)}
    assert all(flags[tuple(h)]!=cfg['free_value'] for h in result.make_miller_array())
    assert not {'phaser','mmtbx','solve_resolve'} & set(sys.modules)


def test_dm_uses_separate_seed_but_does_not_accumulate_observations(own,tmp_path):
    cfg=pipeline.validate(own);prepared=tmp_path/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
    out=tmp_path/'phase';out.mkdir();internal.phase(cfg,out)
    data,_,_=read_sad_pairs(prepared,free_label=FLAG,free_value=cfg['free_value']);prior=read_phase_prior(out/'observations.mtz',data);before=prior.hl.copy()
    seed=gemmi.read_mtz_file(str(out/'experimental.mtz'));ph=np.deg2rad(seed.column_with_label('PHIB').array);fom=seed.column_with_label('FOM').array
    a,ca=internal.seeded_dm(data,prior,ph,fom,.45,cfg['reference'],cycles=1,random_seed=1)
    b,cb=internal.seeded_dm(data,prior,ph+.8*np.sin(np.arange(len(ph))),fom,.45,cfg['reference'],cycles=1,random_seed=1)
    np.testing.assert_array_equal(prior.hl,before);assert np.max(abs(ca-cb))>1e-4
    assert 'Original' in a.metadata['fixed_HL'] and 'never accumulated' in b.metadata['fixed_HL']


def _seed(own, tmp_path):
    cfg=pipeline.validate(own);prepared=tmp_path/'prepared.mtz';prepare_sad(cfg['sad'],cfg['native'],prepared,cfg['free_value']);cfg['prepared_sad']=str(prepared)
    out=tmp_path/'phase';out.mkdir();internal.phase(cfg,out)
    data,_,_=read_sad_pairs(prepared,free_label=FLAG,free_value=cfg['free_value']);prior=read_phase_prior(out/'observations.mtz',data)
    seed=gemmi.read_mtz_file(str(out/'experimental.mtz'));return cfg,data,prior,np.deg2rad(seed.column_with_label('PHIB').array),seed.column_with_label('FOM').array


def test_dm_ablation_switches_default_to_shipping_behaviour(monkeypatch):
    for name in ['OASIS_INTERNAL_DM_COMBINE','OASIS_INTERNAL_DM_GAMMA','OASIS_INTERNAL_DM_MASK',
                 'OASIS_INTERNAL_DM_WANG_RADIUS','OASIS_INTERNAL_DM_WANG_EDGE','OASIS_INTERNAL_DM_MASK_RADIUS',
                 'OASIS_INTERNAL_DM_LADDER','OASIS_INTERNAL_DM_MAP_WEIGHT','OASIS_INTERNAL_DM_MASK_NPY']:
        monkeypatch.delenv(name,raising=False)
    assert internal._dm_options()=={'combine':'mlhl','gamma_scale':1.0,'mask_kind':'variance',
        'wang_radius':None,'wang_edge':1.0,'mask_radius':4.0,'map_weight':1.0}


def test_dm_ladder_and_wang_mask_are_reproducible_and_keep_observation_prior(own,tmp_path,monkeypatch):
    cfg,data,prior,ph,fom=_seed(own,tmp_path)
    monkeypatch.setenv('OASIS_INTERNAL_DM_MASK','wang');monkeypatch.setenv('OASIS_INTERNAL_DM_WANG_RADIUS','3.0')
    monkeypatch.setenv('OASIS_INTERNAL_DM_LADDER','6.0,3.0');monkeypatch.setenv('OASIS_INTERNAL_DM_COMBINE','map')
    monkeypatch.setenv('OASIS_INTERNAL_DM_GAMMA','0')
    before=prior.hl.copy();trace=tmp_path/'trace'
    ladder=internal._ladder_from_env(data)
    a,ca=internal.seeded_dm(data,prior,ph,fom,.45,cfg['reference'],cycles=4,random_seed=5,trace_dir=trace,
        ladder=ladder,**internal._dm_options())
    b,cb=internal.seeded_dm(data,prior,ph,fom,.45,cfg['reference'],cycles=4,random_seed=5,
        ladder=ladder,**internal._dm_options())
    np.testing.assert_array_equal(prior.hl,before)
    np.testing.assert_allclose(ca,cb,atol=0,rtol=0)
    assert [row['resolution_cutoff'] for row in a.metadata['cycles']]==[6.0,3.0,None,None]
    assert np.isfinite(ca).all() and np.isfinite(a.fom).all() and (a.fom>=0).all() and (a.fom<=1).all()
    assert (trace/'cycle_04.mtz').is_file() and (trace/'mask_last.npy').is_file()


def test_invalid_dm_ablation_switches_are_rejected(own,tmp_path):
    cfg,data,prior,ph,fom=_seed(own,tmp_path)
    for kwargs,match in [({'combine':'automatic'},'combination'),({'mask_kind':'wang2'},'mask kind'),
                         ({'ladder':[0.0,None]},'resolution ladder'),({'ladder':[]},'resolution ladder')]:
        with pytest.raises(ValueError,match=match):
            internal.seeded_dm(data,prior,ph,fom,.45,cfg['reference'],cycles=1,random_seed=1,**kwargs)


def test_wang_mask_is_calibrated_to_the_requested_solvent_fraction():
    from oasis_py.native_density import metric_filter, wang_mask
    cell=gemmi.UnitCell(40,40,40,90,90,90);shape=(32,32,32)
    grid=gemmi.FloatGrid(np.zeros(shape,dtype=np.float32),cell,gemmi.SpaceGroup('P1'))
    rng=np.random.default_rng(3);values=np.zeros(shape)
    values[8:24,8:24,8:24]=1.0;values+=rng.normal(scale=.02,size=shape)
    grid.array[:]=values
    kernel=metric_filter(cell,shape,2.0)
    for solvent in (0.4,0.5,0.6):
        mask=wang_mask(grid,solvent,kernel,cell,radius=None)
        assert abs(float(mask.mean())-(1-solvent))<0.05
        assert float(mask.min())>=0 and float(mask.max())<=1.01
