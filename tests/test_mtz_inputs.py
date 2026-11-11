"""Selection, original flags and single-file input regression coverage."""
import json
from pathlib import Path
import gemmi
import numpy as np
import pytest
from phasebuild import mtz_inputs as mod
from phasebuild.refinement import reflection_partition
from test_sad_pipeline import mtz


def config(tmp_path,free=True,mean=True):
    h=gemmi.make_miller_array(gemmi.UnitCell(16,16,16,90,90,90),gemmi.SpaceGroup('P 1'),3.)
    columns=[('F(+)','G'),('SIGF(+)','L'),('F(-)','G'),('SIGF(-)','L')]
    arrays=[np.full(len(h),20.),np.full(len(h),2.),np.full(len(h),16.),np.full(len(h),3.)]
    if mean:columns += [('F','F'),('SIGF','Q')];arrays += [np.full(len(h),19.),np.full(len(h),1.5)]
    if free:columns += [('R-free-flags','I')];arrays += [np.where(np.arange(len(h))%10==0,7,3)]
    sad=tmp_path/'sad.mtz';mtz(sad,columns,np.column_stack([h,*arrays]))
    return dict(sad=str(sad),native='',labels='',native_labels='',free_value=7,seed=42)


def test_single_file_preserves_selected_flags_and_mean(tmp_path):
    cfg=config(tmp_path);before=Path(cfg['sad']).read_bytes();s,n,a=mod.build(cfg)
    np.testing.assert_array_equal(n.column_with_label('F'),19.)
    original=gemmi.read_mtz_file(cfg['sad']);np.testing.assert_array_equal(n.column_with_label('FreeR_flag'),original.column_with_label('R-free-flags'))
    assert a['native_omitted'] and not a['flags_generated']
    assert a['free_column']['label']=='R-free-flags'
    assert Path(cfg['sad']).read_bytes()==before


def test_derive_means_and_singletons_no_false_bijvoet_pair(tmp_path):
    cfg=config(tmp_path,mean=False);m=gemmi.read_mtz_file(cfg['sad']);v=np.array(m,copy=True)
    v[0,5:7]=np.nan;m.set_data(v);m.write_to_file(cfg['sad'])
    s,n,a=mod.build(cfg)
    assert n.column_with_label('F')[0]==20
    assert n.column_with_label('SIGF')[0]==2
    np.testing.assert_allclose(np.asarray(n.column_with_label('F'))[1:],18.)
    np.testing.assert_allclose(np.asarray(n.column_with_label('SIGF'))[1:],np.hypot(2,3)/2)
    assert np.isnan(s.column_with_label('F(-)')[0])
    assert a['mean_source'].startswith('derived')


def test_generated_flags_deterministic_by_hkl_and_reported(tmp_path):
    cfg=config(tmp_path,free=False);_,n,a=mod.build(cfg)
    assert a['flags_generated'] and a['seed']==42
    expected=dict(zip(map(tuple,n.make_miller_array()),np.asarray(n.column_with_label('FreeR_flag'))))
    m=gemmi.read_mtz_file(cfg['sad']);m.set_data(np.array(m,copy=True)[::-1]);m.write_to_file(cfg['sad'])
    _,n2,a2=mod.build(cfg)
    assert expected==dict(zip(map(tuple,n2.make_miller_array()),np.asarray(n2.column_with_label('FreeR_flag'))))
    out=tmp_path/'native.mtz';n2.write_to_file(str(out));_,part=reflection_partition(out,free_value=7)
    assert 0<part['usable_free']<part['usable_work']


def test_nonstandard_labels_manual_indices_and_duplicate_label_selection(tmp_path):
    cfg=config(tmp_path);m=gemmi.read_mtz_file(cfg['sad'])
    for c in m.columns[3:7]:c.label='duplicate'
    m.write_to_file(cfg['sad']);cfg['labels']='fplus=@3,sigfplus=@4,fminus=@5,sigfminus=@6'
    s,n,a=mod.build(cfg);assert s.column_labels()[3:]==list(mod.SAD_NAMES)
    assert a['sad_columns']['fplus']['index']==3
    with pytest.raises(ValueError,match='唯一'):
        cfg['labels']='fplus=duplicate,sigfplus=@4,fminus=@5,sigfminus=@6';mod.build(cfg)


def test_multiple_means_require_explicit_choice(tmp_path):
    cfg=config(tmp_path);m=gemmi.read_mtz_file(cfg['sad']);arr=np.array(m,copy=True)
    m.add_column('FP','F');m.add_column('SIGFP','Q');m.set_data(np.column_stack([arr,np.full(len(arr),80),np.ones(len(arr))]));m.write_to_file(cfg['sad'])
    with pytest.raises(ValueError,match='多组'):mod.build(cfg)
    cfg['native_labels']='f=FP,sigf=SIGFP';_,n,a=mod.build(cfg);np.testing.assert_array_equal(n.column_with_label('F'),80.)
    cfg['native_labels']='f=__derive__,sigf=__derive__';_,n,a=mod.build(cfg);np.testing.assert_array_equal(n.column_with_label('F'),18.)


def test_no_flags_in_mean_file_preserves_existing_sad_flags_by_hkl(tmp_path):
    cfg=config(tmp_path);sad=gemmi.read_mtz_file(cfg['sad']);raw=np.array(sad,copy=True)
    obs=tmp_path/'obs.mtz';mtz(obs,[('FOBS','F'),('SIGFOBS','Q')],raw[:,[0,1,2,7,8]][::-1])
    cfg['native']=str(obs);_,n,a=mod.build(cfg)
    np.testing.assert_array_equal(n.column_with_label('FreeR_flag'),raw[::-1,9]);assert not a['flags_generated']


def test_preparation_and_readonly_validation_canonicalize_without_changes(tmp_path):
    cfg=config(tmp_path);before=Path(cfg['sad']).read_bytes()
    with mod.checked_config(cfg) as (checked,a):
        assert Path(checked['native']).is_file() and not checked['input_columns']
        assert checked['labels']==''
    assert not Path(checked['native']).exists()
    folder=tmp_path/'prepare';folder.mkdir();report=mod.prepare(cfg,folder)
    assert not report['flags_regenerated'] and (folder/'columns.json').is_file()
    assert Path(cfg['sad']).read_bytes()==before
    prepared=gemmi.read_mtz_file(str(folder/'sad.mtz'));assert 'PIPELINE_FREE' in prepared.column_labels()


def test_intensities_are_not_silently_used_as_amplitudes(tmp_path):
    cfg=config(tmp_path);m=gemmi.read_mtz_file(cfg['sad']);m.columns[3].type='K';m.write_to_file(cfg['sad'])
    with pytest.raises(ValueError,match='振幅'):mod.build(cfg)


def test_free_conflict_never_overwritten(tmp_path):
    cfg=config(tmp_path);m=gemmi.read_mtz_file(cfg['sad']);v=np.array(m,copy=True);v[:,9]=np.where(v[:,9]==7,3,7);m.set_data(v)
    other=tmp_path/'native.mtz';m.write_to_file(str(other));cfg['native']=str(other)
    with pytest.raises(ValueError,match='不一致'):mod.build(cfg)


def test_known_existing_free_set_cannot_be_regenerated(tmp_path):
    cfg=config(tmp_path);cfg['native_labels']='free=__generate__'
    with pytest.raises(ValueError,match='禁止重新生成'):mod.build(cfg)


def test_same_column_cannot_be_used_for_both_bijvoet_sides(tmp_path):
    cfg=config(tmp_path);cfg['labels']='fplus=@3,sigfplus=@4,fminus=@3,sigfminus=@6'
    with pytest.raises(ValueError,match='不同的列'):mod.build(cfg)


def test_selected_views_preserve_dataset_metadata_without_guessing_sad_wavelength(tmp_path):
    cfg=config(tmp_path);m=gemmi.read_mtz_file(cfg['sad']);ds=m.datasets[1]
    ds.wavelength=1.234;ds.project_name='project';ds.crystal_name='crystal';ds.dataset_name='observed';m.title='original observations';m.write_to_file(cfg['sad'])
    s,n,a=mod.build(cfg)
    for view in (s,n):
        out=view.datasets[1]
        assert out.wavelength==pytest.approx(1.234) and out.dataset_name=='observed'
        assert out.project_name=='project' and out.crystal_name=='crystal'
        assert view.title=='original observations'
    assert 'wavelength' not in cfg
