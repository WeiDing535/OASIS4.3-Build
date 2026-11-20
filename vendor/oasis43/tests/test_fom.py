import numpy as np
import pytest
from oasis_py.fom import FomCalibration,isotonic
from oasis_py.phasing import PhaseResult
from test_science import dataset


def test_weighted_pool_adjacent_violators():
    np.testing.assert_allclose(isotonic([.1,.8,.2,.9],[1,1,2,1]),[.1,.4,.4,.9])


def test_calibration_recovers_known_curve_on_independent_cases():
    rng=np.random.default_rng(83);records=[]
    for i in range(8):
        f=rng.uniform(.01,.99,2000);c=rng.random(len(f))<.2
        records.append({'id':str(i),'fom':f,'cosine':.6*f,'centric':c})
    model=FomCalibration.fit(records,'native_classical_dm_experimental')
    f=np.linspace(.05,.95,100);prediction=model.predict(f,np.zeros(100,bool))
    np.testing.assert_allclose(prediction,.6*f,atol=.01)
    assert np.all(np.diff(prediction)>=0)


def test_calibration_application_keeps_phase_and_refuses_training_target():
    d=dataset();f=np.linspace(.01,.99,len(d.f));c=d.centric
    m=FomCalibration.fit([{'id':'train','fom':f,'cosine':f*.5,'centric':c}],'test')
    result=PhaseResult(f.copy(),f,None,np.zeros((len(f),4)),{'method':'test'})
    with pytest.raises(ValueError,match='training set'):m.apply(d,result,'train')
    out=m.apply(d,result,'new_target')
    np.testing.assert_array_equal(out.phase,result.phase)
    assert out.hl is None and out.probability is None and np.max(out.fom)<.6


def test_model_cannot_be_applied_to_different_density_modification_configuration():
    d=dataset();f=np.linspace(.01,.99,len(d.f))
    model=FomCalibration.fit([{'id':'train','fom':f,'cosine':f*.5,'centric':d.centric}],'native_classical_dm_experimental')
    result=PhaseResult(f,f,None,None,{'method':'native_classical_dm_experimental','anisotropy':{'status':'fitted'}})
    with pytest.raises(ValueError,match='configuration mismatch'):model.apply(d,result,'new')


def test_old_fom_model_refuses_corrected_scale_phases_until_revalidated():
    d=dataset();f=np.linspace(.01,.99,len(d.f))
    model=FomCalibration.fit([{'id':'train','fom':f,'cosine':f*.5,'centric':d.centric}],'native_classical_dm_experimental')
    result=PhaseResult(f,f,None,None,{'method':'native_classical_dm_experimental','initial_scale_fit_version':'wide_grid_null_v2'})
    with pytest.raises(ValueError,match='scale-fit versions differ'):model.apply(d,result,'new')


def test_version_matched_retrained_model_accepts_corrected_scale_results():
    d=dataset();f=np.linspace(.01,.99,len(d.f))
    model=FomCalibration.fit([{'id':'train','fom':f,'cosine':f*.5,'centric':d.centric}],'native_classical_dm_experimental',initial_scale_fit_version='wide_grid_null_v2')
    result=PhaseResult(f,f,None,None,{'method':'native_classical_dm_experimental','initial_method':'sad_likelihood','initial_scale_fit_version':'wide_grid_null_v2'})
    calibrated=model.apply(d,result,'new')
    np.testing.assert_array_equal(calibrated.phase,result.phase)
    assert calibrated.metadata['fom_calibration']['initial_scale_fit_version']=='wide_grid_null_v2'


def test_ensemble_calibration_cannot_be_applied_before_density_modification():
    d=dataset();f=np.linspace(.01,.99,len(d.f))
    model=FomCalibration.fit([{'id':'train','fom':f,'cosine':f*.5,'centric':d.centric}],'phase_mixture',initial_scale_fit_version='wide_grid_null_v2')
    model.model['ensemble_stage']='after_separate_native_dm'
    result=PhaseResult(f,f,None,None,{'method':'phase_mixture','initial_scale_fit_version':'wide_grid_null_v2','ensemble_stage':'initial_phases'})
    with pytest.raises(ValueError,match='ensemble stage mismatch'):model.apply(d,result,'new')


def test_numerically_cancelled_mixture_cannot_acquire_phase_confidence():
    f=np.linspace(.01,.99,100)
    model=FomCalibration.fit([{'id':'train','fom':f,'cosine':.5*f,'centric':np.zeros(100,bool)}],'test')
    prediction=model.predict(np.array([0.,1e-16,1e-10]),np.zeros(3,bool))
    np.testing.assert_array_equal(prediction,0)
