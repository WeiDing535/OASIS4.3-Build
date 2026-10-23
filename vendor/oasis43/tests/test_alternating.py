import numpy as np
import gemmi
import pytest
from test_science import dataset
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py.native_density import density_modify_native
from oasis_py.alternating import density_modify_alternating


def setup(tmp_path):
    data=dataset('P 21 21 21',extent=4);rng=np.random.default_rng(902)
    data.f=rng.uniform(30,90,len(data.f))
    def initial():
        hl=rng.normal(size=(len(data.f),4))*.35
        p,f=posterior_moments(data,hl)
        return PhaseResult(p,f,None,hl,{'method':'sad_likelihood','anomalous_scale_fit_version':'wide_grid_null_v2'})
    a,b=initial(),initial();ref=tmp_path/'generic.mtz'
    m=gemmi.Mtz(with_base=True);m.cell=data.cell;m.spacegroup=data.spacegroup;m.add_dataset('reference')
    m.add_column('FC_BEST.F_phi.F','F');m.add_column('FC_BEST.F_phi.phi','P')
    m.set_data(np.column_stack([data.hkl,data.f,np.rad2deg(a.phase)]).astype(np.float32));m.write_to_file(str(ref))
    return data,a,b,ref


@pytest.mark.parametrize('pattern',['A','B'])
def test_single_method_exactly_matches_existing_engine(tmp_path,pattern):
    d,a,b,ref=setup(tmp_path)
    old=density_modify_native(d,a if pattern=='A' else b,.5,ref,cycles=3)
    new=density_modify_alternating(d,a,b,.5,ref,cycles=3,pattern=pattern)
    np.testing.assert_array_equal(old.hl,new.hl)
    assert [h['active_prior'] for h in new.metadata['cycles']]==[pattern]*3


def test_identical_priors_do_not_gain_confidence_by_alternating(tmp_path):
    d,a,b,ref=setup(tmp_path)
    old=density_modify_native(d,a,.5,ref,cycles=3)
    new=density_modify_alternating(d,a,a,.5,ref,cycles=3,pattern='AB')
    np.testing.assert_array_equal(old.hl,new.hl)


def test_fixed_priors_replace_evidence_and_callback_cannot_change_state(tmp_path,monkeypatch):
    import oasis_py.native_density as native
    d,a,b,ref=setup(tmp_path);captured=[]
    def fit(data,fc,prior):
        captured.append(prior.copy());return np.full(len(data.f),.2+.1j),[]
    monkeypatch.setattr(native,'fit_shell_likelihood',fit)
    trace=[]
    def callback(cycle,p,f):trace.append(cycle);p[:]=0;f[:]=0
    new=density_modify_alternating(d,a,b,.5,ref,cycles=4,pattern='BA',cycle_callback=callback)
    for got,want in zip(captured,[b.hl,a.hl,b.hl,a.hl]):np.testing.assert_array_equal(got,want)
    expected=a.hl.copy();expected[:,0]+=.2;expected[:,1]+=.1
    np.testing.assert_array_equal(new.hl,expected)
    assert trace==[1,2,3,4] and new.metadata['final_active_prior']=='A'
    assert new.metadata['prior_repeatedly_accumulated'] is False
    p,f=posterior_moments(d,new.hl)
    np.testing.assert_array_equal(p,new.phase);np.testing.assert_array_equal(f,new.fom)
    with pytest.raises(ValueError,match='pattern'):
        density_modify_alternating(d,a,b,.5,ref,pattern='best')


@pytest.mark.parametrize('options,message',[
    (['--sad-method','alternating'],'requires --native-dm'),
    (['--dm-pattern','AB'],'requires --sad-method alternating'),
    (['--sad-method','alternating','--native-dm','--solvent-fraction','.5',
      '--reference-mtz','generic.mtz','--dm-anisotropy'],'no validated anisotropy'),
])
def test_cli_rejects_unsupported_combinations(tmp_path,monkeypatch,capsys,options,message):
    import oasis_py.cli as cli
    monkeypatch.setattr(cli,'load',lambda args:(dataset(),{}))
    assert cli.main(['phase','input.mtz','--out',str(tmp_path/'out.mtz'),*options])==2
    assert message in capsys.readouterr().err


def test_solve_sad_rejects_alternation_before_search_or_output(tmp_path,capsys):
    from oasis_py.cli import main
    folder=tmp_path/'output'
    assert main(['solve-sad','missing.mtz','--sad-method','alternating',
        '--sites','2','--fpp','Se=3.8','--out-dir',str(folder)])==2
    assert 'only via phase --native-dm' in capsys.readouterr().err
    assert not folder.exists()
