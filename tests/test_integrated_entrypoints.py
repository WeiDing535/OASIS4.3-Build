"""Both public HL CLIs preserve the original free-set value and seed contract."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('script,seeded',[
    ('run_phasebuild_hybrid_hl.py',True),('run_build_from_phases.py',False)])
@pytest.mark.parametrize('free_value',[0,7])
def test_cli_passes_declared_free_value_to_fixed_prior(monkeypatch,tmp_path,script,seeded,free_value):
    spec=importlib.util.spec_from_file_location('entry_'+Path(script).stem,ROOT/script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    calls={}
    def factory(path,**kwargs):
        calls['prior']=(path,kwargs)
        return lambda *a,**k: None
    def recycle(*args,**kwargs):
        calls['recycle']=(args,kwargs)
        return {'primary_metric':{'r_work':.2,'r_free':.3},'primary_model':{'path':'primary_model.cif'}}
    monkeypatch.setattr(module,'make_fixed_prior_refiner',factory)
    # Restore the process-global adapter after the CLI's existing injection.
    monkeypatch.setattr(module.hybrid_recycle,'refmac_refine',module.hybrid_recycle.refmac_refine)
    monkeypatch.setattr(module.hybrid_recycle,'recycle',recycle)
    args=['--map','map.mtz','--native','native.mtz','--sequence','sequence.fasta',
          '--phase-prior','original.mtz','--phase-prior-sha256','declared-hash',
          '--output',str(tmp_path/'new-job'),'--runtime','runtime','--free-value',str(free_value)]
    if seeded:args+=['--model','seed.cif']
    module.main(args)
    assert calls['prior']==('original.mtz',{'expected_sha256':'declared-hash','free_value':free_value})
    positional,options=calls['recycle']
    assert positional[2]==('seed.cif' if seeded else None)
    assert options['ccp4']=='runtime'
    assert not (tmp_path/'new-job').exists()
