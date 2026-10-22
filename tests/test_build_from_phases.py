"""No-seed construction and unchanged original nonzero free-value contracts."""
import json
from pathlib import Path
from types import SimpleNamespace
import shutil

import gemmi
import numpy as np
import pytest

from phasebuild import hybrid_recycle
from phasebuild.fixed_phase_prior import make_fixed_prior_refiner
from phasebuild.phase_refinement import merge_experimental_phases
from phasebuild.geometry import template
from phasebuild.trace import write_backbone
from test_phase_refinement import make_mtz


@pytest.mark.parametrize('seeded',[False,True])
def test_unseeded_first_build_and_seeded_compatibility_fixed_last_model(tmp_path,monkeypatch,seeded):
    root=tmp_path/'runtime';(root/'bin').mkdir(parents=True)
    for name in ['cbuccaneer','refmac5']:
        path=root/'bin'/name;path.write_text('synthetic engine');path.chmod(0o755)
    library=root/'lib/data/reference_structures/reference-1tqw';library.parent.mkdir(parents=True)
    library.with_suffix('.mtz').write_text('synthetic resource');library.with_suffix('.pdb').write_text('synthetic resource')
    maps=tmp_path/'map.mtz';source=make_mtz(maps,[('F','F'),('SIGF','Q'),('FWT','F'),('PHWT','P'),('FOM','W')],
        [[1,0,0,10,1,5,30,.5],[2,0,0,20,1,10,60,.5]])
    native=tmp_path/'native.mtz';native.write_bytes(maps.read_bytes())
    sequence=tmp_path/'seq.fasta';sequence.write_text('>A\nAA\n')
    seed=tmp_path/'seed.cif';write_backbone(seed,[template(-63.,-42.,length=2)+[20,20,20]],SimpleNamespace(cell=source.cell,sg=source.spacegroup))
    model=gemmi.read_structure(str(seed));model.write_pdb(str(tmp_path/'builder.pdb'))
    commands=[];refinements=[]
    class FakeProcess:
        def __init__(self,command,**kwargs):
            commands.append(command)
            shutil.copyfile(tmp_path/'builder.pdb',command[command.index('-pdbout')+1]);self.pid=999999
        def wait(self,timeout=None):return 0
        def poll(self):return 0
    def refine(native,model,output,cycles,ccp4=None,timeout=None):
        out=Path(output);out.mkdir();shutil.copyfile(model,out/'model.cif');shutil.copyfile(maps,out/'refined.mtz')
        n=len(refinements)+1;refinements.append((native,model,cycles))
        return dict(statistics=dict(final=dict(cycle=cycles,r_work=.2,r_free=.15+n*.01)),free_flags_preserved=True)
    monkeypatch.setattr(hybrid_recycle.subprocess,'Popen',FakeProcess)
    monkeypatch.setattr(hybrid_recycle,'refmac_refine',refine)
    report=hybrid_recycle.recycle(maps,native,seed if seeded else None,sequence,tmp_path/'output',ccp4=root,progress=lambda _:None)
    assert len(commands)==5 and len(refinements)==5
    assert [c[c.index('-cycles')+1] for c in commands]==['3','2','2','2','2']
    assert [n for _,_,n in refinements]==[10]*5
    if seeded:
        assert commands[0][commands[0].index('-pdbin')+1]==str(seed)
        assert report['engine']=='hybrid-phasebuild-seeded-external-buccaneer-refmac'
        assert report['inputs']['phasebuild_seed_model']['path']==str(seed)
    else:
        assert '-pdbin' not in commands[0]
        assert 'phasebuild_seed_model' not in report['inputs']
        assert report['engine']=='phasebuild-integrated-external-buccaneer-refmac'
        assert 'None' in report['starting_model_source']
    for n,c in enumerate(commands[1:],2):
        assert c[c.index('-pdbin')+1]==str(tmp_path/'output'/('round_%02d_refine'%(n-1))/'model.cif')
    assert report['primary_metric']==report['stages'][-1]['statistics']
    assert report['primary_metric']['r_free']>.16
    assert (tmp_path/'output/primary_model.cif').read_bytes()==(tmp_path/'output/round_05_refine/model.cif').read_bytes()


def test_nonzero_test_value_keeps_original_flags_and_fixed_prior(tmp_path,monkeypatch):
    native=tmp_path/'native.mtz';prior=tmp_path/'prior.mtz'
    make_mtz(native,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],[[1,0,0,10,1,7],[2,0,0,20,1,3]])
    make_mtz(prior,[(k,'A') for k in ['HLA','HLB','HLC','HLD']],[[1,0,0,2,3,1,-1]])
    merged=tmp_path/'merged.mtz';audit=merge_experimental_phases(native,prior,merged,free_value=7)
    assert audit['matched_free_prior_rows']==1 and audit['native_partition']['free_value']==7
    np.testing.assert_array_equal(gemmi.read_mtz_file(str(merged)).column_with_label('FreeR_flag').array,[7,3])
    import phasebuild.fixed_phase_prior as adapter
    def fake_refine(merged,model,output,**kwargs):
        assert kwargs['free_value']==7
        assert kwargs['use_prior'] is True
        Path(output).mkdir()
        return dict(status='success',statistics=dict(final=dict(cycle=10)))
    monkeypatch.setattr(adapter,'refmac_with_experimental_prior',fake_refine)
    refiner=make_fixed_prior_refiner(prior,free_value=7)
    assert refiner.free_value==7
    result=refiner(native,tmp_path/'unused_model.cif',tmp_path/'refine')
    assert result['phase_merge_audit']['native_partition']['free_value']==7


def test_refmac_receives_original_nonzero_test_value_without_reencoding(tmp_path,monkeypatch):
    import phasebuild.phase_refinement as adapter
    native=tmp_path/'native.mtz'
    original=make_mtz(native,[('F','F'),('SIGF','Q'),('FreeR_flag','I'),('SAD_HLA','A'),('SAD_HLB','A'),('SAD_HLC','A'),('SAD_HLD','A')],
        [[1,0,0,10,1,7,2,3,1,-1],[2,0,0,20,1,3,np.nan,np.nan,np.nan,np.nan]])
    model=tmp_path/'model.cif';write_backbone(model,[template(-63,-42,length=2)+[20,20,20]],SimpleNamespace(cell=original.cell,sg=original.spacegroup))
    root=tmp_path/'runtime';(root/'bin').mkdir(parents=True);(root/'bin/refmac5').write_text('stub');(root/'bin/refmac5').chmod(0o755)
    def fake_run(command,**kwargs):
        assert '\nFREE 7\n' in kwargs['input']
        out=Path(kwargs['cwd']);shutil.copyfile(out/'input.pdb',out/'refined.pdb')
        shutil.copyfile(out/'input.mtz',out/'refined.mtz')
        (out/'refmac.xml').write_text('<REFMAC><err_level>0</err_level><Overall_stats>'
            '<n_reflections_all>2</n_reflections_all><n_reflections_free>1</n_reflections_free>'
            '<resolution_low>40</resolution_low><resolution_high>20</resolution_high>'
            '<stats_vs_cycle><new_cycle><cycle>0</cycle><r_factor>0.437</r_factor><r_free>0.433</r_free></new_cycle></stats_vs_cycle>'
            '</Overall_stats></REFMAC>')
        kwargs['stdout'].write('Residual : Phase modified Maximum Likelihood for Fs\n'
            'Ncyc    Rfact    Rfree\n0 0.4375 0.4333 0.58 3000 150 0.0148 1.27 1.56 0.89 0.1\n$TEXT:')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(adapter.subprocess,'run',fake_run)
    result=adapter.refmac_with_experimental_prior(native,model,tmp_path/'out',cycles=0,ccp4=root,free_value=7)
    assert result['partition']['free_value']==7 and result['free_flags_preserved']
    np.testing.assert_array_equal(gemmi.read_mtz_file(str(tmp_path/'out/refined.mtz')).column_with_label('FreeR_flag').array,[7,3])
