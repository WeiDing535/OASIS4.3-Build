import json
from pathlib import Path
from types import SimpleNamespace
import pytest
import phasebuild.iterative_cleanup as module


def setup(tmp_path,monkeypatch,*,refused=False):
 paths={name:tmp_path/(name+'.dat') for name in ['dm_map','native','sequence','model']}
 for name,p in paths.items():p.write_text(name)
 monkeypatch.setattr(module,'read_sequences',lambda path:[('one','A'*80)])
 monkeypatch.setattr(module,'work_density',lambda *a:('grid',{'measured_work_coefficients':100,'free_coefficients_removed':10}))
 class Document:
  def write_file(self,path):Path(path).write_text('selected atoms')
 clean=SimpleNamespace(make_mmcif_document=lambda:Document())
 monkeypatch.setattr(module.gemmi,'read_structure',lambda path:'constructed model')
 def choose(*args,**kwargs):
  if refused:raise ValueError('Insufficient uniquely anchored sequence coverage')
  return clean,dict(operation='test-selection',selection_uses_Rfree=False,coordinate_generation=False)
 monkeypatch.setattr(module,'cleanup_partial_single_copy',choose)
 return paths

@pytest.mark.parametrize('refused',[False,True])
def test_wrap_keeps_original_and_declares_gate_decision(tmp_path,monkeypatch,refused):
 paths=setup(tmp_path,monkeypatch,refused=refused);seen=[]
 def refiner(native,model,output,*args,**kwargs):
  seen.append(Path(model));Path(output).mkdir();return {'statistics':{'final':{'r_free':.99}}}
 wrapped=module.wrap_refiner(refiner,paths['dm_map'],paths['native'],paths['sequence'],asu_copies=1)
 out=tmp_path/'round_01_refine';report=wrapped(paths['native'],paths['model'],out)
 assert paths['model'].read_text()=='model'
 assert (seen[0]==paths['model'])==refused
 assert report['pre_refinement_cleanup']['status']==('refused' if refused else 'completed')
 assert json.loads((out/'refinement_report.json').read_text())['pre_refinement_cleanup']['selection_uses_Rfree'] is False


def test_changed_reference_or_wrong_native_fails_before_output(tmp_path,monkeypatch):
 paths=setup(tmp_path,monkeypatch)
 wrapped=module.wrap_refiner(lambda *a:pytest.fail('must not refine'),paths['dm_map'],paths['native'],paths['sequence'],asu_copies=1)
 wrong=tmp_path/'different.mtz';wrong.write_text('different measurements');out=tmp_path/'out'
 with pytest.raises(ValueError,match='native differs'):wrapped(wrong,paths['model'],out)
 assert not out.exists() and not out.with_name(out.name+'_sequence_cleanup').exists()
 paths['dm_map'].write_text('changed')
 with pytest.raises(RuntimeError,match='input changed'):wrapped(paths['native'],paths['model'],out)
 assert not out.exists()

@pytest.mark.parametrize('copies',[0,2,True,1.0])
def test_ambiguous_copy_count_rejected(tmp_path,copies):
 with pytest.raises(ValueError,match='single-copy'):module.wrap_refiner(None,'unused','unused','unused',asu_copies=copies)
