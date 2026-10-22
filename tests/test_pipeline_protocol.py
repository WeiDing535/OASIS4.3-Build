import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from run_buccaneer_pipeline import main
from phasebuild.data import fingerprint


def test_pipeline_rejects_reusing_a_first_round_with_a_different_free_set(tmp_path,monkeypatch):
    native=tmp_path/'native.mtz';native.write_bytes(b'input fingerprint only')
    seq=tmp_path/'sequence.fasta';seq.write_text('>test\nACDEFG\n')
    build=tmp_path/'build';build.mkdir()
    refine=tmp_path/'refine';refine.mkdir()
    (build/'build_report.json').write_text(json.dumps({'sequence':fingerprint(seq),
        'output':{'sha256':'same model'},'cycles':3,'correlation_mode':True}))
    (refine/'refinement_report.json').write_text(json.dumps({'reflections':{'input':fingerprint(native),
        'fp':'F','sigfp':'SIGF','free_column':'FreeR_flag','free_value':1},
        'model_input':{'sha256':'same model'},'requested_cycles':10}))
    monkeypatch.setattr(sys,'argv',['pipeline','--native',str(native),'--sequence',str(seq),
        '--first-build',str(build),'--first-refinement',str(refine),'--output',str(tmp_path/'out')])
    with pytest.raises(ValueError,match='free-set protocol'):
        main()
