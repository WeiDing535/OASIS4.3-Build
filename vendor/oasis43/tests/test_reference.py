import sys
from pathlib import Path
import pytest
from oasis_py.reference import run_reference
pytestmark=pytest.mark.skipif(sys.platform=="win32",reason="POSIX fake executable adapter tests")


def script(tmp_path,body):
    path=tmp_path/'fake_engine'
    path.write_text('#!'+sys.executable+'\n'+body)
    path.chmod(0o755)
    key=tmp_path/'input.key';key.write_text('END\n')
    return path,key


def test_reference_scratch_is_isolated_and_input_is_logged(tmp_path):
    exe,key=script(tmp_path,"import os,sys\nfrom pathlib import Path\nprint(os.environ['CCP4_SCR'])\nprint(sys.stdin.read())\n")
    out=tmp_path/'job'
    result=run_reference(exe,key,out)
    assert result['returncode']==0
    log=(out/'reference.log').read_text()
    assert str(out/'scratch') in log and 'END' in log
    with pytest.raises(ValueError,match='empty'):run_reference(exe,key,out)


def test_reference_missing_output_is_a_failure(tmp_path):
    exe,key=script(tmp_path,'pass\n')
    with pytest.raises(RuntimeError,match='expected output'):
        run_reference(exe,key,tmp_path/'job',hklin=tmp_path/'unused.mtz')


def test_reference_timeout_is_bounded_and_reports_log(tmp_path):
    exe,key=script(tmp_path,'import time\ntime.sleep(30)\n')
    with pytest.raises(RuntimeError,match='was terminated'):
        run_reference(exe,key,tmp_path/'job',timeout=.2)
