"""Command routing must keep the map and profile arguments intact after wrapping."""
import importlib.util
from pathlib import Path
from unittest.mock import patch
import pytest
spec=importlib.util.spec_from_file_location('sad_development',Path(__file__).resolve().parents[1]/'run_sad_development.py')
cli=importlib.util.module_from_spec(spec);spec.loader.exec_module(cli)
@pytest.mark.parametrize('script',['run_context_seeds','run_weak_trace','run_complete','run_rebuild'])
def test_profile_resolution_selection_and_adaptive_stage_preserve_map_argument(script):
 command=['python','-u',str(cli.ROOT/(script+'.py')),'input.mtz','--profile','lowres.npz','--output','output']
 if script in ['run_complete','run_rebuild']:command+=['--residue-profile','residue-lowres.npz']
 class Mtz:
  def resolution_high(self):return 1.9
 with patch('gemmi.read_mtz_file',return_value=Mtz()) as read:
  result=cli.stage_command({'command':command})
 read.assert_called_once_with('input.mtz')
 assert result[result.index('--profile')+1].endswith('1tqw_1p9_symmetric.npz')
 if script in ['run_complete','run_rebuild']:
  index=result.index(str(cli.ROOT/'run_adaptive_stage.py'))
  assert result[index+1:index+3]==[script[4:],'input.mtz']
  assert result[result.index('--residue-profile')+1].endswith('residue_1tqw_1p9.npz')
 else:assert 'input.mtz'==result[3]
 assert command[command.index('--profile')+1]=='lowres.npz'
