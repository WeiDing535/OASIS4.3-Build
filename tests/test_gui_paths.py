from urllib.parse import quote
import pytest
from phasebuild.gui_paths import browse_path


def test_pasted_file_url_and_finder_quotes(tmp_path):
    folder=tmp_path/'space folder';folder.mkdir();file=folder/'sample.mtz';file.touch()
    for value in [str(file),'"'+str(file)+'"','file://'+quote(str(file))]:
        directory,selected,fallback=browse_path(value)
        assert directory==folder and selected==str(file) and not fallback


def test_explicit_go_does_not_silently_fall_back(tmp_path):
    with pytest.raises(ValueError,match='路径不存在'):
        browse_path(tmp_path/'typo')
    directory,selected,fallback=browse_path(tmp_path/'new_output',allow_parent=True)
    assert directory==tmp_path and selected is None and fallback
