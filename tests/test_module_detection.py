"""CCP4 module detection: found, partly found, absent, and "installed elsewhere"."""
import os
from pathlib import Path
import pytest
from phasebuild import module_detection


def ccp4(tmp_path,names=('cparrot','cbuccaneer','refmac5'),root='ccp4'):
    """A minimal install tree: only the executables matter to detection."""
    home=tmp_path/root
    (home/'bin').mkdir(parents=True)
    for name in names:
        path=home/'bin'/name
        path.write_text('#!/bin/sh\nexit 0\n')
        path.chmod(0o755)
    return home


@pytest.fixture(autouse=True)
def no_host_install(monkeypatch):
    """Ignore this machine's real CCP4 so each case sees only its own tree."""
    monkeypatch.setattr(module_detection,'INSTALL_PATTERNS',())


def detect(**kwargs):
    kwargs.setdefault('environ',{})
    return module_detection.detect(**kwargs)


def test_complete_install_switches_both_modules_on(tmp_path):
    home=ccp4(tmp_path)
    report=detect(location=str(home))
    assert report['available'] is True and report['ccp4_ready'] is True
    assert report['modules']['parrot']['available'] and report['modules']['build']['available']
    assert report['modules']['build']['refmac'].endswith('/bin/refmac5')
    assert report['location']['resolved']==str(home)
    assert '已检测到 CCP4' in report['summary']


def test_parrot_without_buccaneer_leaves_building_off(tmp_path):
    home=ccp4(tmp_path,names=('cparrot',))
    report=detect(location=str(home))
    assert report['modules']['parrot']['available'] is True
    assert report['modules']['build']['available'] is False
    assert report['ccp4_ready'] is False
    assert '只检测到 CCP4 的一部分' in report['summary']
    assert 'refmac5' in report['modules']['build']['reason']


def test_nothing_installed_is_a_plain_answer_not_an_error(tmp_path):
    report=detect(location=str(tmp_path/'nowhere'))
    assert report['available'] is False and report['ccp4_ready'] is False
    assert '未检测到 CCP4' in report['summary']
    assert 'CCP4 位置' in report['summary']      # tells the user where to fix it


def test_bin_directory_is_accepted_and_reported_as_its_root(tmp_path):
    home=ccp4(tmp_path)
    report=detect(location=str(home/'bin'))
    assert report['modules']['parrot']['available'] is True
    assert report['location']['resolved']==str(home)


def test_a_setup_script_is_read_for_the_installation_root(tmp_path):
    home=ccp4(tmp_path)
    script=tmp_path/'ccp4.setup'
    script.write_text('# generated\nexport CCP4="%s"\nexport CBIN="%s/bin"\nexport PATH="$CBIN:$PATH"\n'
                      % (home,home))
    report=detect(location=str(script))
    assert report['modules']['build']['available'] is True
    assert report['location']['resolved']==str(home)
    assert '已从 ccp4.setup 解析出 CCP4 位置' in report['location']['note']


def test_a_script_without_ccp4_assignments_says_so(tmp_path):
    script=tmp_path/'bashrc'
    script.write_text('export EDITOR=vim\n')
    report=detect(location=str(script))
    assert report['available'] is False
    assert '没有找到 CCP4' in report['location']['note']


def test_environment_variables_still_work_without_a_typed_location(tmp_path):
    home=ccp4(tmp_path)
    report=detect(environ={'CCP4':str(home)})
    assert report['ccp4_ready'] is True
    assert report['location']['value']=='' and report['location']['resolved']==str(home)


def test_the_run_directory_is_searched_before_the_environment(tmp_path):
    bundled=ccp4(tmp_path,names=('cparrot',),root='bundled')
    other=ccp4(tmp_path,names=('cparrot','cbuccaneer','refmac5'),root='other')
    report=detect(runtime=str(bundled),environ={'CCP4':str(other)})
    assert report['modules']['parrot']['path'].startswith(str(bundled))


def test_config_summary_uses_the_run_directory(tmp_path):
    home=ccp4(tmp_path)
    report=module_detection.summary({'runtime':str(home)})
    assert report['ccp4_ready'] is True

def test_a_wrong_typed_path_is_not_reported_as_the_place_that_answered(tmp_path,monkeypatch):
    """A typo must not make the window claim the run will use that directory."""
    other=ccp4(tmp_path,names=('cparrot','cbuccaneer','refmac5'),root='installed')
    monkeypatch.setattr(module_detection,'INSTALL_PATTERNS',(str(other/'bin'),))
    report=detect(location=str(tmp_path/'typo'))
    assert report['available'] is True                    # found somewhere else
    assert report['location']['resolved']==''             # but not where we looked
    assert report['location']['found_at']==str(other)


def test_found_at_names_the_install_even_without_a_typed_location(tmp_path):
    home=ccp4(tmp_path)
    report=detect(environ={'CBIN':str(home/'bin')})
    assert report['location']['found_at']==str(home)
    assert report['location']['resolved']==str(home)
