"""Contracts for the user-pinned external CCP4 engine discovery.

OASIS ships no third-party engines. These tests lock in fail-closed behaviour
so a missing/mistyped path never silently falls back to another installation.
"""
from pathlib import Path
import pytest
from phasebuild import engine_paths as ep


def _exe(root, name):
    p = root / 'bin' / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('#!/bin/sh\nexit 0\n')
    p.chmod(0o755)
    return p


def test_ccp4_explicit_path_validates_binaries_and_data(tmp_path):
    root = tmp_path / 'ccp4'
    for name in ep.CCP4_PROGRAMS:
        _exe(root, name)
    for rel in ep.CCP4_DATA:
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
    assert ep.ccp4_directory(str(root)) == str(root.resolve())


def test_ccp4_explicit_incomplete_path_fails_closed(tmp_path, monkeypatch):
    # Even if a complete CCP4 exists elsewhere, an explicit incomplete path
    # must not silently fall back to it.
    other = tmp_path / 'other'
    for name in ep.CCP4_PROGRAMS:
        _exe(other, name)
    monkeypatch.setenv('CCP4', str(other))
    broken = tmp_path / 'broken'
    _exe(broken, 'cbuccaneer')  # only one program
    with pytest.raises(ValueError, match='missing required item'):
        ep.ccp4_directory(str(broken), required=('phaser', 'cbuccaneer'), data=False)


def test_detect_ccp4_reports_incomplete_install_without_raising(tmp_path):
    root = tmp_path / 'ccp4'
    _exe(root, 'phaser')
    found, message = ep.detect_ccp4(str(root))
    assert found == ''
    assert 'missing required item' in message
