"""Linux loader auditing without requiring Linux binaries on the test host."""
import importlib.util
from pathlib import Path
import pytest
SPEC = importlib.util.spec_from_file_location('linux_runtime_audit', Path(__file__).resolve().parents[1] / 'packaging/audit_linux_runtime.py')
AUDIT = importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(AUDIT)


def test_ldd_parser_rejects_missing_libraries():
    with pytest.raises(ValueError, match='Missing Linux library'):
        AUDIT.parse_ldd('libclipper.so => not found')
    result = AUDIT.parse_ldd('linux-vdso.so.1 (0x1234)\nlibc.so.6 => /lib/libc.so.6 (0x12)\n/lib64/ld-linux-x86-64.so.2 (0x34)')
    assert result == [('libc.so.6', Path('/lib/libc.so.6')), ('ld-linux-x86-64.so.2', Path('/lib64/ld-linux-x86-64.so.2'))]


def test_custom_libraries_cannot_fall_back_to_host(tmp_path):
    root = tmp_path / 'runtime';root.mkdir()
    AUDIT.check_dependencies(root, [('libc.so.6', Path('/usr/lib/libc.so.6'))])
    with pytest.raises(ValueError, match='host installation'):
        AUDIT.check_dependencies(root, [('libclipper.so', Path('/host/libclipper.so'))])
    with pytest.raises(ValueError, match='outside OS library'):
        AUDIT.check_dependencies(root, [('libc.so.6', Path('/host/libc.so.6'))])


def test_rpaths_must_be_relative_and_stay_inside_prefix(tmp_path):
    root = tmp_path / 'runtime';(root / 'bin').mkdir(parents=True)
    owner = root / 'bin/tool'
    AUDIT.check_rpaths(root, owner, '(RUNPATH) Library runpath: [$ORIGIN/../lib]')
    for value in ['/opt/ccp4/lib', '$ORIGIN/../../outside', ':$ORIGIN']:
        with pytest.raises(ValueError):
            AUDIT.check_rpaths(root, owner, '(RPATH) Library rpath: [' + value + ']')
