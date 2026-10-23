"""Assemble a complete local-evaluation bundle from verified existing runtimes.

Never installs dependencies, edits the source runtime, or authorizes public
redistribution. The destination and optional ZIP must be new paths.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify_runtime(root, manifest_name):
    root = Path(root).resolve()
    manifest = root / manifest_name
    records = json.loads(manifest.read_text())['files']
    for record in records:
        relative = record['relative_path']
        path = root / relative
        try:
            path.resolve().relative_to(root)
        except ValueError:
            raise ValueError('Runtime file escapes bundle: ' + relative)
        if not path.is_file():
            raise ValueError('Missing runtime file: ' + relative)
        if digest(path) != record['sha256']:
            raise ValueError('Runtime hash differs: ' + relative)
    # Check all links, including directory links not covered by regular files.
    for path in root.rglob('*'):
        if path.is_symlink():
            try:
                path.resolve(strict=True).relative_to(root)
            except (ValueError, OSError, RuntimeError):
                raise ValueError('Runtime symlink escapes or is broken: ' + str(path))
    return {'files_checked': len(records), 'manifest_sha256': digest(manifest)}



def target_platform():
    machine = platform.machine().lower()
    if machine not in ('x86_64', 'amd64'):
        raise RuntimeError('Only x86_64 bundle assembly is implemented')
    if platform.system() == 'Linux':
        return 'linux-x86_64'
    if platform.system() == 'Darwin':
        return 'macos-x86_64'
    raise RuntimeError('Bundle assembly requires Linux or macOS')


def verify_linux_engines(root):
    for name in ('cbuccaneer', 'refmac5'):
        with (Path(root) / 'bin' / name).open('rb') as handle:
            header = handle.read(20)
        if (len(header) < 20 or header[:6] != b'\x7fELF\x02\x01'
                or int.from_bytes(header[18:20], 'little') != 62):
            raise ValueError('Expected a Linux x86_64 ELF engine: ' + name)



def write_tar(output, destination):
    output, destination = Path(output), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('xb') as handle:
        with tarfile.open(fileobj=handle, mode='w:gz', dereference=False) as archive:
            archive.add(output, arcname=output.name)
    destination.with_suffix(destination.suffix + '.sha256').write_text(
        digest(destination) + '  ' + destination.name + '\n')



def verify_doctor_paths(report, output):
    output = Path(output).resolve()
    python_root = output / 'python'
    paths = [report['python'], report['python_prefix']]
    paths.extend(item['path'] for item in report['packages'].values())
    for path in paths:
        try:
            Path(path).resolve().relative_to(python_root)
        except ValueError:
            raise ValueError('Bundled Python loaded outside its prefix: ' + str(path))
    if Path(report['runtime']['root']).resolve() != output / 'runtime':
        raise ValueError('Doctor selected a host engine runtime')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-source', type=Path, required=True)
    parser.add_argument('--runtime-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    archives = parser.add_mutually_exclusive_group()
    archives.add_argument('--zip', type=Path)
    archives.add_argument('--tar', type=Path, help='New .tar.gz archive (recommended on Linux)')
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Refusing to overwrite bundle: ' + str(output))
    if args.zip and args.zip.exists():
        raise FileExistsError('Refusing to overwrite ZIP: ' + str(args.zip))
    if args.tar and args.tar.exists():
        raise FileExistsError('Refusing to overwrite archive: ' + str(args.tar))
    target = target_platform()
    if target == 'linux-x86_64' and args.zip:
        raise ValueError('Use --tar for Linux to preserve executable modes and symlinks')
    if target == 'linux-x86_64':
        verify_linux_engines(args.runtime_source)
    for source in [args.python_source.resolve(), args.runtime_source.resolve()]:
        try:
            output.relative_to(source)
        except ValueError:
            continue
        raise ValueError('Output must not be inside a source runtime')
    report = {'product': 'OASIS4.3-Build', 'version': '0.8.3',
              'platform': target, 'distribution_status': 'local-evaluation',
              'public_redistribution_ready': False, 'dependencies_installed': False,
              'requires_separate_ccp4_installation': False, 'validation': {},
              'scope': 'Complete GUI application runtime; user supplies experimental inputs. '
                       'External Buccaneer/REFMAC retained. Public redistribution rights unresolved.'}
    for label, source, name in [('python', args.python_source, 'python_manifest.json'),
                                ('runtime', args.runtime_source, 'runtime_manifest.json')]:
        print('Verify source ' + label, flush=True)
        report['validation'][label + '_source'] = verify_runtime(source, name)
    spec = importlib.util.spec_from_file_location('phasebuild_assembler', ROOT / 'packaging/assemble_bundle.py')
    assembler = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(assembler)
    assembler.main(['--platform', target, '--output', str(output)])
    for label, source, name in [('python', args.python_source, 'python_manifest.json'),
                                ('runtime', args.runtime_source, 'runtime_manifest.json')]:
        print('Copy and verify ' + label, flush=True)
        shutil.copytree(source, output / label, symlinks=True)
        report['validation'][label + '_copied'] = verify_runtime(output / label, name)
    if target == 'linux-x86_64':
        audit_spec = importlib.util.spec_from_file_location('linux_runtime_audit', ROOT / 'packaging/audit_linux_runtime.py')
        audit = importlib.util.module_from_spec(audit_spec)
        audit_spec.loader.exec_module(audit)
        for directory, kind in [('python', 'python'), ('runtime', 'engines')]:
            inspected = audit.audit(output / directory, kind, allow_manifest=True)
            report['validation'][directory + '_linux_loader'] = {
                'elf_files': len(inspected['elf_checks']),
                'glibc_minimum_from_symbols': inspected['glibc_minimum_from_symbols']}
    # A product-named launcher beside the original compatibility launcher.
    launcher = 'OASIS4.3-Build.sh' if target == 'linux-x86_64' else 'OASIS4.3-Build.command'
    shutil.copy2(output / 'PhaseBuild.command', output / launcher)
    (output / '使用说明.txt').write_text(
        'OASIS4.3-Build 0.8.3 — ' + target + ' 完整本机试用包\n\n'
        '解压整个目录到可写位置，运行 ' + launcher + '。\n'
        '无需另装 CCP4、Python、NumPy、SciPy 或 Gemmi。保持整个目录结构。\n'
        '界面内全部六种模式所需运行组件已随包提供；用户另行提供实验输入文件。\n'
        '已包含外部 Buccaneer/REFMAC，并不表示其算法由本项目原创。\n'
        '平台及验证范围见 bundle.json；未签名/公证，不代表支持其他系统或架构。\n'
        'Linux 无桌面服务器可运行 ./phasebuild gui --no-browser --port 8765，使用 SSH 本地端口转发访问。\n'
        '关闭网页不会停止计算；界面停止任务，终端 Ctrl+C 关闭服务。\n'
        '结果写入新目录，不覆盖已有结果。内置示例按钮仍需用户自行提供对应示例数据。\n\n'
        '这是从本机现有授权安装构建的本机评估包。\n'
        '完整引擎、关联库和资源的公开再分发条件尚未落实，不是公开发行包。\n')
    shutil.copy2(ROOT / 'LICENSE', output / 'LICENSE')
    env = {'PATH': '/usr/bin:/bin:/usr/sbin:/sbin', 'LANG': 'C',
           'HOME': str(output.parent / (output.name + '-validation-home'))}
    print('Verify unified Python 3.12 and OASIS engine', flush=True)
    unified = subprocess.run([str(output/'python/bin/python3'), '-c',
        "import sys; sys.path.insert(0,sys.argv[1]); assert sys.version_info[:2]==(3,12), 'Unified package requires Python 3.12'; import oasis_py.cli",
        str(output/'app/src')], env=env, capture_output=True, text=True, timeout=60)
    if unified.returncode:
        raise RuntimeError('Incomplete unified Python runtime: '+unified.stderr)
    report['validation']['unified_python_312'] = True
    print('Verify bundled doctor and CLI entry points', flush=True)
    doctor = subprocess.run([str(output / 'phasebuild'), 'doctor', '--verify'],
                            env=env, capture_output=True, text=True, timeout=180)
    if doctor.returncode:
        raise RuntimeError('Bundled doctor failed: ' + doctor.stdout + doctor.stderr)
    report['validation']['doctor'] = json.loads(doctor.stdout)
    verify_doctor_paths(report['validation']['doctor'], output)
    modes = ['sad-to-model', 'gui', 'native', 'hybrid-hl', 'build-hl', 'hybrid', 'refine', 'rfactor', 'phase-refine', 'audit', 'evaluate']
    report['validation']['cli_help'] = {}
    for mode in modes:
        result = subprocess.run([str(output / 'phasebuild'), mode, '--help'], env=env,
                                capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise RuntimeError('Bundled CLI failed: ' + mode + '\n' + result.stderr)
        report['validation']['cli_help'][mode] = True
    report['status'] = 'assembled-and-verified'
    (output / 'bundle.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    if args.tar:
        print('Create tar.gz', flush=True)
        write_tar(output, args.tar)
    if args.zip:
        print('Create ZIP', flush=True)
        args.zip.parent.mkdir(parents=True, exist_ok=True)
        if target == 'macos-x86_64':
            subprocess.run(['/usr/bin/ditto', '-c', '-k', '--keepParent', str(output), str(args.zip.resolve())], check=True)
        else:
            raise ValueError('Use --tar for Linux to preserve executable modes and symlinks')
        with zipfile.ZipFile(args.zip) as archive:
            bad = archive.testzip()
            if bad:
                raise RuntimeError('ZIP CRC failed: ' + bad)
        args.zip.with_suffix(args.zip.suffix + '.sha256').write_text(digest(args.zip) + '  ' + args.zip.name + '\n')
    print(json.dumps({'bundle': str(output), 'status': report['status'],
                      'public_redistribution_ready': False}, indent=2))


if __name__ == '__main__':
    main()
