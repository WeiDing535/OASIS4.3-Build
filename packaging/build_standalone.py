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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-source', type=Path, required=True)
    parser.add_argument('--runtime-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--zip', type=Path)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Refusing to overwrite bundle: ' + str(output))
    if args.zip and args.zip.exists():
        raise FileExistsError('Refusing to overwrite ZIP: ' + str(args.zip))
    if platform.system() != 'Darwin' or platform.machine() != 'x86_64':
        raise RuntimeError('This runtime combination is validated for Intel macOS only')
    for source in [args.python_source.resolve(), args.runtime_source.resolve()]:
        try:
            output.relative_to(source)
        except ValueError:
            continue
        raise ValueError('Output must not be inside a source runtime')
    report = {'product': 'OASIS4.3-Build', 'version': '0.8.3',
              'platform': 'macOS-x86_64', 'distribution_status': 'local-evaluation',
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
    assembler.main(['--output', str(output)])
    for label, source, name in [('python', args.python_source, 'python_manifest.json'),
                                ('runtime', args.runtime_source, 'runtime_manifest.json')]:
        print('Copy and verify ' + label, flush=True)
        shutil.copytree(source, output / label, symlinks=True)
        report['validation'][label + '_copied'] = verify_runtime(output / label, name)
    # A product-named launcher beside the original compatibility launcher.
    shutil.copy2(output / 'PhaseBuild.command', output / 'OASIS4.3-Build.command')
    (output / '使用说明.txt').write_text(
        'OASIS4.3-Build 0.8.3 — Intel Mac 完整本机试用包\n\n'
        '解压整个目录到可写位置，双击 OASIS4.3-Build.command。\n'
        '无需另装 CCP4、Python、NumPy、SciPy 或 Gemmi。保持整个目录结构。\n'
        '界面内全部六种模式所需运行组件已随包提供；用户另行提供实验输入文件。\n'
        '已包含外部 Buccaneer/REFMAC，并不表示其算法由本项目原创。\n'
        '当前仅验证 Intel Mac；未签名/公证，未验证 Apple Silicon、Windows 或 Linux。\n'
        '关闭网页不会停止计算；界面停止任务，终端 Ctrl+C 关闭服务。\n'
        '结果写入新目录，不覆盖已有结果。内置示例按钮仍需用户自行提供对应示例数据。\n\n'
        '这是从本机现有授权安装构建的本机评估包。\n'
        '完整引擎、关联库和资源的公开再分发条件尚未落实，不是公开发行包。\n')
    shutil.copy2(ROOT / 'LICENSE', output / 'LICENSE')
    env = {'PATH': '/usr/bin:/bin:/usr/sbin:/sbin', 'LANG': 'C',
           'HOME': str(output.parent / (output.name + '-validation-home'))}
    print('Verify bundled doctor and CLI entry points', flush=True)
    doctor = subprocess.run([str(output / 'phasebuild'), 'doctor', '--verify'],
                            env=env, capture_output=True, text=True, timeout=180)
    if doctor.returncode:
        raise RuntimeError('Bundled doctor failed: ' + doctor.stdout + doctor.stderr)
    report['validation']['doctor'] = json.loads(doctor.stdout)
    modes = ['gui', 'native', 'hybrid-hl', 'build-hl', 'hybrid', 'refine', 'rfactor', 'phase-refine', 'audit', 'evaluate']
    report['validation']['cli_help'] = {}
    for mode in modes:
        result = subprocess.run([str(output / 'phasebuild'), mode, '--help'], env=env,
                                capture_output=True, text=True, timeout=60)
        if result.returncode:
            raise RuntimeError('Bundled CLI failed: ' + mode + '\n' + result.stderr)
        report['validation']['cli_help'][mode] = True
    report['status'] = 'assembled-and-verified'
    (output / 'bundle.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    if args.zip:
        print('Create ZIP', flush=True)
        subprocess.run(['/usr/bin/ditto', '-c', '-k', '--keepParent', str(output), str(args.zip.resolve())], check=True)
        with zipfile.ZipFile(args.zip) as archive:
            bad = archive.testzip()
            if bad:
                raise RuntimeError('ZIP CRC failed: ' + bad)
        args.zip.with_suffix(args.zip.suffix + '.sha256').write_text(digest(args.zip) + '  ' + args.zip.name + '\n')
    print(json.dumps({'bundle': str(output), 'status': report['status'],
                      'public_redistribution_ready': False}, indent=2))


if __name__ == '__main__':
    main()
