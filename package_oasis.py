"""Create a focused OASIS4.3 source package, excluding retired product entry points."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parent
SCRIPTS=['run_gui.py','run_sad_to_model.py','run_oasis_stage.py','run_build_from_phases.py']
# This compatibility adapter is only used by the old workbench in the development
# repository. The focused server always receives oasis_workflow explicitly.
LEGACY_ONLY={'phasebuild.gui_workflow'}


def source_files(root=ROOT):
    root=Path(root);src=root/'src';paths=set();pending=[root/p for p in SCRIPTS]
    def module_path(name):
        path=src/Path(*name.split('.'))
        return path.with_suffix('.py') if path.with_suffix('.py').is_file() else path/'__init__.py' if (path/'__init__.py').is_file() else None
    def include(name):
        if name in LEGACY_ONLY:return
        path=module_path(name)
        if path and path not in paths:
            pending.append(path)
        parts=name.split('.')
        for i in range(1,len(parts)):
            init=src/Path(*parts[:i])/'__init__.py'
            if init.is_file() and init not in paths:pending.append(init)
    while pending:
        path=pending.pop()
        if path in paths:continue
        paths.add(path)
        package='.'.join(path.relative_to(src).parts[:-1]) if src in path.parents else ''
        for item in ast.walk(ast.parse(path.read_text())):
            if isinstance(item,ast.ImportFrom):
                if item.level:
                    components=package.split('.')
                    prefix='.'.join(components[:len(components)-item.level+1])
                    name=prefix+('.'+item.module if item.module else '')
                else:name=item.module or ''
                if name.split('.')[0] not in ('phasebuild','oasis_py'):continue
                include(name)
                for alias in item.names:include(name+'.'+alias.name)
            elif isinstance(item,ast.Import):
                for alias in item.names:
                    if alias.name.split('.')[0] in ('phasebuild','oasis_py'):include(alias.name)
    paths.update((src/'phasebuild/oasis_assets').glob('*'))
    paths.update(p for p in (src/'oasis_py/resources').rglob('*') if p.is_file())
    paths.update(root/p for p in ['OASIS4.3.command','environment.yml','LICENSE','AUTHORS.md','CITATION.cff','docs/OASIS4.3.md',
                   'vendor/oasis43/LICENSE','vendor/oasis43/UPSTREAM.json','vendor/oasis43/TEST_PATCHES.md'])
    return sorted(paths)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    target=parser.parse_args(argv).output
    target.parent.mkdir(parents=True,exist_ok=True)
    prefix='OASIS4.3-focused-development'
    files=source_files()
    manifest=dict(product='OASIS4.3',distribution='focused source development',
                  entrypoints=SCRIPTS,methods=['likelihood','pplus'],external_runtime_bundled=False,files=[])
    metadata='''[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"
[project]
name = "oasis43"
version = "4.3.0.dev1"
description = "OASIS4.3: two SAD phase methods, density modification, Buccaneer and REFMAC"
requires-python = ">=3.12,<3.13"
dependencies = ["gemmi==0.7.5", "numpy>=1.26,<2", "scipy>=1.12"]
[tool.setuptools.packages.find]
where = ["src"]
[tool.setuptools.package-data]
phasebuild = ["oasis_assets/*"]
oasis_py = ["resources/*.json", "resources/*.md", "resources/*.gz", "resources/licenses/*.txt"]
'''
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED) as z:
        for path in files:
            if not path.is_file():raise FileNotFoundError(path)
            relative=str(path.relative_to(ROOT));payload=path.read_bytes()
            manifest['files'].append(dict(path=relative,sha256=hashlib.sha256(payload).hexdigest()))
            z.write(path,prefix+'/'+relative)
        z.writestr(prefix+'/pyproject.toml',metadata)
        z.writestr(prefix+'/README.md',(ROOT/'docs/OASIS4.3.md').read_text())
        z.writestr(prefix+'/PACKAGE_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        assert z.getinfo(prefix+'/OASIS4.3.command').external_attr>>16&0o111
    print(json.dumps(dict(path=str(target),files=len(files),bytes=target.stat().st_size),indent=2))

if __name__=='__main__':main()
