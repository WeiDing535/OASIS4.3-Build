"""Package a source distribution without local data, run history or CCP4 binaries."""
import hashlib
import argparse
import json
import zipfile
from pathlib import Path


def main(argv=None):
    root=Path(__file__).resolve().parent
    version='0.8.3'
    prefix='phasebuild-gui-'+version
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=root/'dist'/(prefix+'-source.zip'))
    target=parser.parse_args(argv).output
    target.parent.mkdir(parents=True,exist_ok=True)
    files=list((root/'src').rglob('*.py'))+list((root/'src/phasebuild/gui_assets').glob('*'))
    files+=list(root.glob('run_*.py'))+list(root.glob('train_*.py'))
    files+=[root/name for name in ['pyproject.toml','README.md','GUI_GUIDE.md','SAD_RESEARCH_GUIDE.md','WEAK_PHASE_DEVELOPMENT.md','PhaseBuild.command','package_gui.py','LICENSE','AUTHORS.md','CITATION.cff','README_EN.md','CHANGELOG.md','environment.yml']]
    files+=list((root/'packaging').glob('*.py'))
    files+=list((root/'docs').glob('*.md'))
    files += [p for p in (root/'src/oasis_py/resources').rglob('*') if p.is_file()]
    files += [p for p in (root/'vendor/oasis43').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    files+=[root/'development/INTEGRATION_REPORT.md']
    files+=list((root/'tests').glob('*.py'))
    files+=[root/'models'/name for name in ['1tqw_1p9_symmetric.npz','residue_1tqw_1p9.npz','rotamers_1tqw_v08.json']]
    files+=[root/'runs/sad_four_optimization/models'/name for name in ['backbone_2p8.npz','residue_2p8.npz']]
    manifest=dict(version=version,distribution='source',ccp4_bundled=False,target_data_bundled=False,files=[])
    with zipfile.ZipFile(str(target),'x',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(set(files)):
            if '__pycache__' in path.parts or not path.is_file():
                continue
            relative=str(path.relative_to(root))
            manifest['files'].append(dict(path=relative,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            archive.write(str(path),prefix+'/'+relative)
        archive.writestr(prefix+'/PACKAGE_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(str(target)) as archive:
        assert archive.testzip() is None
        assert archive.getinfo(prefix+'/PhaseBuild.command').external_attr >> 16 & 0o111
    print(str(target))
    print('Files:',len(manifest['files']),'Bytes:',target.stat().st_size)


if __name__=='__main__':
    main()
