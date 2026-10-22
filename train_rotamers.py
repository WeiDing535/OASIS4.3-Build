import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.rotamers import train_rotamers,save_library


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--structure',required=True)
    p.add_argument('--monomers',default='/Applications/ccp4-9/lib/data/monomers')
    p.add_argument('--output',required=True)
    args=p.parse_args()
    if Path(args.output).exists():
        raise ValueError('Output exists')
    library=train_rotamers(args.structure,args.monomers)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    save_library(library,args.output)
    print({name:len(entry['conformers']) for name,entry in library['entries'].items()})


if __name__=='__main__':
    main()
