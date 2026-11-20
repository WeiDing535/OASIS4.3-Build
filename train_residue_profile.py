import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density
from phasebuild.residue_profile import train_residue_profile


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mtz',required=True)
    p.add_argument('--structure',required=True)
    p.add_argument('--resolution',type=float,default=1.9)
    p.add_argument('--output',required=True)
    args = p.parse_args()
    if Path(args.output).exists():
        raise ValueError('Output already exists')
    density = Density(args.mtz,'FC_BEST.F_phi.F','FC_BEST.F_phi.phi',resolution=args.resolution)
    target = train_residue_profile(density,args.structure,args.mtz)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    target.save(args.output)
    Path(args.output+'.json').write_text(json.dumps(target.metadata,indent=2)+'\n')
    print(json.dumps(target.metadata,indent=2))


if __name__=='__main__':
    main()
