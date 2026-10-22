"""Evaluate fixed-coordinate Rwork/Rfree without running REFMAC."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.rfactors import calculate_rfactors


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--output',required=True)
    p.add_argument('--amplitude',default='F');p.add_argument('--sigma',default='SIGF')
    p.add_argument('--free-column',default='FreeR_flag');p.add_argument('--free-value',type=int,default=0)
    p.add_argument('--method',choices=['fft','direct'],default='fft');p.add_argument('--no-solvent',action='store_true')
    a=p.parse_args()
    calculate_rfactors(a.mtz,a.model,a.output,a.amplitude,a.sigma,a.free_column,a.free_value,
                      bulk_solvent=not a.no_solvent,method=a.method)


if __name__=='__main__':
    main()
