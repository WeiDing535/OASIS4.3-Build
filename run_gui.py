"""Launch PhaseBuild's local browser GUI from a source checkout."""
import sys
from pathlib import Path
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_thread_key]='1'
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.gui import main

if __name__=='__main__':
    main()
