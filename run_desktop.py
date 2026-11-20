"""Launch the OASIS4.3 desktop window from a source checkout."""
import sys
from pathlib import Path
import os
import platform
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_thread_key]='1'
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))

# Platform note (entry layer only, does not touch core algorithms): the GUI is pure Tkinter and runs on both macOS and Linux;
# the external CCP4 engine is user-supplied. Linux needs no special configuration beyond a display server.
if sys.platform=='darwin' and platform.machine()=='arm64':
    print('Apple Silicon detected; make sure Rosetta 2 is installed, otherwise the external engine cannot run',file=sys.stderr)

from phasebuild.desktop import main

if __name__=='__main__':
    main()
