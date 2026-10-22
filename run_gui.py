"""Launch PhaseBuild's local browser GUI from a source checkout."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.gui import main

if __name__=='__main__':
    main()
