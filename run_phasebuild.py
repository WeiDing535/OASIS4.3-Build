"""Run from a source checkout, including with CCP4's isolated Python launcher."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.cli import main

if __name__ == '__main__':
    main()
