# Phases-to-model integration entry points

This repository provides two entry points: `run_build_from_phases.py` (no initial model) and `run_phasebuild_hybrid_hl.py` (with a seed model). Run them in an existing Python environment and add `--help` to list arguments.

Both entries require an existing phase map, native data, sequence, original experimental HL and a new output directory; the seed entry additionally requires a model. Original free-set values are declared as-is via `--free-value`. The experimental prior stays fixed; missing values are not filled in with model phases. See GUI_GUIDE.md for column names and index constraints.

External Buccaneer performs building and REFMAC performs refinement. The native algorithm research module is kept separately — do not label external integration results as native performance. This source release includes no historical real-data run records and makes no new structure-solving claims for the new entry points.
