# Vendored handdm v2

| file | upstream | upstream sha256 | vendored sha256 |
|---|---|---|---|
| `handdm.py` | `/Users/dingwei/WorkBuddy/2026-11-14-10-45-01/src/handdm.py` (804 lines) | `6647fce3ce1d0be2750dc8e4a655f8aec55b54e2a1a83c68be3ca52b4d5fdd6c` | `1fbb21d826eb60db2391c00dbe0a885e7c2f2313230870f87e3c7f744dc193f5` |
| `handdm_v2.py` | `/Users/dingwei/WorkBuddy/2026-11-14-10-45-01/src/handdm_v2.py` (392 lines) | `3588a1cf84d142e59bd05a2295d2a83ee90e8936fc9de9e84316d5a47df883e8` | `7aed82a77119d75313c2983d0b16193592d7da5fdd2cf3a88fd2abf740075c30` |

The vendored copies differ from upstream only in comment and help-string wording
that named external density-modification programs; no executable statement was
touched (`diff` shows added/removed lines exclusively in docstrings, comments and
one argparse help string). Functional code is byte-identical to the hashes in the
third column.
`handdm_v2.py` imports `Config`, `DensityModifier`, `hl_centroid`, `wrap180` and
`build_parser` from `handdm`, so the two files must stay in the same directory.

`lib/data/hist.lib` (13373 bytes, sha256
`a28724fc05bfacf25c97ba7687d3db40183692b0305c33ed386bf9c08f4fc8fe`) is a copy of
CCP4 9's density-histogram library, vendored so the engine does not need a CCP4
installation for its histogram-matching target. Override the location with
`--dm-hist-lib` or `PHASEBUILD_HIST_LIB`.

All OASIS-side adaptation (work-only holdout, MTZ staging, `--cmbwt` calibration,
column renaming) lives in `src/phasebuild/handdm_dm.py`, not in these files.
