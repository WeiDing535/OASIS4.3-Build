# Vendored oasis-dm v2

| file | upstream | upstream sha256 | vendored sha256 |
|---|---|---|---|
| `oasis_dm_core.py` | `/Users/dingwei/WorkBuddy/2026-11-14-10-45-01/src/handdm.py` (804 lines) | `6647fce3ce1d0be2750dc8e4a655f8aec55b54e2a1a83c68be3ca52b4d5fdd6c` | `12d31984bf85a034531d89d268df189df3fe6de625aeda1844a308027589ae1c` |
| `oasis_dm_v2.py` | `/Users/dingwei/WorkBuddy/2026-11-14-10-45-01/src/handdm_v2.py` (392 lines) | `3588a1cf84d142e59bd05a2295d2a83ee90e8936fc9de9e84316d5a47df883e8` | `1d18cce5c66e2997533c8128ca189f730723411c783c31dcb1c20b03d38ad7b1` |

The engine was written under the working name `handdm`; inside OASIS4.3 it is
called **oasis-dm** and its two files are `oasis_dm_core.py` (the base modifier)
and `oasis_dm_v2.py` (the omit-based recombination and free-set early stop). The
upstream paths in the table are the real files it was copied from, so they keep
the original names.

The vendored copies differ from those upstream files only in

* the module rename — `oasis_dm_v2.py` does `from oasis_dm_core import ...`
  instead of `from handdm import ...`;
* docstrings, comments and one argparse help string, where the old program name
  or the name of an external density-modification program appeared.

No numerical statement was touched: `diff -u` against the upstream files shows
only the lines above, and every one of them is inside a comment, a docstring, a
help string or the import target. The hashes in the third and fourth columns
record both sides so a future re-vendoring shows a deliberate difference.
`oasis_dm_v2.py` imports `Config`, `DensityModifier`, `hl_centroid`, `wrap180`
and `build_parser` from `oasis_dm_core`, so the two files must stay in the same
directory.

`lib/data/hist.lib` (13373 bytes, sha256
`a28724fc05bfacf25c97ba7687d3db40183692b0305c33ed386bf9c08f4fc8fe`) is a copy of
CCP4 9's density-histogram library, vendored so the engine does not need a CCP4
installation for its histogram-matching target. Override the location with
`--dm-hist-lib` or `PHASEBUILD_HIST_LIB`.

All OASIS-side adaptation (work-only holdout, MTZ staging, `--cmbwt` calibration,
column renaming) lives in `src/phasebuild/oasis_dm.py`, not in these files.
