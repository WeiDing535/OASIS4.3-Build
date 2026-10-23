# Release checks for OASIS4.3-GUI-0.1

Checked on 2026-10-17 on Intel macOS with Python 3.12.

| Check | Result |
|---|---|
| Offline unit and regression suite | 157 passed; 1 real-reference test deselected |
| Pinned real 1TQW reference integration | 1 passed; 328,518 reflections; prepared SHA-256 matched |
| GUI construction without a display server | Passed; product title, sequence, copies, solvent and help controls present |
| Conda environment update and GUI launch | Passed |
| One-command installer in a fresh prefix | Passed |
| Wheel and source distribution metadata (`twine check`) | Passed |
| `CITATION.cff` schema 1.2.0 validation | Passed |
| Public-tree audit | Passed; no MTZ/PDB/mmCIF/NPZ target data or personal absolute paths |
| GitHub Actions: Ubuntu, Windows, Intel macOS, Apple Silicon macOS | Passed ([run 39048389968](https://github.com/WeiDing535/OASIS4.3-GUI/actions/runs/39048389968)) |

The private 100-case benchmark is intentionally absent from this public repository. Historical aggregate benchmark observations are development evidence, not a blind external validation claim for this release.
