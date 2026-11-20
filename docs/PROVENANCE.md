# Source provenance

This public source export was prepared on 2026-10-22 from the locally verified PhaseBuild 0.8.3 trial copy. Numerical modules, GUI assets and tests were retained without algorithm changes. Repository documentation, portable-source packaging metadata and the source launcher were adapted for publication.

The small profile/conformer assets retain their original metadata. They were developed from the unrelated 1TQW reference, with separate historical validation metadata where available. Metadata may describe earlier local training locations; those paths are provenance, not required downloads. The repository does not bundle the original training structures or an external-engine runtime.

The source GUI is a loopback web application. Python packages (Gemmi, NumPy, SciPy) and optional external engines (CCP4 Buccaneer/REFMAC and their dictionaries/reference resources) remain separately supplied dependencies, subject to their own licensing. This source export is not a redistribution clearance for external runtimes.

The earlier OASIS4.3-GUI source and v0.1 release remain in their original repository. OASIS4.3-Build is published separately and retains PhaseBuild as its internal package identity.

The subsequent unified development integration adds an unchanged `oasis_py` numerical-source snapshot and new orchestration/UI adapters. The original source-copy equivalence statement above refers to the first publication baseline, not later interface additions. Upstream hashes and one explicitly documented test-only portability patch are in `vendor/oasis43`.

## Public example and validation data

`validation_ten_cases/inputs/` contains structure-factor, sequence and heavy-atom inputs derived from the public PDB entries 3O2E, 5D9D, 5HHK, 6HF7, 6ZPP, 7E1D, 7ENM, 7TBS, 7TFQ and 7TJ1, plus the insulin S-SAD teaching example in `ins/`. The wwPDB usage policy states that data files in the PDB archive are released under CC0 1.0 and encourages attribution of the original structure authors; this repository does not claim any new licence over those data. `models/` contains only small derived statistical parameters trained or validated on the public 1TQW/1AJR references; the original training structures are not redistributed. Reference-data provenance is recorded in `src/oasis_py/resources/reference-1tqw.json` and `vendor/oasis43/docs/DATA_SOURCES.md`.

Generated validation products under `validation_ten_cases/campaign_*` and `validation_ten_cases/verify_*`, local work outputs, and the third-party runtime are excluded from the public tree by `.gitignore`.
