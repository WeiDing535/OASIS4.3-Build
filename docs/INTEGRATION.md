# Unified SAD to model-building development version

The browser UI now includes **SAD 定相到模型**. The independent phasing and building entry points remain available. The upstream `oasis_py` numerical source is preserved without edits; source hashes, revision and license materials are recorded in `vendor/oasis43/UPSTREAM.json` and its documentation. The optional original Qt GUI is not required by the unified browser UI.

## Environment

Use the project conda environment (`environment.yml`): Python 3.12, NumPy 1.26, SciPy 1.13, Gemmi 0.7.5. A build operator can create it with `conda env create -f environment.yml`, activate `oasis43-unified`, and run `python -B run_gui.py --project .`. The local development copy also discovers its project `.conda-unified` environment through `PhaseBuild.command`.

New complete binary bundles must contain the unified Python 3.12 runtime and pass its import check. The previously uploaded Intel Mac 0.8.3 draft remains an earlier model-building package with Python 3.9; it is not silently replaced or described as the new integrated package. Linux remains the primary next binary target, pending supplied Linux runtime inputs.

## First integration scope

Inputs are an existing heavy-atom model in the correct hand/setting, SAD amplitudes with anomalous differences or complete Friedel pairs, original merged native observations with `F/SIGF/FreeR_flag`, standard-amino-acid FASTA/PIR, ASU copy count, explicit element/f-double-prime and an unrelated density-modification reference. Initial integration uses the existing continuous SAD likelihood and MLHL-best density-modification path. Heavy-atom search, automatic hand selection and automatic SeMet composition conversion are not added to this first combined entry point. Their independent workflows remain available.

Native/SAD indices must already be merged reciprocal-ASU HKLs in the same setting. The adapter refuses implicit reindexing and conflicting original free flags. Native observations are never overwritten or rescaled by this integration layer. Native may be the same observation source as SAD after providing the required observation columns, but the original native input must not already contain downstream map columns.

## Separate immutable artifacts

1. `01_inputs/sad.mtz`: original SAD columns plus authoritative native free flags joined by exact HKL. Original free values are retained; unmatched SAD rows and free reflections are excluded by upstream phasing.
2. `02_sad/experimental.mtz`: pre-density-modification SAD result with HLA/HLB/HLC/HLD. This file alone is the fixed experimental prior for later refinement.
3. `03_dm/map_phases.mtz`: final density-modified phases and FOM, with their own provenance. Its HL is not substituted for the original prior.
4. `04_bridge/building.mtz`: original native rows/columns and unchanged flags, plus `FWT = native F * final DM FOM` once and `PHWT = final DM PHIB`. Missing phases and initial free-reflection map coefficients remain NaN. This is the explicit initial mFo-style coefficient convention, not a claim to reproduce every upstream DM feedback coefficient.
5. `05_build/`: existing no-seed Buccaneer construction and fixed-original-HL REFMAC recycling. The last planned round is the endpoint. The final model is copied to `primary_model.cif` after its engine report hash is checked.

The complete sequence file is repeated by ASU copy count for building, preserving its distinct records rather than concatenating chains. No target-reference coordinates enter this workflow.

## CLI and restart

`python -B run_sad_to_model.py --help` lists the required inputs. Add `--check` for read-only validation and a plan; this creates no output directory. A new run always requires a new output directory.

After an interrupted run, repeat the same command with `--resume`. The controller checks configuration, source and dependency versions, input/runtime identities and all completed-stage artifact hashes. It does not repeat verified completed SAD or DM stages. Partial stage directories are renamed and preserved before rerunning that stage. An output lock prevents concurrent runners. GUI execution exposes five-stage progress and preserves failed outputs; recovery is currently through the explicit CLI `--resume`, not a GUI restart button.

## Scientific and distribution limits

Upstream SAD and DM exclude the original free reflections. Later model-derived maps from external engines have not been certified work-only, and Rfree is not an independent end-to-end blind score. This software integration does not establish real-data structure-solving performance. Component tests, actual synthetic SAD/DM checks and mock-builder orchestration checks are reported separately. Complete public binary redistribution still requires the external-engine and resource permissions already documented for standalone packages.

## Integration verification on 2026-10-23

Combined local conda regression: 366 passed, no failures/errors; one explicitly reference-data-dependent case was deselected. The suite includes actual OASIS SAD and DM on synthetic observations feeding a test construction substitute, plus exact-HKL flag/map tests, resume/tamper/resource checks, origin-choice rejection and real process-group cancellation cleanup. External model-building/refinement algorithms retain their existing component tests; no new real-data end-to-end success is claimed.

A real browser preflight with synthetic inputs and original free value 7 passed without starting calculations or creating the requested output directory. Upstream source/resource hashes (31 files) match the preserved OASIS commit. One downstream test-only numerical reproducibility adjustment is documented in `vendor/oasis43/TEST_PATCHES.md`; numerical engine files are unchanged.
