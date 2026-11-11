# Product and component versions

The focused OASIS4.3 product version is `4.3.0.dev7`, declared once in
`pyproject.toml` under `[tool.oasis43]`. `package_oasis.py` uses that value in
the focused archive's project metadata and manifest, and copies Python and
dependency requirements from the repository project metadata.

The development repository still packages the historical `phasebuild` library
as `0.8.3`; its legacy source and standalone package names retain that version.
The unchanged upstream `oasis_py` snapshot reports version `0.1`. These identify
components, not the focused product release. They must not be rewritten merely
to make their strings agree: the upstream hashes and historical release
identities are retained for reproducibility.

The focused source ZIP is a runnable product subset. It intentionally excludes
the repository test suite, Git history, local validation records and external
engines. Preserve the full repository and a separate backup for development.
A complete binary bundle for the focused product has not been validated or
released. The existing Intel Mac 0.8.3 bundle remains historical.

Version 4.3.0.dev2 separates `02_sad/observations.mtz` (fixed observation-only refinement prior) from `experimental.mtz` (unchanged phase result consumed by DM). This changes the focused refinement protocol; old checkpoints require their original code and must not be resumed under this version.

Version 4.3.0.dev3 adds the explicit `pplus_trust` GUI route. It imports original Phaser SAD observation HL, runs nonzero KL-bounded P+ branch updates and separate map-only RESOLVE/Wang14 plus classic40 DM, then builds from no protein model. Optional single-copy cleanup is followed by a fixed amplitude-only endpoint. Old method IDs remain available. The supplied prior must come from the same work data, hand and heavy model, before DM/model improvement; the GUI does not automatically run Phaser. New code requires new run directories; original source identities remain necessary for resuming old checkpoints.

Version 4.3.0.dev4 adds pplus_internal with internally generated observation likelihood and seeded internal MLHL DM. The local installation-free bundle includes scientific Python and authorized local Buccaneer/REFMAC components, forbids Phenix execution and does not require separately installed CCP4. This changes the scientific route; the dev3 external-branch Rfree claim is not inherited. Public redistribution conditions for bundled engines remain unresolved.
