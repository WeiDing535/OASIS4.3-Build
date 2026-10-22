# OASIS4.3-Build

[中文](README.md) | English

Research software for building protein models from existing crystallographic phase maps. It includes the PhaseBuild 0.8.3 native engine, a local browser GUI, and explicit external Buccaneer/REFMAC integration. The earlier [OASIS4.3-GUI](https://github.com/WeiDing535/OASIS4.3-GUI) repository remains unchanged.

## Linux-first complete bundle

The next complete runtime target is Linux x86_64, with Python and scientific engines included. The Linux launcher, ELF/dependency audit and tar.gz assembler are implemented; an actual Linux runtime bundle still requires native runtime inputs and clean-system validation. See [Linux preparation](docs/LINUX.md). The existing Intel Mac local-evaluation bundle is retained. No public complete binary release has been made.

## Run from source

Use an existing Python ≥3.9 environment with Gemmi, NumPy and SciPy. External construction/refinement additionally requires a user-provided CCP4 installation and reference resources. No Python or CCP4 runtime is bundled.

On macOS with CCP4 at `/Applications/ccp4-9`, double-click `PhaseBuild.command` or run `./PhaseBuild.command`. For other installations set `PHASEBUILD_PYTHON` and `PHASEBUILD_RUNTIME` before launching. Alternatively run `python -B run_gui.py --project .` from the repository root.

Open the complete localhost URL printed by the launcher. Keep the service running; closing the page does not stop calculations. Stop jobs in the GUI and press Ctrl+C in the terminal to exit. Always use a new output directory.

## Scope and validation

Modes include native model building, model extension, external REFMAC refinement, standalone fixed-coordinate R factors, and seeded or unseeded external Buccaneer/REFMAC recycling with a fixed original experimental HL prior. This is downstream model building from supplied phases, not a complete upstream SAD phasing engine.

The source baseline passed 186 regression cases on macOS Intel with CCP4 Python 3.9.18. GUI loading and its local API were checked. External-engine contract tests include mocked processes. The new unseeded integration has no newly established real-data structure-solving result. Cross-platform runtime support is not established by this local check.

Run tests using an existing pytest installation: `python -B -m pytest -q -p no:cacheprovider tests`.

Maintained by DingWei (丁玮), Institute of Physics, Chinese Academy of Sciences. LGPL-2.1-or-later; see LICENSE. Third-party packages, engines and reference resources retain their own terms. See [provenance](docs/PROVENANCE.md).
