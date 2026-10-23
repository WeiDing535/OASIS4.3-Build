# OASIS4.3-GUI-0.1

English | [中文](README.md)

Project home: [github.com/WeiDing535/OASIS4.3-GUI](https://github.com/WeiDing535/OASIS4.3-GUI)

OASIS4.3-GUI-0.1 is independent research software for protein X-ray SAD phasing. It brings together a modern Python implementation of the OASIS/P+ method, continuous SAD phase likelihood, and density modification that does not depend on an external Parrot runtime, with both command-line and desktop interfaces. The program reads an experimental MTZ file, a protein sequence, and heavy-atom coordinates; estimates the asymmetric-unit composition and solvent content automatically; and writes an MTZ containing PHIB, FOM, and Hendrickson–Lattman coefficients.

The project is maintained by DingWei (丁玮), Institute of Physics, Chinese Academy of Sciences, and is licensed under LGPL-2.1-or-later. The algorithms remain at the research-validation stage. FOM is a model confidence estimate and cannot, by itself, establish that a structure has been solved.

## One-command installation

A 64-bit Python 3.12 installation, an internet connection, and approximately 1 GB of free space are required. From the source directory, run:

```bash
python3.12 install.py
```

On Windows, use:

```powershell
py -3.12 install.py
```

The installer creates an isolated environment in the user's home directory, installs the GUI dependencies, and downloads a generic reference MTZ once from a fixed HTTPS location with SHA-256 verification. For offline installation, first obtain the original 1TQW reference file identified in the data-source documentation, then run:

```bash
python3.12 install.py --reference-source /path/to/reference-1tqw.mtz
```

The default installation location is `~/.oasis-sad`. Use `--prefix` to select another location, `--no-launch` to install without opening the GUI, or `--cli-only` to install only the command-line application.

## Conda environment

The repository provides a reproducible Conda environment file. This installs from source; the project has not yet been published on conda-forge.

```bash
conda env create -f environment.yml
conda activate oasis43-gui
oasis-py reference-fetch
oasis-gui
```

## GUI inputs

The current GUI focuses on protein SAD. It requires:

- an MTZ containing anomalous amplitudes and their standard uncertainties;
- a protein sequence in FASTA, PIR, or plain-text format;
- the number of copies of the complete sequence file in the crystallographic asymmetric unit;
- a heavy-atom PDB or mmCIF file in the same unit cell, space group, and origin as the experimental data;
- a positive f″ value appropriate for the element and experimental wavelength.

The program calculates a Matthews solvent-content estimate from the sequence, copy number, unit cell, and space group. The `?` button beside each input provides detailed guidance. SIR, DMR, and additional experimental options are available through `oasis-py --help`.

## Testing

The standard test suite is fully offline and uses synthetic reference data:

```bash
python -m pip install -e '.[test,gui]' -c constraints.txt
pytest -m 'not reference_data'
```

GitHub Actions has passed on Linux, Windows, Intel macOS, and Apple Silicon macOS; see the [release-check record](docs/RELEASE_CHECKS.md). To run the real-reference integration test, first execute `oasis-py reference-fetch`, then explicitly run `pytest -m reference_data`.

## Citation and method background

Software citation metadata are provided in [CITATION.cff](CITATION.cff). Background on the OASIS/IPCAS methods includes: Ding W. et al., *IPCAS: a direct-method-based pipeline from phasing to model building and refinement for macromolecular structure determination*, J. Appl. Cryst. 53, 253–261 (2020), DOI: [10.1107/S1600576719015115](https://doi.org/10.1107/S1600576719015115).

For the Clipper-derived numerical kernel, corresponding source, and additional permission, see [Third-party code and licences](docs/THIRD_PARTY.md). The reference MTZ provenance, checksums, and distribution policy are documented in [Reference-data provenance](docs/DATA_SOURCES.md).

## Contact

DingWei (丁玮): dingwei@iphy.ac.cn
