# Third-party code and licences

## Clipper

`src/oasis_py/clipper_likelihood.py` adapts the `SFweight_spline<T>::TargetHL` numerical kernel from Clipper 2.1 (`clipper/contrib/sfweight.cpp`). Copyright © 2000–2006 Kevin Cowtan and the University of York. It is licensed under LGPL-2.1-or-later with the additional linking permission reproduced in `src/oasis_py/resources/THIRD_PARTY_NOTICES.md`.

The wheel and source distribution include the complete, unchanged corresponding archive `clipper-2.1.20201109.tar.gz`, its SHA-256 record and the LGPL text. No Clipper binary is needed at runtime.

## Qt and Python dependencies

The optional GUI uses PySide6-Essentials. NumPy, SciPy, Gemmi, PySide6 and their dependencies are installed from their package channels and retain their own licences. Their binaries are not copied into this repository.

## Parrot and historical OASIS sources

The runtime does not execute or distribute Parrot. Historical Parrot application source and the original Fortran programs are deliberately absent from this public repository because this release does not assert redistribution rights for them. The independent density-modification implementation is supplied as Python source here.

The project-level LGPL licence applies to this repository's code. Reference diffraction/map data have separate provenance and terms described in `DATA_SOURCES.md`.
