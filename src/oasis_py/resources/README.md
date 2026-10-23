# Provisioned protein histogram reference

The public repository does not contain `reference-1tqw.mtz`. The explicit
`oasis-py reference-fetch` command retrieves the pinned upstream file and
creates a cache file containing exactly H, K, L and the two map-coefficient
columns read by the independent Python histogram algorithm. Values, cell,
space group and reflection order are unchanged from the reference used in
the preceding 100-case experiments. Only unused columns were removed.
`reference-1tqw.json` records the upstream-data URL and both file hashes.

Original reference data: CCP4 reference_structures/reference-1tqw.mtz:
https://ftp.ccp4.ac.uk/ccp4/6.4.0/ccp4-6.4.0/lib/data/reference_structures/reference-1tqw.mtz
Reference structure and its authors: https://www.rcsb.org/structure/1TQW
Reference-use documentation: https://www.ccp4.ac.uk/html/parrot.html

The resource is map data, not Parrot source or executable code. Running the
Python program does not require a CCP4/Parrot installation. PDB 1TQW is
outside our 100 target cases. No identity/homology screening against future
user targets is implemented: for 1TQW itself or closely related targets,
provide an independent reference with --reference-mtz.

Requirements for another reference: valid crystal cell and space group;
finite positive map amplitudes FC_BEST.F_phi.F (MTZ type F), meaningful map
phases FC_BEST.F_phi.phi (type P), and enough reflections in the working
resolution range. It need not share the target's cell or space group.
Use an independently solved, representative protein map. A raw observed
SAD MTZ without map phases is insufficient. References dominated by nucleic
acids, membrane lipids, very unusual chemistry or inadequate resolution
need separate validation. The bundled map extends to approximately 1.2 A;
operation above its resolution coverage has not been validated.

The LGPL Clipper kernel's notice, complete upstream source and licence are
also included here. No Clipper binary library is required at runtime. See
the repository documents `docs/DATA_SOURCES.md` and `docs/THIRD_PARTY.md`.
