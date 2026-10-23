# Reference-data provenance

The default histogram reference is derived from CCP4's processed `reference-1tqw.mtz`. Provisioning downloads the pinned upstream file from:

`https://ftp.ccp4.ac.uk/ccp4/6.4.0/ccp4-6.4.0/lib/data/reference_structures/reference-1tqw.mtz`

The source SHA-256 is `e43c1c7ba791e33b3ce0b36a149425947e2ff2130cc46dece8a5328c51afd3c5`.

The program retains H, K, L and the two map-coefficient columns it uses, without changing coefficient values. The prepared-file SHA-256 is `1683f8b91212dde07c2f98cbc584d5649925e84ed7464b4c3c9b416c1ddd3ea7`.

The processed MTZ is not committed to this repository or included in wheels. `oasis-py reference-fetch` performs the explicit one-time HTTPS download, verifies the source hash, creates the minimal file and verifies its final hash. Upstream data terms apply separately from the software's LGPL licence.

The reference is an independently solved protein map and is outside the private 100-case development benchmark. Targets identical or closely related to 1TQW should use an independent reference. Nucleic-acid, membrane-rich or chemically unusual targets may require a separately validated reference supplied through `--reference-mtz`.

No user benchmark MTZ, complete target structure or unpublished result is included in the public repository.
