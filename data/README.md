# Data provenance

The immutable inputs for the reference prototype are described in [`../config/sources.tsv`](../config/sources.tsv) and verified by [`../config/checksums.sha256`](../config/checksums.sha256). Run `scripts/fetch_reference.sh` from any directory to retrieve missing files. A checksum mismatch means that the remote representation changed and must be investigated rather than silently accepted.

`raw/` contains downloaded evidence and is excluded from Git. `processed/` contains deterministic tables and FASTA files produced only from CDS features. The published GenBank records are not edited in place.

Coordinates in the processed tables are zero-based, half-open coordinates within the submitted BGC record. These are not necessarily coordinates in the full source assembly. An accession without a version suffix reproduces the metadata embedded in the authors' export; it must not be silently upgraded.

