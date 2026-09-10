# Data provenance

The immutable inputs for the reference prototype are described in [`../config/sources.tsv`](../config/sources.tsv) and verified by [`../config/checksums.sha256`](../config/checksums.sha256). Run `scripts/fetch_reference.sh` from any directory to retrieve missing files. A checksum mismatch means that the remote representation changed and must be investigated rather than silently accepted.

`raw/` contains downloaded evidence and is excluded from Git. `processed/` contains deterministic tables and FASTA files produced only from CDS features. The published GenBank records are not edited in place.

For S3, `raw/s3/` contains the pinned MGnify catalogue, combined GFF/DNA files, primary BioSample/ENA XML, and captured pipeline definitions. Their URLs, sizes, retrieval timestamps, and hashes are owned by `config/s3_download_manifest.json` and `config/maize_*_downloads.json`. `interim/maize/` is a rebuildable concatenation of all annotated CDS and is excluded from Git. Derived candidate sequences and tables remain in `processed/` and `results/`.

Coordinates in the processed tables are zero-based, half-open coordinates within the submitted BGC record. These are not necessarily coordinates in the full source assembly. An accession without a version suffix reproduces the metadata embedded in the authors' export; it must not be silently upgraded.
