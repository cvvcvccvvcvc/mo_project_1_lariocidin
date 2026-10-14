#!/usr/bin/env python3
"""Resolve exact selected RefSeq FTP paths from official assembly metadata.

Assembly display names are not FTP directory names. Preserve the selected
accession versions and record this additional metadata snapshot separately.
"""
import csv
import hashlib
import json
from datetime import datetime, timezone
from profile_search import ROOT


def main():
    path = ROOT / "data/raw/diversity_2026-10-12/refseq_bacterial_assembly_summary.txt"
    selected_path = ROOT / "inputs/selections/species/diversity_refseq_species.tsv"
    with selected_path.open() as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    wanted = {r["assembly"] for r in selected}
    resolved = {}
    with path.open() as handle:
        for line in handle:
            if line.lstrip("# ").startswith("assembly_accession\t"):
                fields = line.lstrip("# ").rstrip().split("\t")
                break
        else:
            raise ValueError("No NCBI assembly-summary header")
        for values in csv.reader(handle, delimiter="\t"):
            if values and values[0] in wanted:
                row = dict(zip(fields, values))
                resolved[values[0]] = row["ftp_path"].replace("ftp://", "https://", 1).rstrip("/")
    missing = wanted - set(resolved)
    if missing:
        raise ValueError(f"Selected assembly versions absent from FTP metadata: {sorted(missing)[:10]}")
    for row in selected:
        row["ftp_path"] = resolved[row["assembly"]]
        if not row["ftp_path"].startswith("https://ftp.ncbi.nlm.nih.gov/genomes/all/"):
            raise ValueError(f"No official FTP path: {row['assembly']}")
    temporary = selected_path.with_suffix(".tsv.tmp")
    with temporary.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)
    temporary.replace(selected_path)
    with path.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    record = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
              "local_path": str(path.relative_to(ROOT)), "sha256": sha, "bytes": path.stat().st_size,
              "retrieved_utc": datetime.now(timezone.utc).isoformat(), "selected_accessions_resolved": len(selected)}
    (ROOT / "inputs/provenance/diversity_ftp_download_manifest.json").write_text(json.dumps({"refseq_assembly_summary": record}, indent=2)+"\n")
    print(f"Official FTP paths resolved for {len(selected)} fixed assembly versions")


if __name__ == "__main__":
    main()
