#!/usr/bin/env python3
"""Retrieve pinned combined GFF/DNA files and read their annotated CDS."""
import gzip
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable

from profile_search import ROOT


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fetch_one(row, old):
    identifier = row["Genome"]
    url = row["FTP_download"].replace("ftp://", "https://", 1)
    path = ROOT / f"data/raw/s3/genomes/{identifier}.gff.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = path.with_suffix(".gz.download")
        subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
                        "--retry", "3", "--max-time", "180", url, "--output", str(temporary)], check=True)
        if old and digest(temporary) != old["sha256"]:
            raise ValueError(f"Remote snapshot changed: {identifier}")
        temporary.rename(path)
    checksum = digest(path)
    if old and (checksum != old["sha256"] or url != old["url"]):
        raise ValueError(f"Pinned source mismatch: {identifier}")
    return identifier, old or {"url": url, "local_path": str(path.relative_to(ROOT)),
                              "bytes": path.stat().st_size, "sha256": checksum,
                              "retrieved_utc": datetime.now(timezone.utc).isoformat()}


def fetch(rows):
    manifest = ROOT / "inputs/provenance/maize_genome_downloads.json"
    records = json.loads(manifest.read_text()) if manifest.exists() else {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(fetch_one, row, records.get(row["Genome"])) for row in rows]
        for index, future in enumerate(as_completed(futures), 1):
            identifier, record = future.result()
            records[identifier] = record
            manifest.write_text(json.dumps(dict(sorted(records.items())), indent=2) + "\n")
            if index % 50 == 0 or index == len(futures):
                print(f"GFF cached/verified: {index}/{len(futures)}", flush=True)


def read_genome(identifier):
    features = []
    with gzip.open(ROOT / f"data/raw/s3/genomes/{identifier}.gff.gz", "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip().split("\t")
            if f[2] != "CDS":
                continue
            attributes = dict(item.split("=", 1) for item in f[8].split(";") if "=" in item)
            features.append({"Genome": identifier, "target_id": attributes["ID"], "contig": f[0],
                             "start_0based": int(f[3]) - 1, "end_0based_exclusive": int(f[4]),
                             "strand": f[6], "phase": int(f[7]), "source": f[1],
                             "product": unquote(attributes.get("product", "")),
                             "attributes": attributes})
        sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    assert sequences, f"No embedded DNA: {identifier}"
    for feature in features:
        dna = sequences[feature["contig"]].seq[feature["start_0based"]:feature["end_0based_exclusive"]]
        if feature["strand"] == "-":
            dna = dna.reverse_complement()
        dna = dna[feature["phase"]:]
        table = int(feature["attributes"].get("transl_table", 11))
        assert len(dna) % 3 == 0, f"Incomplete codon: {feature['target_id']}"
        protein = str(dna.translate(table=table)).removesuffix("*")
        # Prokka/Prodigal initiator translation; validate representatives against
        # the released protein catalogue before trusting the genome-wide scan.
        if str(dna[:3]) in CodonTable.unambiguous_dna_by_id[table].start_codons:
            protein = "M" + protein[1:]
        feature["protein"] = protein
        feature["contig_bp"] = len(sequences[feature["contig"]])
    return features, sequences


if __name__ == "__main__":
    table = pd.read_csv(ROOT / "inputs/selections/maize_s3_genomes.tsv", sep="\t")
    fetch(table[table.in_scope].to_dict("records"))
