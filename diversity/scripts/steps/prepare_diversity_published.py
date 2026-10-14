#!/usr/bin/env python3
"""Freeze sequence/provenance comparators from the two published collections."""
import hashlib
import json
import zipfile
import warnings

import pandas as pd
from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.SeqRecord import SeqRecord
from prepare_literature_reference import parse_member
from profile_search import ROOT


def reference_metadata(report):
    info = report["assembly_info"]
    sample = info.get("biosample", {})
    attrs = {r["name"]: r.get("value", "") for r in sample.get("attributes", [])}
    checkm = report.get("checkm_info", {})
    return {"assembly": report["accession"], "organism": report["organism"]["organism_name"],
            "strain": attrs.get("strain", report["organism"].get("infraspecific_names", {}).get("strain", "")),
            "biosample": sample.get("accession", ""), "assembly_level": info["assembly_level"],
            "biosample_title": sample.get("description", {}).get("title", ""),
            "biosample_package": sample.get("package", ""), "metagenomic": attrs.get("metagenomic", ""),
            "sample_type": attrs.get("sample_type", ""), "isolation_source": attrs.get("isolation_source", ""),
            "host": attrs.get("host", ""), "country": attrs.get("geo_loc_name", ""),
            "metagenome_source": attrs.get("metagenome_source", ""),
            "wgs_project": report.get("wgs_info", {}).get("wgs_project_accession", ""),
            "source_completeness": checkm.get("completeness", ""),
            "source_contamination": checkm.get("contamination", "")}


def main():
    processed = ROOT / "data/processed"
    literature = pd.read_csv(processed / "literature_bgc_manifest.tsv", sep="\t", keep_default_na=False)
    original = pd.read_csv(processed / "reference_manifest.tsv", sep="\t", keep_default_na=False)
    records, metadata = [], []
    with zipfile.ZipFile(ROOT / "data/raw/research_2026-09-10/triculamin_lasso_BGCs.zip") as archive:
        for row in literature.to_dict("records"):
            record, _, _ = parse_member(row["archive_member"], archive.read(row["archive_member"]))
            records.append(SeqRecord(record.seq, id=row["bgc_id"], description=""))
            metadata.append({"record_id": row["bgc_id"], "source": "Merrild_2025",
             "accession": row["accession"], "strain": row["strain"],
             "source_group": row["source_group"], "source_file": row["archive_member"]})
    for row in original.to_dict("records"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            record = SeqIO.read(ROOT / "data/raw" / row["relative_raw_path"], "genbank")
        source = next((f for f in record.features if f.type == "source"), None)
        identifier = "lar29__" + row["sample_id"]
        records.append(SeqRecord(record.seq, id=identifier, description=""))
        metadata.append({"record_id": identifier, "source": "Jangra_LAR_BGC_29",
                         "accession": row["accession"], "strain": ";".join(source.qualifiers.get("strain", [])) if source else "",
                         "source_group": "canonical", "source_file": row["file_name"]})
    for row, record in zip(metadata, records):
        row["DNA_sha256"] = hashlib.sha256(str(record.seq).upper().encode()).hexdigest()
    SeqIO.write(records, processed / "discovery/diversity_published_DNA.fna", "fasta")
    pd.DataFrame(metadata).to_csv(processed / "discovery/diversity_published_records.tsv", sep="\t", index=False)
    c_records = list(SeqIO.parse(processed / "literature_unique_cyclases.faa", "fasta"))
    c_records += list(SeqIO.parse(processed / "lrcc_reference.faa", "fasta"))
    SeqIO.write(c_records, processed / "discovery/diversity_published_C.faa", "fasta")
    # Preserve primary metadata already present in the frozen source snapshot.
    source_metadata = []
    selected = set(pd.read_csv(ROOT / "inputs/selections/diversity_refseq_species.tsv", sep="\t").assembly)
    for path in sorted((ROOT / "data/raw/diversity_2026-10-12").glob("refseq_reference_page_*.json")):
        for report in json.loads(path.read_text())["reports"]:
            if report["accession"] not in selected:
                continue
            source_metadata.append(reference_metadata(report))
    if set(r["assembly"] for r in source_metadata) != selected:
        raise ValueError("Primary metadata do not cover the fixed RefSeq selection")
    pd.DataFrame(source_metadata).to_csv(processed / "discovery/diversity_refseq_metadata.tsv", sep="\t", index=False)
    print(f"Published comparator: {len(records)} exported records; {len(set(r['DNA_sha256'] for r in metadata))} exact forward DNA strings; records are not independent BGC counts")


if __name__ == "__main__":
    main()
