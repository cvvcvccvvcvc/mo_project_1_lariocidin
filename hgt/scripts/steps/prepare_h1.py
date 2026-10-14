#!/usr/bin/env python3
"""Prepare sequence/provenance inventories from frozen publication archives."""
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import warnings
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from audit_published_inputs import BROKEN, REPAIRED, REPAIR_MEMBER

ROOT = Path(__file__).resolve().parents[2]
SNAP = ROOT / "inputs/2026-10-13/base/data/raw"
H1 = ROOT / "inputs/H1_2026-10-13"
TABLES = ROOT / "results/tables"
SEQUENCES = ROOT / "results/H1_sequences"


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def published_records():
    for source, filename in [
        ("Nature_2025", SNAP / "LAR-BGC-42280aa.zip"),
        ("Merrild_2025", SNAP / "research_2026-09-10/triculamin_lasso_BGCs.zip"),
    ]:
        with ZipFile(filename) as archive:
            for member in sorted(archive.namelist()):
                if not member.endswith((".gb", ".gbk")):
                    continue
                text = archive.read(member).decode()
                if source == "Merrild_2025" and member == REPAIR_MEMBER:
                    if text.count(BROKEN) != 1:
                        raise ValueError("Known repair target changed")
                    text = text.replace(BROKEN, REPAIRED, 1)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    record = SeqIO.read(io.StringIO(text), "genbank")
                locus_id = source + "_" + hashlib.sha256(member.encode()).hexdigest()[:12]
                yield locus_id, source, member, record


def author_role(feature):
    labels = feature.qualifiers.get("label", [])
    products = feature.qualifiers.get("product", [])
    if any(p.startswith("C_") for p in products) or "C" in labels:
        return "C"
    for role in ["B1", "B2", "A1", "A", "D", "E", "F", "T", "P"]:
        if any(label == role or label.startswith(role + "_") for label in labels):
            return role
    return ""


def wgs_prefix(accession):
    accession = accession.removeprefix("NZ_")
    match = re.match(r"^([A-Z]{4}|[A-Z]{6})\d{2}\d+", accession)
    return match.group(1) if match else ""


def main():
    SEQUENCES.mkdir(parents=True, exist_ok=True)
    inventory = {r["locus_id"]: r for r in read_tsv(TABLES / "published_loci.tsv")}
    nature_metadata = {r["file_name"]: r for r in read_tsv(H1 / "reference_metadata/reference_metadata.tsv")}
    nature_manifest = {r["file_name"]: r for r in read_tsv(H1 / "reference_metadata/reference_manifest.tsv")}
    loci, genes = [], []
    locus_dna, coding_dna, proteins, c_proteins, nature_proteins = [], [], [], [], []
    for locus_id, source, member, record in published_records():
        group = inventory[locus_id]["author_group"]
        if group not in ["canonical", "noncanonical"]:
            continue
        accession = inventory[locus_id]["accessions"]
        name_evidence = record.annotations.get("organism", "")
        sample_id = ""
        if source == "Nature_2025":
            metadata = nature_metadata[Path(member).name]
            name_evidence, sample_id = metadata["display_name"], metadata["sample_id"]
            accession = metadata["accession_override"] or nature_manifest[Path(member).name]["accession"] or accession
        if member == "canonical/Ptar_C.gbk":
            name_evidence = "Paenibacillus tarimensis DSM 19409"
        elif member == "canonical/Bgel_C.gbk":
            name_evidence = "Brevibacillus gelatini DSM 100115"
        elif member == "canonical/Cpal_C.gbk":
            name_evidence = "Chitinasiproducens palmae JS23"
        # These anchor labels are publication evidence, not carrier validation.
        source_feature = next((f for f in record.features if f.type == "source"), None)
        strain = ";".join(source_feature.qualifiers.get("strain", [])) if source_feature else ""
        cds = [f for f in record.features if f.type == "CDS"]
        c_count = sum(author_role(f) == "C" for f in cds) if source == "Merrild_2025" else 0
        precursor_count = sum(author_role(f) in ["A", "A1"] for f in cds) if source == "Merrild_2025" else 0
        loci.append({
            "locus_id": locus_id, "source": source, "lineage": group,
            "archive_member": member, "sample_id": sample_id,
            "accession_evidence": accession, "wgs_prefix": wgs_prefix(accession),
            "name_evidence": name_evidence, "strain_evidence": strain,
            "author_c_features": c_count, "author_precursor_features": precursor_count,
            "record_bp": len(record), "dna_sha256": inventory[locus_id]["dna_sha256"],
        })
        locus_dna.append(SeqRecord(record.seq, id=locus_id, description=""))
        for index, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            gene_id = f"{locus_id}.f{index}"
            protein = "".join(feature.qualifiers.get("translation", [])).replace(" ", "").rstrip("*")
            role = author_role(feature) if source == "Merrild_2025" else ""
            gene = {
                "gene_id": gene_id, "locus_id": locus_id, "feature_index": index,
                "lineage": group, "author_role": role, "location": str(feature.location),
                "start_0based": int(feature.location.start), "end_0based": int(feature.location.end),
                "strand": feature.location.strand,
                "parts": json.dumps([[int(p.start), int(p.end), p.strand] for p in feature.location.parts]),
                "partial_location": any(type(p).__name__ != "ExactPosition" for p in [feature.location.start, feature.location.end]),
                "protein_id": ";".join(feature.qualifiers.get("protein_id", [])),
                "locus_tag": ";".join(feature.qualifiers.get("locus_tag", [])),
                "label": ";".join(feature.qualifiers.get("label", [])),
                "product": ";".join(feature.qualifiers.get("product", [])),
                "protein_aa": len(protein), "protein_sha256": hashlib.sha256(protein.encode()).hexdigest() if protein else "",
                "codon_start": int(feature.qualifiers.get("codon_start", [1])[0]),
                "translation_table": int(feature.qualifiers.get("transl_table", [11])[0]),
            }
            genes.append(gene)
            coding_dna.append(SeqRecord(feature.extract(record.seq), id=gene_id, description=""))
            if protein:
                sequence = SeqRecord(Seq(protein), id=gene_id, description="")
                proteins.append(sequence)
                if role == "C":
                    c_proteins.append(sequence)
                if source == "Nature_2025":
                    nature_proteins.append(sequence)
    write_tsv(TABLES / "H1_publication_loci.tsv", loci)
    write_tsv(TABLES / "H1_publication_genes.tsv", genes)
    for name, rows in [
        ("loci.fna", locus_dna), ("coding_sequences.fna", coding_dna),
        ("proteins.faa", proteins), ("merrild_author_C.faa", c_proteins),
        ("nature_proteins.faa", nature_proteins),
    ]:
        SeqIO.write(rows, SEQUENCES / name, "fasta")
    prefixes = {r["wgs_prefix"] for r in loci} - {""}
    names = {r["name_evidence"].lower() for r in loci}
    genera = {s.split()[0] for s in names if s and s.split()[0] not in ["uncultured", "synthetic", "unspecified"]}
    candidates, relatives, npdc_candidates = [], [], []
    npdc_codes = {x for locus in loci for x in re.findall(r"NPDC\d+", locus["archive_member"])}
    path = H1 / "reference_metadata/refseq_bacterial_assembly_summary.txt"
    with path.open() as handle:
        next(handle)
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["assembly_accession"] = row.pop("#assembly_accession")
            name = row["organism_name"].lower()
            prefix = re.match(r"^[A-Z]+", row["wgs_master"])
            if prefix and prefix.group() in prefixes:
                candidates.append(row)
            if name in names or name.split()[0] in genera:
                relatives.append(row)
            if npdc_codes.intersection(re.findall(r"NPDC\d+", ";".join(row[k] for k in ["infraspecific_name", "isolate", "asm_name"]))):
                npdc_candidates.append({"assembly_accession": row["assembly_accession"],
                                        **{k: v for k, v in row.items() if k != "assembly_accession"}})
    if candidates:
        write_tsv(TABLES / "H1_wgs_assembly_candidates.tsv", candidates)
    if relatives:
        write_tsv(TABLES / "H1_relative_metadata_pool.tsv", relatives)
    if npdc_candidates:
        write_tsv(TABLES / "H1_npdc_assembly_candidates.tsv", npdc_candidates)
    print(json.dumps({"target_exports": len(loci), "genes": len(genes), "author_C_features": len(c_proteins),
                      "wgs_candidate_rows": len(candidates), "relative_metadata_rows": len(relatives),
                      "missing_C_annotations": sum(r["source"] == "Merrild_2025" and r["author_c_features"] != 1 for r in loci)}))


if __name__ == "__main__":
    main()
