#!/usr/bin/env python3
"""Inventory frozen published exports; annotations are evidence, not role calls."""
import csv
import hashlib
import io
import json
from pathlib import Path
import warnings
from zipfile import ZipFile

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "inputs/2026-10-13/base/data/raw"
OUTPUT = ROOT / "results/tables"
REPAIR_MEMBER = "noncanonical/WP_199826149.1.gbk"
BROKEN = (
    "            project_name        :: SeaBioTech actinobacteria experimental_factor\r\n"
    "            :: not applicable"
)
REPAIRED = (
    "            project_name        :: SeaBioTech actinobacteria\r\n"
    "            experimental_factor :: not applicable"
)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def joined(feature, key):
    return ";".join(feature.qualifiers.get(key, []))


def write_table(name, rows):
    with (OUTPUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def check_reference_anchor():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        genome = SeqIO.read(INPUTS / "CP169648.1.gb", "genbank")
        with ZipFile(INPUTS / "LAR-BGC-42280aa.zip") as archive:
            member = next(n for n in archive.namelist() if n.endswith("/Paenibacillus M2.gbk"))
            m2 = SeqIO.read(io.StringIO(archive.read(member).decode()), "genbank")
        with ZipFile(INPUTS / "research_2026-09-10/triculamin_lasso_BGCs.zip") as archive:
            dr312 = SeqIO.read(io.StringIO(archive.read("canonical/WP_222215840.1.gbk").decode()), "genbank")
    sequence, full = str(m2.seq).upper(), str(genome.seq).upper()
    reverse = str(m2.seq.reverse_complement()).upper()
    m2_cds = [f for f in m2.features if f.type == "CDS"]
    # The reference C is the fourth CDS in this published M2 export.
    dr312_c = next(f for f in dr312.features if f.type == "CDS" and "WP_222215840.1" in f.qualifiers.get("protein_id", []))
    result = {
        "reference_genome": genome.id,
        "reference_organism": genome.annotations.get("organism"),
        "m2_export_bp": len(sequence),
        "forward_exact_occurrences": full.count(sequence),
        "reverse_exact_occurrences": full.count(reverse),
        "forward_start_0based": full.find(sequence),
        "reverse_start_0based": full.find(reverse),
        "merrild_export_organism": dr312.annotations.get("organism"),
        "merrild_C_accession": "WP_222215840.1",
        "m2_C_vs_merrild_C_exact_protein": m2_cds[3].qualifiers["translation"] == dr312_c.qualifiers["translation"],
        "m2_export_exact_in_merrild_window": sequence in str(dr312.seq).upper(),
        "limitation": "Sequence correspondence only. Protein identity does not establish carrier identity or HGT. Raw reads not inspected.",
    }
    (OUTPUT / "reference_anchor_checks.json").write_text(json.dumps(result, indent=2) + "\n")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    loci, genes = [], []
    for source, path in [
        ("Nature_2025", INPUTS / "LAR-BGC-42280aa.zip"),
        ("Merrild_2025", INPUTS / "research_2026-09-10/triculamin_lasso_BGCs.zip"),
    ]:
        with ZipFile(path) as archive:
            for member in sorted(archive.namelist()):
                if not member.endswith((".gb", ".gbk")):
                    continue
                raw = archive.read(member)
                text, repair = raw.decode("utf-8"), ""
                if source == "Merrild_2025" and member == REPAIR_MEMBER:
                    if text.count(BROKEN) != 1:
                        raise ValueError("Known malformed comment changed")
                    text = text.replace(BROKEN, REPAIRED, 1)
                    repair = "split_MIGS_comment_fields_only"
                    derived = ROOT / "inputs/2026-10-13/derived" / member
                    derived.parent.mkdir(parents=True, exist_ok=True)
                    derived.write_bytes(text.encode("utf-8"))
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    records = list(SeqIO.parse(io.StringIO(text), "genbank"))
                if len(records) != 1:
                    raise ValueError(f"Expected one record: {member}")
                record = records[0]
                locus_id = source + "_" + sha(member.encode())[:12]
                dna = str(record.seq).upper()
                cds = [f for f in record.features if f.type == "CDS"]
                source_features = [f for f in record.features if f.type == "source"]
                source_notes = ";".join(joined(f, "note") for f in source_features)
                # Hashes group exact content; shorter exports may overlap larger ones.
                rc = str(record.seq.reverse_complement()).upper()
                loci.append({
                    "locus_id": locus_id, "source": source,
                    "author_group": "canonical" if source == "Nature_2025" else member.split("/")[0],
                    "archive_member": member, "raw_sha256": sha(raw),
                    "parsed_sha256": sha(text.encode()), "parse_repair": repair,
                    "parser_warning_count": len(caught),
                    "record_id": record.id,
                    "accessions": ";".join(record.annotations.get("accessions", [])),
                    "organism_in_export": record.annotations.get("organism", ""),
                    "strain_in_source": ";".join(joined(f, "strain") for f in source_features),
                    "taxon_xrefs": ";".join(joined(f, "db_xref") for f in source_features),
                    "source_notes": source_notes,
                    "record_bp": len(dna), "ambiguous_bases": sum(b not in "ACGT" for b in dna),
                    "dna_sha256": sha(dna.encode()),
                    "strand_independent_dna_sha256": sha(min(dna, rc).encode()),
                    "cds_count": len(cds),
                    "translated_cds_count": sum(bool(f.qualifiers.get("translation")) for f in cds),
                    "assembly_mapping_status": "not_checked_in_this_inventory",
                    "assembly_quality_status": "requires_full_assembly",
                })
                for index, feature in enumerate(record.features):
                    if feature.type != "CDS":
                        continue
                    protein = joined(feature, "translation").replace(" ", "").rstrip("*")
                    coding = feature.extract(record.seq)
                    table = int(feature.qualifiers.get("transl_table", [11])[0])
                    offset = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        translated = str(coding[offset:].translate(table=table)).rstrip("*")
                    match = "not_annotated"
                    if protein:
                        match = "exact" if protein == translated else "mismatch_requires_review"
                        if match != "exact" and protein.startswith("M") and protein[1:] == translated[1:]:
                            match = "start_residue_only_requires_review"
                    genes.append({
                        "locus_id": locus_id, "feature_index": index,
                        "location": str(feature.location), "strand": feature.location.strand,
                        "gene": joined(feature, "gene"), "locus_tag": joined(feature, "locus_tag"),
                        "protein_id": joined(feature, "protein_id"), "label": joined(feature, "label"),
                        "product": joined(feature, "product"),
                        "protein_aa": len(protein), "protein_sha256": sha(protein.encode()) if protein else "",
                        "internal_stop_count": protein.count("*"),
                        "translation_vs_export_dna": match,
                        "role_validation_status": "not_assigned_by_inventory",
                    })
    write_table("published_loci.tsv", loci)
    write_table("published_cds.tsv", genes)
    groups = {}
    for row in loci:
        groups.setdefault(row["strand_independent_dna_sha256"], []).append(row)
    duplicates = [
        {"dna_group_sha256": key, "export_count": len(rows),
         "locus_ids": ";".join(r["locus_id"] for r in rows),
         "sources": ";".join(sorted(set(r["source"] for r in rows)))}
        for key, rows in groups.items() if len(rows) > 1
    ]
    if duplicates:
        write_table("exact_export_duplicates.tsv", duplicates)
    else:
        (OUTPUT / "exact_export_duplicates.tsv").write_text("dna_group_sha256\texport_count\tlocus_ids\tsources\n")
    summary = {
        "exports": len(loci), "cds_features": len(genes),
        "author_groups": {g: sum(r["author_group"] == g for r in loci) for g in sorted(set(r["author_group"] for r in loci))},
        "parse_repairs": sum(bool(r["parse_repair"]) for r in loci),
        "parser_warnings": sum(r["parser_warning_count"] for r in loci),
        "distinct_exact_export_dna": len(groups), "duplicate_export_groups": len(duplicates),
        "protein_dna_checks": {s: sum(r["translation_vs_export_dna"] == s for r in genes) for s in sorted(set(r["translation_vs_export_dna"] for r in genes))},
        "limitation": "Export inventory only. Neither unique BGC counts, carrier identity, assembly quality nor HGT established.",
    }
    (OUTPUT / "published_input_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    check_reference_anchor()
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
