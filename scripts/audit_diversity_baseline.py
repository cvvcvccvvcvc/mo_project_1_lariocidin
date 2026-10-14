#!/usr/bin/env python3
"""Audit S3 novelty and evidence without modifying its original classifications."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

from fetch_hgt_assemblies import Cache, ROOT
from maize_genomes import read_genome
from prepare_literature_reference import parse_member


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, data):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


def wgs_project(accession):
    # Compare WGS master/version prefixes, not organism display names.
    match = re.match(r"(?:NZ_)?([A-Z]{4,6}\d{2})\d+", accession)
    return match.group(1) if match else ""


def main():
    processed, tables = ROOT / "data/processed", ROOT / "results/tables"
    raw = ROOT / "data/raw/diversity_2026-10-12"
    library = rows(processed / "literature_bgc_manifest.tsv")
    original = rows(processed / "reference_manifest.tsv")
    candidates = rows(tables / "maize_curated_candidates.tsv")
    domain_rows = {r["protein_id"]: r for r in rows(tables / "diversity_domain_evidence.tsv")}
    cores = {r["candidate_id"]: r for r in rows(tables / "maize_supported_core_comparison.tsv")}
    references = []
    with zipfile.ZipFile(ROOT / "data/raw/research_2026-09-10/triculamin_lasso_BGCs.zip") as archive:
        for row in library:
            record, _, _ = parse_member(row["archive_member"], archive.read(row["archive_member"]))
            references.append((row["archive_member"], row["accession"], row["strain"], str(record.seq).upper()))
    for row in original:
        record = SeqIO.read(ROOT / "data/raw" / row["relative_raw_path"], "genbank")
        source = next((f for f in record.features if f.type == "source"), None)
        strain = ";".join(source.qualifiers.get("strain", [])) if source else ""
        references.append(("LAR-BGC/" + row["file_name"], row["accession"], strain, str(record.seq).upper()))

    novelty, public_proteins = [], []
    candidate_sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "maize_direct_candidates.faa", "fasta")}
    cache = Cache(raw, ROOT / "config/diversity_novelty_download_manifest.json")
    for candidate in candidates:
        identifier = candidate["target_id"]
        if identifier not in cores:
            continue
        acc = candidate["Genome_accession"] + ".1"
        payload = json.loads((raw / f"{acc}__assembly_summary.json").read_text())["result"]
        assembly = payload[payload["uids"][0]]
        if candidate["Sample_accession"] != assembly["biosampleaccn"]:
            raise ValueError(f"BioSample mismatch for {identifier}")
        with gzip.open(raw / f"{acc}__refseq_proteins.faa.gz", "rt") as handle:
            protein_accessions = [r.id for r in SeqIO.parse(handle, "fasta")
                                  if str(r.seq) == candidate_sequences[identifier]]
        for protein in protein_accessions:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=" + protein + "&rettype=ipg&retmode=text"
            ipg = cache.download(protein + "__ipg", url, ".tsv")
            with ipg.open() as handle:
                for occurrence in csv.DictReader(handle, delimiter="\t"):
                    public_proteins.append({"candidate_id": identifier, "exact_RefSeq_protein": protein, **occurrence})
        features, contigs = read_genome(candidate["Genome"])
        contig_digests = sorted(hashlib.sha256(min(str(r.seq).upper(), str(r.seq.reverse_complement()).upper()).encode()).hexdigest()
                                for r in contigs.values())
        c_feature = next(f for f in features if f["target_id"] == identifier)
        core = cores[identifier]
        sequence = contigs[c_feature["contig"]].seq[int(core["core_start_0based"]):int(core["core_end_0based_exclusive"])]
        if c_feature["strand"] == "-":
            sequence = sequence.reverse_complement()
        forward, reverse = str(sequence).upper(), str(sequence.reverse_complement()).upper()
        matches = [r[0] for r in references if forward in r[3] or reverse in r[3]]
        project_matches = [r[0] for r in references if wgs_project(r[1]) == assembly["wgs"]]
        strain_matches = [r[0] for r in references if candidate["strain"] in r[2].split(";")]
        novelty.append({"candidate_id": identifier, "strain": candidate["strain"],
                        "genbank_assembly": assembly["synonym"]["genbank"],
                        "refseq_assembly": assembly["synonym"]["refseq"],
                        "paired_assembly_similarity": assembly["synonym"].get("similarity", ""),
                        "exact_C_protein_accessions_in_RefSeq": ";".join(protein_accessions),
                        "biosample": assembly["biosampleaccn"], "wgs_project": assembly["wgs"],
                        "strain_matches_in_422_records": ";".join(strain_matches),
                        "wgs_matches_in_422_records": ";".join(project_matches),
                        "exact_core_matches_in_422_records": ";".join(matches),
                        "core_sha256": hashlib.sha256(forward.encode()).hexdigest(),
                        "assembly_contig_multiset_sha256": hashlib.sha256("\n".join(contig_digests).encode()).hexdigest(),
                        "C_sha256": core["C_sha256"], "C_changes_vs_nearest": core["C_changes_vs_nearest"],
                        "novelty_scope": "comparison_to_29_LAR_and_393_Merrild_records_only"})
    write_table(tables / "diversity_bacillus_novelty.tsv", novelty)
    if public_proteins:
        write_table(tables / "diversity_bacillus_public_protein_occurrences.tsv", public_proteins)

    features = {r["feature_id"]: r for r in rows(processed / "all_cds.tsv")}
    lit_features = rows(processed / "literature_all_cds.tsv")
    lit_sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "literature_all_cds.faa", "fasta")}
    exact = defaultdict(list)
    for f in lit_features:
        if f["feature_id"] in lit_sequences:
            exact[lit_sequences[f["feature_id"]]].append(f)
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "all_cds.faa", "fasta")}
    evidence = []
    for r in rows(tables / "pseudomonadota_core_roles.tsv"):
        f = features[r["feature_id"]]
        matches = exact[sequences[r["feature_id"]]]
        evidence.append({"cohort": "historical_HGT_14", "locus": r["sample_id"], "protein_id": r["feature_id"],
                         "assigned_role": r["role"], "original_assignment": r["role_assignment_basis"],
                         "annotation_product": f["product"], "annotation_label": f["label"],
                         "existing_domain_annotations": f["sec_met_domain"],
                         "exact_published_CDS_matches": ";".join(x["feature_id"] for x in matches),
                         "exact_match_labels": ";".join(sorted({x["label"] for x in matches if x["label"]})),
                         "exact_match_products": ";".join(sorted({x["product"] for x in matches if x["product"]})),
                         "function_status": "inferred_role_not_per_locus_experimental_validation"})
    neighbors = {r["target_id"]: r for r in rows(tables / "maize_candidate_neighbors.tsv")}
    for r in rows(tables / "maize_supported_core_proteins.tsv"):
        f = neighbors[r["protein_id"]]
        matches = exact[r["sequence"]]
        evidence.append({"cohort": "historical_Maize_3", "locus": r["candidate_id"], "protein_id": r["protein_id"],
                         "assigned_role": r["role"], "original_assignment": "manual_context_curation",
                         "annotation_product": f["product"], "annotation_label": "",
                         "existing_domain_annotations": "not_measured_in_this_audit",
                         "exact_published_CDS_matches": ";".join(x["feature_id"] for x in matches),
                         "exact_match_labels": ";".join(sorted({x["label"] for x in matches if x["label"]})),
                         "exact_match_products": ";".join(sorted({x["product"] for x in matches if x["product"]})),
                         "function_status": "inferred_role_not_per_locus_experimental_validation"})
    for row in evidence:
        domain = domain_rows[row["protein_id"]]
        row["measured_Pfam_panel_domains"] = domain["domains"]
        row["measured_Pfam_panel_accessions"] = domain["domain_accessions"]
        row["domain_panel_status"] = domain["status"]
    write_table(tables / "diversity_role_evidence.tsv", evidence)

    interim = ROOT / "data/interim/diversity_baseline"
    interim.mkdir(parents=True, exist_ok=True)
    database = interim / "reviewed_bacterial.faa"
    with gzip.open(raw / "reviewed_bacterial_proteins.faa.gz", "rb") as source, database.open("wb") as dest:
        shutil.copyfileobj(source, dest)
    db_count = sum(line.startswith(">") for line in database.open())
    blast = ROOT / ".env/bin"
    with (tables / "diversity_baseline_blast.log").open("w") as log:
        subprocess.run([str(blast / "makeblastdb"), "-in", str(database), "-dbtype", "prot"], check=True, stdout=log, stderr=log)
        fields = "6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp"
        subprocess.run([str(blast / "blastp"), "-query", str(processed / "maize_direct_candidates.faa"),
                        "-db", str(database), "-evalue", "1e-5", "-max_target_seqs", "10", "-num_threads", "4",
                        "-outfmt", fields, "-out", str(tables / "diversity_maize_vs_reviewed.tsv")], check=True, stdout=log, stderr=log)
    titles = {r.id: r.description for r in SeqIO.parse(database, "fasta")}
    hits = defaultdict(list)
    with (tables / "diversity_maize_vs_reviewed.tsv").open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            hits[f[0]].append(f)
    comparison = []
    for candidate in candidates:
        found = sorted(hits[candidate["target_id"]], key=lambda r: (-float(r[7]), r[1]))
        best = found[0] if found else None
        comparison.append({"candidate_id": candidate["target_id"], "historical_status": candidate["curation_status"],
                           "annotation_product": neighbors[candidate["target_id"]]["product"],
                           "reviewed_subject": best[1] if best else "", "reviewed_title": titles[best[1]] if best else "",
                           "reviewed_bitscore": best[7] if best else "", "reviewed_query_coverage": best[8] if best else "",
                           "reviewed_aligned_fraction_of_subject": float(best[3])/float(best[5]) if best else "",
                           "historical_lasso_bitscore": candidate["bitscore"],
                           "historical_AsnB_bitscore": candidate["asnb_bitscore"],
                           "measured_Pfam_panel_domains": domain_rows[candidate["target_id"]]["domains"],
                           "interpretation": "nearest_reviewed_homolog_not_function_assignment"})
    write_table(tables / "diversity_broad_comparison.tsv", comparison)
    summary = {"historical_cyclase_database": dict(Counter(r["source_groups"] for r in rows(processed / "literature_unique_cyclases.tsv"))),
               "historical_reviewed_AsnB": 5, "new_reviewed_bacterial_sequences": db_count,
               "candidates_compared": len(comparison), "roles_audited": len(evidence),
               "distinct_Bacillus_biosamples": len({r["biosample"] for r in novelty}),
               "distinct_Bacillus_defined_DNA_cores": len({r["core_sha256"] for r in novelty}),
               "distinct_Bacillus_C_sequences": len({r["C_sha256"] for r in novelty}),
               "published_records_compared": len(references)}
    (tables / "diversity_baseline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
