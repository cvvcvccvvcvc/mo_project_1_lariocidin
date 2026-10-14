#!/usr/bin/env python3
"""Curate physical gene copies using frozen DNA, reference homology and context."""
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_h2 import ROOT, SNAP, H1, OUT, read, write


def digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def domains(path):
    result = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            result[f[0]].append({"accession": f[4].split(".")[0], "model": f[3], "start_aa": int(f[17]) - 1,
                                 "end_aa": int(f[18]), "E_value": f[12], "bitscore": float(f[13])})
    return result


def main():
    hosts = read(OUT / "host_inputs.tsv")
    aliases = {r["genome_id"]: r["host_id"] for r in read(OUT / "host_aliases.tsv")}
    source_inputs = {r["genome_id"]: r for r in read(H1 / "H1_all_assembly_inputs.tsv")}
    catalogue = {r["locus_id"]: r for r in read(H1 / "H1_carrier_catalogue.tsv")}
    genes = {r["gene_id"]: r for r in read(H1 / "H1_publication_genes.tsv")}
    pub_proteins = {r.id: str(r.seq) for r in SeqIO.parse(H1 / "proteins.faa", "fasta")}
    pub_coding = {r.id: str(r.seq).upper() for r in SeqIO.parse(H1 / "coding_sequences.fna", "fasta")}
    refs = {r["gene_id"]: r["family"] for r in read(OUT / "role_references.tsv")}
    publication_domains = domains(OUT / "publication_domains.domtbl")
    mappings = defaultdict(list)
    for row in read(H1 / "H1_locus_genome_correspondence.tsv"):
        if row["carrier_assignment"] == "source_provenance_supported" or row["sequence_evidence"] == "whole_export_exact":
            mappings[aliases[row["genome_id"]]].append(row)
    native_root = ROOT
    native_root = type(ROOT)(os.environ.get("HGT_GENOME_ROOT", str(native_root)))
    variants, components, precursor_rows, errors, diagnostic_background = [], [], [], [], []
    fasta = defaultdict(dict)
    qc = {r["host_id"]: r for r in read(OUT / "marker_genome_QC.tsv")}
    for host in hosts:
        host_id = host["host_id"]
        directory = OUT / "native" / host_id
        native_metadata = read(directory / "proteins.tsv")
        native_sequences = {r.id: str(r.seq) for r in SeqIO.parse(directory / "proteins.faa", "fasta")}
        native_domains = domains(directory / "domains.domtbl")
        representative_dna = {r.id: str(r.seq).upper() for r in SeqIO.parse(directory / "genome.fna", "fasta")}
        contigs = read(directory / "contigs.tsv")
        by_hash = defaultdict(list)
        for row in contigs:
            by_hash[row["canonical_contig_sha256"]].append(row)
        source_dna = {host["genome_id"]: representative_dna}
        c_candidates = []
        for mapping in mappings[host_id]:
            locus = catalogue[mapping["locus_id"]]
            gene = genes[locus["C_gene_id"]]
            origin = mapping["genome_id"]
            if origin not in source_dna:
                path = native_root / source_inputs[origin]["path"]
                opening = gzip.open if path.suffix == ".gz" else open
                with opening(path, "rt") as handle, warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    source_dna[origin] = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, "genbank")}
            original = source_dna[origin][mapping["contig"]]
            reverse = str(Seq(original).reverse_complement())
            canonical = min(original, reverse)
            representative_matches = by_hash[digest(canonical)]
            if len(representative_matches) != 1:
                errors.append({"locus_id": locus["locus_id"], "host_id": host_id, "reason": "ambiguous_or_missing_contig_alias; copies_not_arbitrarily_assigned"})
                continue
            representative = representative_matches[0]
            flip = (original != canonical) != (representative["native_to_canonical_reverse"] == "True")
            sequence = representative_dna[representative["contig"]]
            left = int(mapping["export_start_0based"])
            orientation = int(mapping["orientation"])
            parts, extracted = [], []
            for a, b, strand in json.loads(gene["parts"]):
                a, b = (left + a, left + b) if orientation == 1 else (left + int(locus["record_bp"]) - b, left + int(locus["record_bp"]) - a)
                strand *= orientation
                if flip:
                    a, b, strand = len(sequence) - b, len(sequence) - a, -strand
                part = sequence[a:b]
                extracted.append(part if strand == 1 else str(Seq(part).reverse_complement()))
                parts.append((a, b, strand))
            if "".join(extracted) != pub_coding[gene["gene_id"]]:
                errors.append({"locus_id": locus["locus_id"], "host_id": host_id, "reason": "native_C_coding_DNA_does_not_match_frozen_publication; review_required"})
                continue
            start, end = min(x[0] for x in parts), max(x[1] for x in parts)
            strand = parts[0][2]
            exact = [r for r in native_metadata if r["contig"] == representative["contig"] and int(r["start_0based"]) == start and int(r["end_0based"]) == end and int(r["strand"]) == strand and native_sequences[r["protein_id"]] == pub_proteins[gene["gene_id"]]]
            c_candidates.append({"host_id": host_id, "contig": representative["contig"], "start": start, "end": end, "strand": strand,
                                 "locus_id": locus["locus_id"], "lineage": locus["analytical_lineage"], "gene_id": gene["gene_id"],
                                 "protein": pub_proteins[gene["gene_id"]], "native_protein_id": exact[0]["protein_id"] if len(exact) == 1 else "",
                                 "native_model_exact": len(exact) == 1, "source_assignment": mapping["carrier_assignment"],
                                 "C_contig_edge_bp": min(start, len(sequence) - end), "annotation_variant": len(exact) != 1})
        grouped = defaultdict(list)
        for candidate in c_candidates:
            grouped[(candidate["lineage"], candidate["contig"], candidate["strand"])].append(candidate)
        physical_c = []
        for (lineage, contig, strand), rows in grouped.items():
            blocks = []
            for row in sorted(rows, key=lambda r: (r["start"], r["end"])):
                if blocks and row["start"] < blocks[-1][1]:
                    blocks[-1][1] = max(blocks[-1][1], row["end"])
                    blocks[-1][2].append(row)
                else:
                    blocks.append([row["start"], row["end"], [row]])
            for start, end, members in blocks:
                members.sort(key=lambda r: (not r["native_model_exact"], -len(r["protein"]), r["gene_id"]))
                chosen = members[0]
                key = f"{host_id}|{contig}|{start}-{end}|{strand}"
                tip = host_id + "_C" + digest(key)[:10]
                family = lineage + "_C"
                different_models = len({r["protein"] for r in members}) > 1
                for member in members:
                    variants.append({"tip_id": tip, "family": family, "host_id": host_id, "publication_locus_id": member["locus_id"],
                                     "publication_gene_id": member["gene_id"], "protein_aa": len(member["protein"]),
                                     "native_model_exact": member["native_model_exact"], "selected_model": member is chosen,
                                     "annotation_models_disagree": different_models,
                                     "decision": "prefer_exact_native_model; otherwise_longest_DNA_verified_publication_model; alternatives_retained"})
                domain_hits = publication_domains.get(chosen["gene_id"], [])
                row = {"tip_id": tip, "family": family, "host_id": host_id, "contig": contig,
                       "start_0based": chosen["start"], "end_0based": chosen["end"], "strand": strand,
                       "protein_aa": len(chosen["protein"]), "native_protein_id": chosen["native_protein_id"],
                       "publication_locus_ids": ";".join(sorted({r["locus_id"] for r in members})),
                       "publication_gene_ids": ";".join(sorted({r["gene_id"] for r in members})),
                       "role_basis": "DNA_verified_publication_C_and_context; broad_Asn_GATase_domains_not_unique_function",
                       "domain_accessions": ";".join(sorted({x["accession"] for x in domain_hits})),
                       "native_model_exact": chosen["native_model_exact"], "annotation_models_disagree": different_models,
                       "C_contig_edge_bp": chosen["C_contig_edge_bp"], "host_QC_pass": qc[host_id]["host_tree_admission"] == "True",
                       "orthology_status": "homologous_C_family; one_to_one_orthology_not_assumed_for_multiple_copies"}
                components.append(row)
                fasta[family][tip] = chosen["protein"]
                physical_c.append({**chosen, "family": family, "tip_id": tip})
        # Full-length reference homology plus local C context; retain every copy.
        hits = defaultdict(list)
        with (directory / "role_homology.tsv").open() as handle:
            for line in handle:
                f = line.strip().split("\t")
                qcov = (int(f[5]) - int(f[4]) + 1) / int(f[10])
                scov = (int(f[7]) - int(f[6]) + 1) / int(f[11])
                if min(qcov, scov) >= 0.7:
                    family = refs[f[1]]
                    if family.startswith("canonical_ABC_"):
                        family = "canonical_ABC_family"
                    hits[(f[0], family)].append((float(f[9]), f[1], float(f[2]), qcov, scov))
        proteins_by_id = {r["protein_id"]: r for r in native_metadata}
        for (protein_id, family), matches in hits.items():
            if family.endswith("_C"):
                continue
            gene = proteins_by_id[protein_id]
            lineage = "noncanonical" if family.startswith("noncanonical_") else "canonical"
            nearby = [r for r in physical_c if r["lineage"] == lineage and r["contig"] == gene["contig"] and
                      max(0, int(gene["start_0based"]) - r["end"], r["start"] - int(gene["end_0based"])) <= 5000]
            if not nearby or gene["partial_location"] == "True":
                continue
            best = max(matches)
            tip = host_id + "_G" + digest(f"{gene['contig']}|{gene['start_0based']}|{gene['end_0based']}|{gene['strand']}")[:10]
            domain_hits = native_domains.get(protein_id, [])
            components.append({"tip_id": tip, "family": family, "host_id": host_id, "contig": gene["contig"],
                               "start_0based": int(gene["start_0based"]), "end_0based": int(gene["end_0based"]), "strand": int(gene["strand"]),
                               "protein_aa": gene["protein_aa"], "native_protein_id": protein_id,
                               "publication_locus_ids": ";".join(sorted({r["locus_id"] for r in nearby})), "publication_gene_ids": "",
                               "role_basis": f"full_reference_homology_to_{best[1]}; identity={best[2]}; qcov={best[3]:.3f}; scov={best[4]:.3f}; within5kb_of_verified_C",
                               "domain_accessions": ";".join(sorted({x["accession"] for x in domain_hits})),
                               "native_model_exact": True, "annotation_models_disagree": False, "C_contig_edge_bp": "not_C",
                               "host_QC_pass": qc[host_id]["host_tree_admission"] == "True",
                               "orthology_status": "reference_homolog_in_context; family_tree_and_copy_number_required; function_not_assumed"})
            fasta[family][tip] = native_sequences[protein_id]
        # Independent metabolic-family homologs outside mapped loci for C diagnosis.
        for gene in native_metadata:
            if gene["overlaps_mapped_export_window"] == "False" and any(x["accession"] == "PF00733" for x in native_domains.get(gene["protein_id"], [])):
                diagnostic_background.append({"host_id": host_id, "native_protein_id": gene["protein_id"], "contig": gene["contig"],
                                               "start_0based": gene["start_0based"], "end_0based": gene["end_0based"], "product": gene["product"],
                                               "basis": "Asn_synthase_domain_outside_mapped_loci; metabolic_function_not_proven_by_annotation"})
    # Precursors remain sequences/alleles; their short lengths do not force a tree.
    by_locus = defaultdict(list)
    for gene in genes.values():
        if gene["author_role"] in ["A", "A1"] or gene["gene_id"] == "Nature_2025_91bc539987e9.f1":
            by_locus[gene["locus_id"]].append(gene)
    admitted_precursor_checks = {r["gene_id"]: r for r in read(H1 / "H1_gene_genome_correspondence.tsv")
                                if r["DNA_status"] == "exact_coding_DNA" and r["translation_status"] in ["exact_translation", "exact_with_bacterial_initiator_M"]}
    for locus_id, locus in catalogue.items():
        for gene in by_locus[locus_id]:
            precursor_rows.append({"locus_id": locus_id, "gene_id": gene["gene_id"], "analytical_lineage": locus["analytical_lineage"],
                                   "protein_aa": gene["protein_aa"], "amino_acid_sequence": pub_proteins.get(gene["gene_id"], ""),
                                   "coding_DNA": pub_coding[gene["gene_id"]], "source_genome_inputs": locus["primary_genome_inputs"],
                                   "native_DNA_translation_confirmed": gene["gene_id"] in admitted_precursor_checks,
                                   "published_annotation": gene["product"], "tree_status": "sequence_allele_only; standalone_tree_not_forced"})
    if precursor_rows:
        write(OUT / "precursor_alleles.tsv", precursor_rows)
    gene_directory = OUT / "components"
    gene_directory.mkdir(exist_ok=True)
    for family, seqs in fasta.items():
        SeqIO.write((SeqRecord(Seq(sequence), id=tip, description="") for tip, sequence in sorted(seqs.items())), gene_directory / (family + ".faa"), "fasta")
    write(OUT / "component_leaf_map.tsv", components)
    write(OUT / "C_annotation_variants.tsv", variants)
    if diagnostic_background:
        write(OUT / "C_metabolic_background_candidates.tsv", diagnostic_background)
    (OUT / "component_mapping_errors.json").write_text(json.dumps(errors, indent=2) + "\n")
    summary = {"families": {family: {"physical_gene_copies": len(seqs), "unique_sequences": len(set(seqs.values()))} for family, seqs in fasta.items()},
               "mapping_errors_requiring_review": len(errors), "C_variants": len(variants), "metabolic_background_candidates": len(diagnostic_background)}
    (OUT / "component_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
