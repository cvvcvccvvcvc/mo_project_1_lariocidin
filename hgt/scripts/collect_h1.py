#!/usr/bin/env python3
"""Collect completed H1 shards, retaining source and sequence evidence separately."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

from prepare_h1 import ROOT, TABLES, read_tsv, write_tsv
from audit_h1 import bare


def true(value):
    return str(value) == "True"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-root", type=Path, action="append")
    args = parser.parse_args()
    roots = args.shard_root or [ROOT / "results/H1_shards"]
    quality, mappings, checks, errors = [], [], [], []
    role_rows = None
    for root in roots:
        for directory in sorted(p for p in root.iterdir() if p.is_dir()):
            summary = json.loads((directory / "H1_correspondence_summary.json").read_text())
            if summary["parse_errors"]:
                raise ValueError(f"Failed input in {directory}")
            quality.extend(read_tsv(directory / "H1_genome_quality.tsv"))
            for filename, target in [("H1_locus_genome_correspondence.tsv", mappings),
                                     ("H1_gene_genome_correspondence.tsv", checks)]:
                if (directory / filename).exists():
                    target.extend(read_tsv(directory / filename))
            errors.extend(json.loads((directory / "H1_genome_parse_errors.json").read_text()))
            current_roles = read_tsv(directory / "H1_C_role_evidence.tsv")
            if role_rows is None:
                role_rows = current_roles
            elif current_roles != role_rows:
                raise ValueError("C role inputs differ between runs")
    if len({r["genome_id"] for r in quality}) != len(quality):
        raise ValueError("Repeated genome inputs across shards")
    frozen_taxonomy = {}
    for path in (ROOT / "inputs/H1_2026-10-13/legacy_metadata").glob("*.json"):
        result = json.loads(path.read_text()).get("result", {})
        for uid in result.get("uids", []):
            record = result[uid]
            if record.get("assemblyaccession"):
                frozen_taxonomy[record["assemblyaccession"]] = record
    for record in read_tsv(TABLES / "H1_relative_metadata_pool.tsv"):
        frozen_taxonomy[record["assembly_accession"]] = record
        if record["gbrs_paired_asm"] != "na":
            frozen_taxonomy[record["gbrs_paired_asm"]] = record
    for genome in quality:
        metadata = frozen_taxonomy.get(genome["requested_assembly_accession"], frozen_taxonomy.get(genome["sequence_provider_assembly"], {}))
        genome["species_taxid"] = genome["species_taxid"] or metadata.get("species_taxid", metadata.get("speciestaxid", ""))
        genome["taxid"] = genome["taxid"] or metadata.get("taxid", "")
        genome["taxonomy_basis"] = "frozen_NCBI_assembly_or_native_record; taxid_not_an_independent_host_tree"
    loci = {r["locus_id"]: r for r in read_tsv(TABLES / "H1_publication_loci.tsv")}
    genes = {r["gene_id"]: r for r in read_tsv(TABLES / "H1_publication_genes.tsv")}
    genome_quality = {r["genome_id"]: r for r in quality}
    named_anchors = {r["locus_id"]: r for r in read_tsv(ROOT / "config/H1_named_anchor_review.tsv")}
    roles = {r["locus_id"]: r for r in role_rows}
    lineage_review = {r["locus_id"]: r for r in read_tsv(ROOT / "config/H1_lineage_review.tsv")}
    for locus_id, review in lineage_review.items():
        roles[locus_id]["C_gene_id"] = review["C_candidate_gene_id"]
        roles[locus_id]["role_basis"] = "reviewed_full_length_C_identity_to_noncanonical_reference; lineage_uncertain"
    ipg = read_tsv(ROOT / "inputs/H1_2026-10-13/ncbi/unlinked_C_ipg.tsv")
    ipg_by_protein = defaultdict(list)
    for row in ipg:
        ipg_by_protein[row["Protein"]].append(row)
    by_locus, by_mapping = defaultdict(list), defaultdict(list)
    for row in checks:
        gene = genes[row["gene_id"]]
        touches_edge = int(gene["start_0based"]) == 0 or int(gene["end_0based"]) == int(loci[row["locus_id"]]["record_bp"])
        row["touches_published_export_boundary"] = touches_edge
        row["translation_interpretation"] = "projected_coordinates_or_sequence_change_requires_review" if row["DNA_status"] != "exact_coding_DNA" else \
            "export_boundary_truncation_possible_despite_exact_location_flags" if touches_edge and row["translation_status"] == "translation_mismatch_requires_review" else \
            "annotation_translation_requires_review" if row["translation_status"] == "translation_mismatch_requires_review" else row["translation_status"]
        if row["locus_id"] in lineage_review and row["gene_id"] == lineage_review[row["locus_id"]]["C_candidate_gene_id"]:
            row["role_evidence"] = "C"
        by_mapping[row["mapping_id"]].append(row)
    for row in mappings:
        c_id = roles[row["locus_id"]]["C_gene_id"]
        if row["C_start_0based"] == "":
            c_checks = [g for g in by_mapping[row["mapping_id"]] if g["gene_id"] == c_id]
            if c_checks:
                row["C_start_0based"] = c_checks[0]["projected_start_0based"]
                row["C_end_0based"] = c_checks[0]["projected_end_0based"]
        row["analytical_lineage"] = lineage_review.get(row["locus_id"], {}).get("analytical_lineage", row["lineage"])
        protein_id = genes.get(c_id, {}).get("protein_id", "")
        # WP accessions identify proteins shared by strains, not an original carrier.
        ipg_source = [] if protein_id.startswith("WP_") else ipg_by_protein[protein_id]
        supported = [x for x in ipg_source if row["C_start_0based"] != "" and bare(x["Nucleotide Accession"]) == bare(row["contig"])
                     and int(x["Start"]) - 1 == int(row["C_start_0based"])
                     and int(x["Stop"]) == int(row["C_end_0based"])]
        row["specific_protein_IPG_coordinate_match"] = bool(supported)
        named = named_anchors.get(row["locus_id"], {})
        strain = genome_quality[row["genome_id"]]["strain"].removeprefix("strain=")
        normalize = lambda value: re.sub(r"[^A-Z0-9]", "", value.upper())
        named_support = bool(named) and row["sequence_evidence"] == "whole_export_exact" and normalize(strain) == normalize(named["published_strain"])
        row["reviewed_named_anchor_whole_DNA_and_strain_match"] = named_support
        primary = any(true(row[k]) for k in ["source_nucleotide_base_match", "nucleotide_assembly_link",
                                             "published_NPDC_strain_match"]) or bool(supported) or named_support
        row["carrier_assignment"] = "source_provenance_supported" if primary else "sequence_occurrence_only; published_carrier_not_assigned"
        by_locus[row["locus_id"]].append(row)
    coordinate_key = lambda row: (row["genome_id"], row["contig"], row["C_start_0based"], row["C_end_0based"])
    native_c_evidence = {coordinate_key(row): row for row in mappings if row["C_distance_to_contig_edge_bp"] != ""}
    for row in mappings:
        if row["C_distance_to_contig_edge_bp"] == "" and coordinate_key(row) in native_c_evidence:
            evidence = native_c_evidence[coordinate_key(row)]
            row["C_distance_to_contig_edge_bp"] = evidence["C_distance_to_contig_edge_bp"]
            row["identical_C_protein_annotation_copies_in_genome"] = evidence["identical_C_protein_annotation_copies_in_genome"]
    catalogue = []
    for locus_id, locus in loci.items():
        matches = by_locus[locus_id]
        primary = [r for r in matches if r["carrier_assignment"] == "source_provenance_supported"]
        full = [r for r in primary if r["sequence_evidence"] == "whole_export_exact"]
        native_c = [g for r in primary for g in by_mapping[r["mapping_id"]] if g["role_evidence"] == "C"]
        c_ok = any(g["DNA_status"] == "exact_coding_DNA" and g["translation_status"] in
                   ["exact_translation", "exact_with_bacterial_initiator_M"] for g in native_c)
        catalogue.append({**locus, "lineage_label_basis": "initial_technical_assignment" if locus["source"] == "Nature_2025" else "Merrild_archive_directory",
                          "analytical_lineage": lineage_review.get(locus_id, {}).get("analytical_lineage", locus["lineage"]),
                          "C_gene_id": roles[locus_id]["C_gene_id"], "C_role_basis": roles[locus_id]["role_basis"],
                          "primary_genome_inputs": ";".join(sorted({r["genome_id"] for r in primary})),
                          "other_sequence_occurrences": ";".join(sorted({r["genome_id"] for r in matches if r not in primary})),
                          "mapping_status": "whole_export_and_carrier_supported" if full else
                          "C_and_carrier_supported; context_requires_review" if primary else
                          "sequence_occurrence_without_primary_carrier" if matches else "not_resolved_in_frozen_genomes",
                          "C_DNA_and_translation_confirmed": c_ok,
                          "C_partial_annotation": genes.get(roles[locus_id]["C_gene_id"], {}).get("partial_location", "unassigned"),
                          "H2_admission": "lineage_and_function_review_required" if locus_id in lineage_review else
                          "sequence_and_carrier_ready; component_curation_and_marker_QC_pending" if c_ok else
                          "manual_sequence_or_carrier_review_required",
                          "primary_mapping_count": len(primary), "sequence_occurrence_count": len(matches),
                          "HGT_status": "not_assessed"})
    # Physical intervals use overlap of C, not equality of variable export boundaries.
    buckets = defaultdict(list)
    for row in mappings:
        if row["C_start_0based"] != "":
            buckets[(row["genome_id"], row["contig"])].append(row)
    physical = []
    for key, rows in sorted(buckets.items()):
        groups = []
        for row in sorted(rows, key=lambda r: (int(r["C_start_0based"]), int(r["C_end_0based"]))):
            left, right = int(row["C_start_0based"]), int(row["C_end_0based"])
            if groups and left < groups[-1][1]:
                groups[-1][1] = max(groups[-1][1], right)
                groups[-1][2].append(row)
            else:
                groups.append([left, right, [row]])
        for left, right, members in groups:
            lineages = sorted({r["analytical_lineage"] for r in members})
            physical.append({"physical_locus_id": f"{key[0]}|{key[1]}|C:{left}-{right}",
                             "genome_id": key[0], "contig": key[1], "C_interval_start_0based": left,
                             "C_interval_end_0based": right, "lineages": ";".join(lineages),
                             "lineage_conflict": len(set(lineages) - {"uncertain_Nature_system"}) > 1,
                             "publication_locus_ids": ";".join(sorted({r["locus_id"] for r in members})),
                             "publication_export_count": len({r["locus_id"] for r in members}),
                             "source_supported_export_count": len({r["locus_id"] for r in members if r["carrier_assignment"] == "source_provenance_supported"}),
                             "limitation": "exact_C_matches_only; divergent_copies_not_excluded"})
    duplicate_rows = []
    for basis in ["sequence_content_sha256", "biosample", "sequence_provider_assembly"]:
        grouped = defaultdict(list)
        for row in quality:
            if row[basis]:
                grouped[row[basis]].append(row)
        for value, members in sorted(grouped.items()):
            if len(members) > 1:
                duplicate_rows.append({"basis": basis, "value": value,
                                       "genome_inputs": ";".join(sorted(r["genome_id"] for r in members)),
                                       "biosamples": ";".join(sorted({r["biosample"] for r in members})),
                                       "interpretation": "identical_genomic_sequence; preserve_samples; hash_alone_does_not_establish_isolate_identity" if basis == "sequence_content_sha256" else
                                       "same_provider_assembly; technical_redundancy" if basis == "sequence_provider_assembly" else
                                       "same_sample; assembly_versions_are_not_independent_carriers"})
    component_rows = []
    for row in mappings:
        if row["carrier_assignment"] != "source_provenance_supported":
            continue
        for gene in by_mapping[row["mapping_id"]]:
            if not gene["role_evidence"]:
                continue
            component_rows.append({**gene, "lineage": row["lineage"], "sequence_evidence": row["sequence_evidence"],
                                   "functional_curation": "published_role_or_C_homology; domains_and_orthology_pending"})
    candidate_rows = []
    b_hits = TABLES / "H1_B_homology.tsv"
    if b_hits.exists():
        best = {}
        for line in b_hits.open():
            fields = line.rstrip().split("\t")
            query, reference = fields[:2]
            qcov = (int(fields[5]) - int(fields[4]) + 1) / int(fields[10])
            scov = (int(fields[7]) - int(fields[6]) + 1) / int(fields[11])
            if min(qcov, scov) < 0.7:
                continue
            key = query, genes[reference]["author_role"]
            hit = {"gene_id": query, "locus_id": genes[query]["locus_id"], "component_candidate": key[1],
                   "reference_gene_id": reference, "identity_percent": fields[2], "bitscore": float(fields[9]),
                   "E_value": fields[8], "query_coverage": qcov, "reference_coverage": scov,
                   "functional_curation": "homology_candidate; orthology_domains_and_context_require_review"}
            if key not in best or hit["bitscore"] > best[key]["bitscore"]:
                best[key] = hit
        for hit in sorted(best.values(), key=lambda r: (r["gene_id"], r["component_candidate"])):
            native = [g for m in by_locus[hit["locus_id"]] if m["carrier_assignment"] == "source_provenance_supported"
                      for g in by_mapping[m["mapping_id"]] if g["gene_id"] == hit["gene_id"]]
            hit["source_DNA_and_translation_confirmed"] = any(g["DNA_status"] == "exact_coding_DNA" and
                g["translation_status"] in ["exact_translation", "exact_with_bacterial_initiator_M"] for g in native)
            candidate_rows.append(hit)
        if candidate_rows:
            write_tsv(TABLES / "H1_component_candidates.tsv", candidate_rows)
    for name, rows in [("genome_quality", quality), ("locus_genome_correspondence", mappings),
                       ("gene_genome_correspondence", checks), ("C_role_evidence", role_rows),
                       ("carrier_catalogue", catalogue), ("physical_loci", physical),
                       ("genome_redundancy", duplicate_rows), ("component_correspondence", component_rows)]:
        if rows:
            write_tsv(TABLES / f"H1_{name}.tsv", rows)
    for basis, rows, id_field, filename in [("dna_sha256", catalogue, "locus_id", "export_redundancy"),
                                           ("protein_sha256", [dict(genes[r["C_gene_id"]], locus_id=r["locus_id"],
                                            analytical_lineage=r["analytical_lineage"]) for r in catalogue if r["C_gene_id"]],
                                            "locus_id", "C_redundancy")]:
        groups = defaultdict(list)
        for row in rows:
            groups[row[basis]].append(row)
        repeated = [{"sequence_sha256": key, "publication_export_count": len(members),
                     "publication_locus_ids": ";".join(sorted(r[id_field] for r in members)),
                     "analytical_lineages": ";".join(sorted({r["analytical_lineage"] for r in members})),
                     "interpretation": "sequence_identity; carrier_occurrences_retained; not_an_HGT_event_count"}
                    for key, members in sorted(groups.items()) if key and len(members) > 1]
        if repeated:
            write_tsv(TABLES / f"H1_{filename}.tsv", repeated)
    quality_flags = []
    for row in mappings:
        if row["carrier_assignment"] != "source_provenance_supported":
            continue
        flags = []
        if row["C_distance_to_contig_edge_bp"] != "" and int(row["C_distance_to_contig_edge_bp"]) <= 1000:
            flags.append("C_within_1kb_of_native_contig_edge; neighbourhood_or_C_completeness_requires_review")
        if int(row["identical_C_protein_annotation_copies_in_genome"]) > 1:
            flags.append("multiple_identical_C_annotations; retain_separate_coordinates")
        if row["sequence_evidence"] != "whole_export_exact":
            flags.append("C_match_only; projected_flanks_not_assumed_identical")
        if flags:
            quality_flags.append({"mapping_id": row["mapping_id"], "locus_id": row["locus_id"], "genome_id": row["genome_id"],
                                  "flags": ";".join(flags), "C_contig_edge_bp": row["C_distance_to_contig_edge_bp"],
                                  "C_exact_annotation_copies": row["identical_C_protein_annotation_copies_in_genome"]})
    if quality_flags:
        write_tsv(TABLES / "H1_mapping_review_flags.tsv", quality_flags)
    primary_mappings = [r for r in mappings if r["carrier_assignment"] == "source_provenance_supported"]
    source_genomes = {r["genome_id"] for r in primary_mappings}
    source_hashes = {r["sequence_content_sha256"] for r in quality if r["genome_id"] in source_genomes}
    by_accession = {r["requested_assembly_accession"]: r for r in quality}
    relative_panel = []
    selection_path = TABLES / "H1_relative_selection.tsv"
    for row in read_tsv(selection_path) if selection_path.exists() else []:
        genome = by_accession.get(row["assembly_accession"], {})
        relative_panel.append({**row, "downloaded_genome_input": genome.get("genome_id", ""),
                               "realised_role": "published_carrier_recovered; retain_as_target" if genome.get("genome_id") in source_genomes else "taxonomic_relative_control",
                               "native_genome_N50_bp": genome.get("N50_bp", ""), "native_contig_records": genome.get("records", ""),
                               "sequence_identical_to_a_source_input": genome.get("sequence_content_sha256") in source_hashes,
                               "marker_QC": genome.get("marker_tree_input_status", "not_available"),
                               "BGC_presence_absence": "not_assessed_by_sensitive_search"})
    if relative_panel:
        write_tsv(TABLES / "H1_relative_panel.tsv", relative_panel)
    primary_mapping_ids = {r["mapping_id"] for r in primary_mappings}
    primary_gene_checks = [r for r in checks if r["mapping_id"] in primary_mapping_ids]
    summary = {"genome_inputs": len(quality), "publication_exports": len(loci), "sequence_mappings": len(mappings),
               "gene_checks": len(checks), "physical_locus_occurrences": len(physical), "parse_errors": len(errors),
               "C_role_unassigned_exports": sum(not r["C_gene_id"] for r in catalogue),
               "lineages": {lineage: {"exports": sum(r["analytical_lineage"] == lineage for r in catalogue),
                             "statuses": dict(Counter(r["mapping_status"] for r in catalogue if r["analytical_lineage"] == lineage)),
                             "C_sequence_and_carrier_confirmed": sum(r["analytical_lineage"] == lineage and r["C_DNA_and_translation_confirmed"] for r in catalogue)}
                            for lineage in ["canonical", "noncanonical", "uncertain_Nature_system"]},
               "DNA_check_statuses": dict(Counter(r["DNA_status"] for r in checks)),
               "translation_check_statuses": dict(Counter(r["translation_status"] for r in checks))}
    summary.update({"source_supported_exports": sum(bool(r["primary_genome_inputs"]) for r in catalogue),
                    "source_genome_inputs": len(source_genomes),
                    "source_physical_locus_occurrences": sum(int(r["source_supported_export_count"]) > 0 for r in physical),
                    "source_mapping_count": len(primary_mappings),
                    "source_mapping_evidence": dict(Counter(r["sequence_evidence"] for r in primary_mappings)),
                    "source_gene_checks": len(primary_gene_checks),
                    "source_gene_translation_interpretation": dict(Counter(r["translation_interpretation"] for r in primary_gene_checks)),
                    "source_assembly_levels": dict(Counter(r["assembly_level"] for r in quality if r["genome_id"] in source_genomes)),
                    "source_min_N50": min(int(r["N50_bp"]) for r in quality if r["genome_id"] in source_genomes),
                    "source_max_records": max(int(r["records"]) for r in quality if r["genome_id"] in source_genomes),
                    "source_no_native_translations": [r["genome_id"] for r in quality if r["genome_id"] in source_genomes and not int(r["translated_cds_count"])],
                    "relative_selection_count": len(relative_panel),
                    "relative_realised_roles": dict(Counter(r["realised_role"] for r in relative_panel)),
                    "relative_sequence_identical_source_count": sum(r["sequence_identical_to_a_source_input"] for r in relative_panel),
                    "B_candidates": {role: {"genes": sum(r["component_candidate"] == role for r in candidate_rows),
                                            "exports": len({r["locus_id"] for r in candidate_rows if r["component_candidate"] == role}),
                                            "source_DNA_translation_confirmed": sum(r["component_candidate"] == role and r["source_DNA_and_translation_confirmed"] for r in candidate_rows)} for role in ["B1", "B2"]},
                    "unique_export_sequences": len({r["dna_sha256"] for r in catalogue}),
                    "unique_C_protein_sequences": len({genes[r["C_gene_id"]]["protein_sha256"] for r in catalogue if r["C_gene_id"]}),
                    "source_edge_flags": sum(r["C_distance_to_contig_edge_bp"] != "" and int(r["C_distance_to_contig_edge_bp"]) <= 1000 for r in primary_mappings)})
    (TABLES / "H1_collected_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (TABLES / "H1_review_metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
