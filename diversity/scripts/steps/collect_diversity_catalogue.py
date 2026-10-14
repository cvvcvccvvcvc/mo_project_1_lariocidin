#!/usr/bin/env python3
"""Combine completed species/context shards and compare defined sequences.

This consolidates evidence; it does not infer biochemical function from a
closed comparator list or claim global novelty from an absent exact match.
"""
import argparse
import hashlib
import json
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from profile_search import ROOT


def read(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def sequence_index(paths):
    index = defaultdict(list)
    for path in paths:
        for record in SeqIO.parse(path, "fasta"):
            index[hashlib.sha256(str(record.seq).encode()).hexdigest()].append(record.id)
    return index


def wgs_project(accession):
    master = re.fullmatch(r"(?:NZ_)?([A-Z]{4,6})0+\.(\d+)", accession)
    if master:
        return master[1] + f"{int(master[2]):02d}"
    match = re.match(r"(?:NZ_)?([A-Z]{4,6}\d{2})\d*", accession)
    return match.group(1) if match else ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--prefix", default="diversity_species")
    parser.add_argument("--output-prefix", help="Read separately recurated shard tables and preserve the original stage outputs")
    parser.add_argument("--selection", type=str)
    parser.add_argument("--refseq-metadata", default="data/processed/discovery/diversity_refseq_metadata.tsv")
    parser.add_argument("--mgnify-metadata", default="inputs/selections/diversity_soil_species.tsv")
    args = parser.parse_args()
    output_prefix = args.output_prefix or args.prefix
    tables, processed = ROOT / "results/tables", ROOT / "data/processed/discovery"
    metadata = read(ROOT / args.refseq_metadata).set_index("assembly")
    soil = read(ROOT / args.mgnify_metadata).set_index("Genome")
    published = read(processed / "diversity_published_records.tsv").set_index("record_id")
    by_contig, by_project, by_strain = defaultdict(list), defaultdict(list), defaultdict(list)
    for identifier, row in published.to_dict("index").items():
        contig = row["accession"].removeprefix("NZ_").split(".")[0]
        project = wgs_project(row["accession"])
        if contig:
            by_contig[contig].append(identifier)
        if project:
            by_project[project].append(identifier)
        for strain in row["strain"].split(";"):
            if strain:
                by_strain[strain].append(identifier)
    published_dna = {r.id: str(r.seq).upper() for r in SeqIO.parse(processed / "diversity_published_DNA.fna", "fasta")}
    published_c = sequence_index([processed / "diversity_published_C.faa"])
    published_c_sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "diversity_published_C.faa", "fasta")}
    published_a = sequence_index([processed.parent / "literature_unique_precursors.faa", processed.parent / "all_cds.faa"])
    reviewed_titles = {r.id: r.description for r in SeqIO.parse(ROOT / "data/interim/diversity_baseline/reviewed_bacterial.faa", "fasta")}
    catalogues, proteomes, hit_tables, role_tables, validation_tables = [], [], [], [], []
    fasta_records, precursor_evidence, core_records = [], [], []
    for shard in range(args.shards):
        name = f"{args.prefix}_{shard}"
        result_name = f"{output_prefix}_{shard}"
        context = processed / "contexts" / f"{name}_contexts"
        proteomes.append(read(tables / f"{name}_proteomes.tsv"))
        hits = read(tables / f"{name}_profile_hits.tsv")
        hit_tables.append(hits)
        validation = read(context / "validation.tsv")
        validation_tables.append(validation)
        if set(hits.target_id) != set(validation.target_id):
            raise ValueError(f"Context validation does not account for every nomination in {name}")
        data = read(tables / f"{result_name}_candidate_catalogue.tsv")
        unresolved = validation[validation.annotated_exact_matches == 0].merge(
            hits.sort_values("sequence_score", ascending=False).drop_duplicates("target_id"),
            on=["target_id", "genome_id"], validate="one_to_one")
        if not unresolved.empty:
            unresolved = unresolved.merge(proteomes[-1][["genome_id", "source", "assembly", "biosample"]],
                                           on="genome_id", validate="many_to_one")
            unresolved["status"] = "context_unresolved"
            unresolved["family"] = "unresolved"
            unresolved["interpretation_reason"] = "Nominated protein not exactly reconciled with downloaded assembly annotation"
            unresolved["locus_id"] = unresolved.target_id
            data = pd.concat([data, unresolved], ignore_index=True).fillna("")
        catalogues.append(data)
        role_tables.append(read(tables / f"{result_name}_role_evidence.tsv"))
        precursor_evidence.append(read(tables / f"{result_name}_precursor_evidence.tsv" if args.output_prefix else context / "precursor_evidence.tsv"))
        sequences = {r.id: r for path in (context / "neighbors.faa", context / "precursor_queries.faa") for r in SeqIO.parse(path, "fasta")}
        windows = SeqIO.to_dict(SeqIO.parse(context / "windows.fna", "fasta"))
        for row in data[data.status == "family_bgc_candidate"].to_dict("records"):
            fasta_records.extend(sequences[x] for x in row["chosen_core_gene_ids"].split(";"))
            window = windows[row["locus_id"]]
            core = window[int(row["core_start"])-int(row["window_start"]):int(row["core_end"])-int(row["window_start"])]
            if int(row["C_strand"]) == -1:
                core = core.reverse_complement(id=True, description=True)
            core_records.append(core)
    proteomes = pd.concat(proteomes, ignore_index=True)
    hits = pd.concat(hit_tables, ignore_index=True)
    data = pd.concat(catalogues, ignore_index=True).fillna("")
    expected = (set(read(ROOT / args.selection).genome_id) if args.selection else
                set(read(ROOT / "inputs/selections/diversity_refseq_species.tsv").assembly) | set(soil.index))
    if len(proteomes) != len(expected) or set(proteomes.genome_id) != expected:
        raise ValueError("Completed proteomes differ from the fixed source selection")
    core_by_locus = {r.id: str(r.seq).upper() for r in core_records}
    core_gene_sequences = {r.id: str(r.seq) for r in fasta_records}
    additions = []
    for row in data.to_dict("records"):
        genome = row["genome_id"]
        if genome in metadata.index:
            source = metadata.loc[genome].to_dict()
            source["taxonomy_source"] = "NCBI_RefSeq"
            source["source_genome_type"] = "primary_metadata_in_sample_type_and_metagenome_source"
        else:
            r = soil.loc[genome]
            labels = dict(x.split("__", 1) for x in r.Lineage.split(";"))
            organism = labels["s"] or (labels["g"] + " sp. (species unassigned)" if labels["g"] else "Unclassified " + (labels["f"] or labels["p"]))
            source = {"organism": organism, "strain": "",
                      "biosample": r.Sample_accession, "wgs_project": wgs_project(r.Genome_accession),
                      "source_genome_type": r.Genome_type, "taxonomy_source": "MGnify_GTDB_lineage",
                      "lineage": r.Lineage, "source_completeness": r.Completeness,
                      "source_contamination": r.Contamination, "N_contigs": r.N_contigs,
                      "N50": r.N50, "country": r.Country, "study": r.Study_accession}
        c_exact = published_c.get(row.get("C_sha256", ""), [])
        a_exact = [published_a.get(x, []) for x in str(row.get("precursor_sha256", "")).split(";") if x]
        core = core_by_locus.get(row["locus_id"], "")
        reverse = str(Seq(core).reverse_complement()) if core else ""
        core_matches = [identifier for identifier, dna in published_dna.items() if core and (core in dna or reverse in dna)]
        query_c = core_gene_sequences.get(row.get("C_query_id", ""), "")
        reference_c = published_c_sequences.get(row.get("C_nearest_reference", ""), "")
        # A shifted CDS start changes the FASTA string even on identical DNA.
        # Compare the shared C-terminal segment, excluding the initiating residue.
        overlap = max(0, min(len(query_c), len(reference_c))-1)
        same_tail = bool(overlap and query_c[-overlap:] == reference_c[-overlap:])
        c_interpretation = "not_evaluated_for_uncertain_context"
        if core:
            if c_exact:
                c_interpretation = "exact_published_C"
            elif core_matches and same_tail and len(query_c) != len(reference_c):
                c_interpretation = "start_annotation_difference_on_published_DNA"
            elif not core_matches and not same_tail:
                c_interpretation = "different_C_amino_acids_and_unrepresented_DNA_span"
            else:
                c_interpretation = "nonexact_C_representation_requires_review"
        contig = str(row.get("contig", "")).removeprefix("NZ_").split(".")[0]
        accession_matches = by_contig.get(contig, [])
        project_matches = by_project.get(source.get("wgs_project", ""), [])
        strain_matches = by_strain.get(source.get("strain", ""), [])
        carrier_comparison = "not_matched_in_exported_carrier_metadata"
        if accession_matches or project_matches:
            carrier_comparison = "represented_contig_or_WGS"
        elif strain_matches:
            carrier_comparison = "represented_strain_label_locus_relationship_unresolved"
        elif core_matches:
            carrier_comparison = "exact_exported_locus_span_carrier_metadata_unresolved"
        additions.append({**source, "published_exact_C_ids": ";".join(c_exact),
                          "C_reviewed_title": reviewed_titles.get(row.get("C_reviewed_reference", ""), ""),
                          "published_exact_A_ids": "|".join(";".join(x) for x in a_exact),
                          "published_exact_core_records": ";".join(core_matches),
                          "published_contig_records": ";".join(accession_matches),
                          "published_WGS_records": ";".join(project_matches),
                          "published_strain_records": ";".join(strain_matches),
                          "C_sequence_comparison": "exact_published" if c_exact else "not_exact_in_frozen_C_set" if row.get("C_sha256") else "unresolved",
                          "C_difference_interpretation": c_interpretation,
                          "C_shared_terminal_segment_aa": overlap,
                          "C_shared_terminal_segment_identical": same_tail,
                          "C_nearest_reference_length": len(reference_c),
                          "A_sequence_comparison": "all_exact_published" if a_exact and all(a_exact) else "at_least_one_not_exact_published" if a_exact else "no_admissible_A",
                          "core_comparison": "exact_published_defined_span" if core_matches else "not_exact_in_exported_records" if core else "not_evaluated_for_uncertain_context",
                          "carrier_comparison": carrier_comparison,
                          "novelty_scope": "29_LAR_plus_393_Merrild_exports_only"})
    extra = pd.DataFrame(additions).fillna("")
    # Existing accession/BioSample columns are checked against the frozen
    # primary source metadata rather than silently replaced.
    if not (data.biosample.astype(str).reset_index(drop=True) == extra.biosample.astype(str)).all():
        raise ValueError("Candidate and primary-source BioSample metadata disagree")
    data = pd.concat([data, extra.drop(columns=[x for x in extra if x in data])], axis=1)
    primary = data[data.status == "family_bgc_candidate"].copy()
    uncertain = data[data.status != "family_bgc_candidate"].copy()
    for table in (data, primary, uncertain):
        table.to_csv(tables / (f"{output_prefix}_all_contexts.tsv" if table is data else f"{output_prefix}_BGC_catalogue.tsv" if table is primary else f"{output_prefix}_uncertain.tsv"), sep="\t", index=False)
    roles = pd.concat(role_tables, ignore_index=True)
    precursors = pd.concat(precursor_evidence, ignore_index=True)
    roles.to_csv(tables / f"{output_prefix}_role_evidence.tsv", sep="\t", index=False)
    precursors.to_csv(tables / f"{output_prefix}_precursor_evidence.tsv", sep="\t", index=False)
    role_lookup = {(r["query_id"], r["proposed_role"]): r for r in roles.to_dict("records")}
    a_lookup = precursors.set_index("query_id").to_dict("index")
    core_evidence = []
    for row in primary.to_dict("records"):
        for gene, assignment in zip(row["chosen_core_gene_ids"].split(";"), row["architecture"].split("-")):
            for role in assignment.split("+"):
                evidence = a_lookup[gene] if role == "A" else role_lookup[(gene, role)]
                core_evidence.append({"genome_id": row["genome_id"], "organism": row["organism"],
                                      "family": row["family"], "architecture_role": role, **evidence,
                                      "query_id": gene,
                                      "reference_id": evidence.get("reference_id", evidence.get("subject_id", "")),
                                      "reference_assignment_basis": evidence.get("reference_assignment_basis", "explicit_published_precursor_label"),
                                      "domain_panel": evidence.get("domain_panel", ""),
                                      "position_evidence": ("same_C_strand_selected_family_architecture" if evidence["strand"] == row["C_strand"]
                                                            else "opposite_A_with_intervening_C_strand_B1_B2_pair"),
                                      "function_status": "predicted_role_not_per_locus_function_experiment"})
    pd.DataFrame(core_evidence).fillna("").to_csv(tables / f"{output_prefix}_core_gene_evidence.tsv", sep="\t", index=False)
    roles[roles.locus_id.isin(primary.locus_id) & roles.proposed_role.isin(["ABC", "GNAT", "C11"])].to_csv(
        tables / f"{output_prefix}_auxiliary_gene_evidence.tsv", sep="\t", index=False)
    pd.concat(validation_tables, ignore_index=True).to_csv(tables / f"{output_prefix}_context_validation.tsv", sep="\t", index=False)
    proteomes.to_csv(tables / f"{output_prefix}_proteomes.tsv", sep="\t", index=False)
    hits.to_csv(tables / f"{output_prefix}_profile_hits.tsv", sep="\t", index=False)
    SeqIO.write({r.id: r for r in fasta_records}.values(), processed / f"{output_prefix}_BGC_genes.faa", "fasta")
    SeqIO.write(core_records, processed / f"{output_prefix}_BGC_cores.fna", "fasta")
    families = {}
    for family, group in primary.groupby("family"):
        a_hashes = {x for cell in group.precursor_sha256 for x in cell.split(";") if x}
        families[family] = {"source_locus_rows": len(group), "assembly_content_groups": group.assembly_contig_multiset_sha256.nunique(),
                            "C_variants": group.C_sha256.nunique(), "A_variants": len(a_hashes), "defined_core_variants": group.core_sha256.nunique(),
                            "C_variants_not_exact_published": group.loc[group.C_sequence_comparison != "exact_published", "C_sha256"].nunique(),
                            "C_start_annotation_differences_on_published_DNA": int((group.C_difference_interpretation == "start_annotation_difference_on_published_DNA").sum()),
                            "C_amino_acid_variants_with_unrepresented_DNA": group.loc[group.C_difference_interpretation == "different_C_amino_acids_and_unrepresented_DNA_span", "C_sha256"].nunique(),
                            "defined_cores_not_exact_published": group.loc[group.core_comparison != "exact_published_defined_span", "core_sha256"].nunique(),
                            "architectures": group.architecture.value_counts().to_dict(),
                            "sources": group.source.value_counts().to_dict()}
    summary = {"status": "completed_species_search_and_context_catalogue", "proteomes": len(proteomes),
               "proteins": int(proteomes.protein_count.sum()), "raw_proteome_bytes": int(proteomes.bytes.sum()),
               "source_counts": proteomes.source.value_counts().to_dict(), "nominated_proteins": hits.target_id.nunique(),
               "nominated_genomes": hits.genome_id.nunique(), "profile_hit_rows": len(hits),
               "context_status_counts": data.status.value_counts().to_dict(), "families": families,
               "primary_assembly_content_groups": primary.assembly_contig_multiset_sha256.nunique(),
               "primary_biosamples": primary.biosample.nunique(), "primary_organism_labels": primary.organism.nunique(),
               "published_export_records": len(published), "known_cross_source_GCA_base_overlaps": 3 if args.prefix == "diversity_species" else "see_selection_overlap_table",
               "limitations": ["Source rows are not all independent carriers; exact assembly-content hashes reconcile candidate duplicates.",
                               "Species references omit unclassified species and within-species accessory diversity.",
                               "Role assignments and supported BGCs are predictions, not demonstrated products.",
                               "Absent exact matches establish novelty only relative to the frozen exported collections."]}
    (tables / f"{output_prefix}_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
