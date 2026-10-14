#!/usr/bin/env python3
"""Select source recovery and taxonomic relatives from frozen primary metadata."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import re
from pathlib import Path

from prepare_h1 import ROOT, H1, TABLES, read_tsv, write_tsv
from fetch_h1 import Snapshot, fetch_genome, atomic_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--basis", type=Path, default=H1 / "relative_selection_basis")
    args = parser.parse_args()
    catalogue = read_tsv(args.basis / "H1_carrier_catalogue.tsv")
    quality = read_tsv(args.basis / "H1_genome_quality.tsv")
    source_ids = {x for r in catalogue for x in r["primary_genome_inputs"].split(";") if x}
    source = [r for r in quality if r["genome_id"] in source_ids]
    existing_accessions = {r[k] for r in quality for k in ["requested_assembly_accession", "sequence_provider_assembly"] if r[k]}
    existing_biosamples = {r["biosample"] for r in quality if r["biosample"]}
    primary_sample_ids = {r["biosample"] for r in source if r["biosample"]}
    unresolved = [r for r in catalogue if not r["primary_genome_inputs"]]
    unresolved_samples = {x for r in unresolved for x in re.findall(r"SAM[END][A-Z]*\d+", r["archive_member"])}
    unresolved_npdc = {x for r in unresolved for x in re.findall(r"NPDC\d+", r["archive_member"])}
    metadata, recovery = {}, {}
    file = H1 / "reference_metadata/refseq_bacterial_assembly_summary.txt"
    with file.open() as handle:
        header = next(line[1:].strip().split("\t") for line in handle if line.startswith("#assembly_accession"))
        for row in csv.DictReader(handle, fieldnames=header, delimiter="\t"):
            accession = row["assembly_accession"]
            if row["biosample"] in unresolved_samples or any(x in row["infraspecific_name"] + ";" + row["isolate"] for x in unresolved_npdc):
                recovery[accession] = row
            if accession in existing_accessions or row["gbrs_paired_asm"] in existing_accessions:
                metadata[accession] = row
    for genome in source:
        candidates = [r for r in metadata.values() if r["assembly_accession"] in
                      [genome["requested_assembly_accession"], genome["sequence_provider_assembly"]] or
                      r["gbrs_paired_asm"] in [genome["requested_assembly_accession"], genome["sequence_provider_assembly"]]]
        if candidates:
            genome["species_taxid"] = candidates[0]["species_taxid"]
            genome["strain"] = genome["strain"] or candidates[0]["infraspecific_name"]
    species_sources = defaultdict(list)
    for genome in source:
        if genome["species_taxid"]:
            species_sources[genome["species_taxid"]].append(genome)
    pool = read_tsv(TABLES / "H1_relative_metadata_pool.tsv")
    by_species = defaultdict(list)
    for row in pool:
        by_species[row["species_taxid"]].append(row)
    choices, selection_gaps = [], []
    def rank(row):
        return (row["relation_to_type_material"] in ["na", ""], row["assembly_level"] != "Complete Genome",
                int(row["contig_count"]) if row["contig_count"].isdigit() else 10**9, row["assembly_accession"])
    for taxid, focal in sorted(species_sources.items()):
        used_samples, chosen = set(primary_sample_ids), []
        focal_strains = {re.sub(r"[^A-Z0-9]", "", g["strain"].upper().removeprefix("STRAIN=")) for g in focal if g["strain"]}
        for row in sorted(by_species[taxid], key=rank):
            strain = re.sub(r"[^A-Z0-9]", "", row["infraspecific_name"].upper().removeprefix("STRAIN="))
            if row["assembly_accession"] in existing_accessions or row["gbrs_paired_asm"] in existing_accessions or row["biosample"] in used_samples:
                continue
            if strain and strain != "NA" and strain in focal_strains:
                continue
            if "metagenome" in row["excluded_from_refseq"].lower() or "MAG" in row["isolate"] or row["genome_rep"] != "Full" or not row["biosample"]:
                continue
            chosen.append(row)
            used_samples.add(row["biosample"])
            choices.append({**row, "focal_species_taxid": taxid, "focal_genome_inputs": ";".join(sorted(g["genome_id"] for g in focal)),
                            "selection_reason": "same_NCBI_species; distinct_sample_and_reported_strain; prefer_type_material_then_complete_then_fewer_contigs",
                            "role": "vertical_inheritance_and_distribution_control",
                            "closeness_limit": "taxonomy_based_initial_panel; genetic_nearest_neighbour_not_yet_established",
                            "BGC_status": "not_assessed; exact_match_absence_is_not_BGC_absence"})
            if len(chosen) == 2:
                break
        if len(chosen) < 2:
            selection_gaps.append({"species_taxid": taxid, "focal_genome_inputs": ";".join(g["genome_id"] for g in focal),
                                   "selected_distinct_relatives": len(chosen),
                                   "reason": "fewer_than_two_eligible_same_species_inputs_in_frozen_pool; genus_nearest_requires_marker_tree"})
    ipg = read_tsv(H1 / "ncbi/unlinked_C_ipg.tsv")
    genes = {r["gene_id"]: r for r in read_tsv(TABLES / "H1_publication_genes.tsv")}
    specific_ids = {genes[r["C_gene_id"]]["protein_id"] for r in catalogue if r["C_gene_id"] and
                    not genes[r["C_gene_id"]]["protein_id"].startswith("WP_")}
    ipg_accessions = sorted({r["Assembly"] for r in ipg if r["Protein"] in specific_ids and r["Assembly"] and r["Assembly"] not in existing_accessions})
    snapshot = Snapshot()
    if ipg_accessions:
        query = " OR ".join(f'"{x}"[Assembly Accession]' for x in ipg_accessions)
        found = snapshot.json("H1_supplement_IPG_assembly_search", "esearch.fcgi", {"db": "assembly", "term": query, "retmax": 100, "retmode": "json"})
        ids = found["esearchresult"]["idlist"]
        if ids:
            summaries = snapshot.json("H1_supplement_IPG_assembly_summaries", "esummary.fcgi", {"db": "assembly", "id": ",".join(ids), "retmode": "json"})
            for uid in summaries["result"]["uids"]:
                row = summaries["result"][uid]
                synonym = row.get("synonym", {})
                requested = [x for x in ipg_accessions if x in [row["assemblyaccession"], synonym.get("genbank"), synonym.get("refseq")]]
                if requested:
                    recovery[row["assemblyaccession"]] = row
    resolved_ipg = {r.get("assemblyaccession", k) for k, r in recovery.items()}
    resolved_ipg.update(x for row in recovery.values() for x in row.get("synonym", {}).values())
    missing_ipg = sorted(set(ipg_accessions) - resolved_ipg)
    write_tsv(TABLES / "H1_relative_selection.tsv", choices)
    if selection_gaps:
        write_tsv(TABLES / "H1_relative_selection_gaps.tsv", selection_gaps)
    metadata_all = {r["assembly_accession"]: r for r in choices}
    metadata_all.update(recovery)
    metadata_all = {k: v for k, v in metadata_all.items() if k not in existing_accessions}
    report = {"primary_genome_inputs": len(source), "focal_species_with_taxid": len(species_sources),
              "focal_genomes_without_species_taxid": [g["genome_id"] for g in source if not g["species_taxid"]],
              "relative_choices": len(choices), "species_with_relative_selection_gaps": len(selection_gaps),
              "source_recovery_candidates": len(recovery), "missing_specific_IPG_assemblies": missing_ipg,
              "new_genome_inputs": len(metadata_all)}
    atomic_json(TABLES / "H1_supplement_selection_summary.json", report)
    atomic_json(H1 / "ncbi/H1_supplement_candidate_metadata.json", metadata_all)
    print(json.dumps(report), flush=True)
    if not args.download:
        return
    supplemental, errors = [], []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(fetch_genome, pair) for pair in sorted(metadata_all.items())]
        for index, future in enumerate(as_completed(futures), 1):
            accession, downloaded, error = future.result()
            if error:
                errors.append({"assembly_accession": accession, "error": error})
                continue
            path, url = downloaded
            snapshot.record("genome__" + accession, path, url)
            row = metadata_all[accession]
            supplemental.append({"genome_id": accession, "assembly_accession": accession, "path": str(path.relative_to(ROOT)),
                                 "input_scope": "source_recovery_candidate" if accession in recovery else "taxonomic_relative",
                                 "biosample": row.get("biosample", row.get("biosampleaccn", "")),
                                 "organism": row.get("organism_name", row.get("organism", "")),
                                 "assembly_level": row.get("assembly_level", row.get("assemblystatus", "")),
                                 "source_strain": row.get("infraspecific_name", ""), "taxid": row.get("taxid", "")})
            if index % 50 == 0:
                print(f"Supplemental genomes {index}/{len(futures)}; errors {len(errors)}", flush=True)
    write_tsv(TABLES / "H1_supplement_assembly_inputs.tsv", sorted(supplemental, key=lambda r: r["genome_id"]))
    atomic_json(TABLES / "H1_supplement_download_errors.json", errors)
    print(json.dumps({"downloaded_inputs": len(supplemental), "errors": len(errors)}), flush=True)


if __name__ == "__main__":
    main()
