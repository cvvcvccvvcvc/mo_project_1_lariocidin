#!/usr/bin/env python3
"""Freeze metadata and reviewed comparators for the species-wide diversity stage.

This prepares inputs, not a completed genome screen. Existing S3/S5 snapshots
and results are left intact. Reuse the project's pinned download cache.
"""
from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
from collections import Counter
from pathlib import Path

from source_cache import Cache, ROOT, result_uids, assembly_summaries
from fetch_research_inputs import restore_snapshots


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", default="2026-10-12")
    args = parser.parse_args()
    if (ROOT / "inputs/archives/diversity_metadata_snapshots.zip").exists():
        restore_snapshots(diversity=True)
    raw = ROOT / "data/raw" / f"diversity_{args.snapshot}"
    cache = Cache(raw, ROOT / "inputs/provenance/diversity_download_manifest.json")
    soil_base = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/soil/v1.0"
    soil = cache.download("soil_metadata", soil_base + "/genomes-all_metadata.tsv", ".tsv")
    cache.download("soil_readme", soil_base + "/README_v1.0.txt", ".txt")
    cache.download("soil_protein_listing", soil_base + "/protein_catalogue/", ".html")
    query = urllib.parse.quote("(reviewed:true) AND (taxonomy_id:2)")
    cache.download("reviewed_bacterial_proteins", "https://rest.uniprot.org/uniprotkb/stream?"
                   f"query={query}&format=fasta&compressed=true", ".faa.gz")
    cache.download("reviewed_bacterial_count", "https://rest.uniprot.org/uniprotkb/search?"
                   f"query={query}&format=tsv&fields=accession&size=1", ".tsv")
    novelty_cache = Cache(raw, ROOT / "inputs/provenance/diversity_novelty_download_manifest.json")
    # Exact assembly versions of the three historical S3 findings, not a new
    # taxon selection. Preserve each primary BioSample/paired-assembly link.
    for accession in ("GCA_002582335.1", "GCA_002569665.1", "GCA_002572225.1"):
        result = novelty_cache.get_json(accession + "__search", "esearch.fcgi",
                    {"db": "assembly", "term": accession + "[Assembly Accession]", "retmode": "json"})
        summaries = assembly_summaries(novelty_cache, accession, result_uids(result))
        if len(summaries) != 1 or summaries[0]["synonym"]["genbank"] != accession:
            raise ValueError(f"Ambiguous historical assembly: {accession}")
        base = summaries[0]["ftppath_refseq"].replace("ftp://", "https://")
        name = base.rsplit("/", 1)[1]
        novelty_cache.download(accession + "__refseq_proteins", base + "/" + name + "_protein.faa.gz", ".faa.gz")
    novelty_cache.get_json("bacillus_C_protein_summaries", "esummary.fcgi",
                  {"db": "protein", "id": "WP_176544197.1,WP_176547010.1", "retmode": "json"})

    reports, token, page = [], "", 0
    while True:
        params = {"filters.reference_only": "true", "filters.assembly_source": "refseq",
                  "page_size": 1000}
        if token:
            params["page_token"] = token
        url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/2/dataset_report?" + urllib.parse.urlencode(params)
        path = cache.download(f"refseq_reference_page_{page:03d}", url, ".json")
        payload = json.loads(path.read_text())
        reports.extend(payload.get("reports", []))
        token = payload.get("next_page_token", "")
        if not token:
            break
        page += 1
    if not reports or len({r["accession"] for r in reports}) != len(reports):
        raise ValueError("Reference report is empty or has repeated accessions")

    selected = []
    for report in reports:
        info, organism = report["assembly_info"], report["organism"]
        if info["assembly_status"] != "current" or not report.get("annotation_info"):
            continue
        selected.append({"source": "RefSeq", "assembly": report["accession"],
                         "paired_assembly": report.get("paired_accession", ""),
                         "biosample": info.get("biosample", {}).get("accession", ""),
                         "tax_id": organism["tax_id"], "organism": organism["organism_name"],
                         "assembly_name": info["assembly_name"],
                         "annotation": report["annotation_info"].get("name", ""),
                         "selection": "current_annotated_bacterial_refseq_reference"})
    with (ROOT / "inputs/selections/species/diversity_refseq_species.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    with soil.open() as handle:
        genomes = list(csv.DictReader(handle, delimiter="\t"))
    selected_soil = [r for r in genomes if r["Genome"] == r["Species_rep"]
                     and r["Lineage"].startswith("d__Bacteria;")
                     and float(r["Completeness"]) >= 50 and float(r["Contamination"]) < 10]
    with (ROOT / "inputs/selections/species/diversity_soil_species.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(genomes[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_soil)
    by_genbank_base = {r["paired_assembly"].split(".")[0]: r for r in selected if r["paired_assembly"]}
    overlap = []
    for row in selected_soil:
        reference = by_genbank_base.get(row["Genome_accession"].split(".")[0])
        if reference:
            overlap.append({"soil_genome": row["Genome"], "soil_assembly": row["Genome_accession"],
                            "soil_biosample": row["Sample_accession"], "refseq_assembly": reference["assembly"],
                            "paired_genbank_assembly": reference["paired_assembly"],
                            "refseq_biosample": reference["biosample"],
                            "status": "same_GCA_base_verify_version_before_collapsing_occurrences"})
    with (ROOT / "results/tables/diversity_source_overlap.tsv").open("w") as handle:
        fields = ["soil_genome", "soil_assembly", "soil_biosample", "refseq_assembly",
                  "paired_genbank_assembly", "refseq_biosample", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(overlap)
    summary = {"stage": "source_selection_not_search_result", "snapshot": args.snapshot,
               "refseq_report_total": len(reports), "refseq_selected": len(selected),
               "soil_total": len(genomes), "soil_selected_species_representatives": len(selected_soil),
               "soil_selected_source_labels": dict(Counter(r["Genome_type"] for r in selected_soil)),
               "cross_source_shared_GCA_base": len(overlap),
               "soil_selected_phyla": dict(Counter(r["Lineage"].split(";")[1] for r in selected_soil)),
               "limitations": ["RefSeq references omit undefined species and within-species accessory diversity.",
                               "Soil MAG/isolate labels require primary provenance review for any finding.",
                               "The 50/10 numerical MAG filter does not establish MIMAG high quality.",
                               "Cross-source assembly aliases must be reconciled before counting carriers."]}
    out = ROOT / "results/tables/diversity_source_selection.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("limitations", "soil_selected_phyla")}))


if __name__ == "__main__":
    main()
