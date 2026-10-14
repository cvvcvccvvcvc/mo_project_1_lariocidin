#!/usr/bin/env python3
"""Collect HMMER hits, mask duplicates/BGC overlaps, and prepare marker inputs."""
import argparse
from collections import Counter, defaultdict
import math
import hashlib
import json
import statistics
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_h2 import ROOT, SNAP, H1, OUT, read, write


def hmm_models(path):
    models = []
    for block in path.read_text().split("//"):
        fields = {}
        for line in block.splitlines():
            if line.startswith(("NAME ", "ACC ", "LENG ")):
                key, value = line.split(maxsplit=1)
                fields[key] = value.strip()
        if fields:
            models.append({"name": fields["NAME"], "accession": fields["ACC"], "marker_id": fields["ACC"].split(".")[0], "hmm_length": int(fields["LENG"])})
    return models


def table_hits(path):
    result = defaultdict(dict)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            result[fields[2]][fields[0]] = {"protein_id": fields[0], "model": fields[2], "E_value": fields[4], "bitscore": float(fields[5])}
    return result


def taxonomy():
    rows = {}
    for path in (SNAP / "references").glob("NCBI_taxonomy_*.xml"):
        for taxon in ET.parse(path).getroot().findall("Taxon"):
            ranks = {x.findtext("Rank"): x.findtext("ScientificName") for x in taxon.findall("LineageEx/Taxon")}
            ranks[taxon.findtext("Rank")] = taxon.findtext("ScientificName")
            rows[taxon.findtext("TaxId")] = {"taxonomy_name": taxon.findtext("ScientificName"), "phylum": ranks.get("phylum", "unresolved"),
                                           "class": ranks.get("class", "unresolved"), "order": ranks.get("order", "unresolved"),
                                           "family": ranks.get("family", "unresolved"), "genus": ranks.get("genus", "unresolved")}
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--actino", action="store_true", help="Additional published 138-marker Actinomycetota panel")
    args = parser.parse_args()
    destination = OUT / "actino" if args.actino else OUT
    destination.mkdir(exist_ok=True)
    tax = taxonomy()
    hosts = read(OUT / "host_inputs.tsv")
    if args.actino:
        hosts = [r for r in hosts if tax.get(r["taxid"].split(";")[0], {}).get("phylum") == "Actinomycetota"]
    models = hmm_models(SNAP / ("references/Actinobacteria.hmm" if args.actino else "references/Bacteria.hmm"))
    cutoff = math.ceil(0.8 * len(models))
    lengths = defaultdict(list)
    native, hits = {}, {}
    for host in hosts:
        directory = OUT / "native" / host["host_id"]
        if not (directory / "complete.json").exists():
            raise ValueError(f"Native stage incomplete: {host['host_id']}")
        metadata = {r["protein_id"]: r for r in read(directory / "proteins.tsv")}
        native[host["host_id"]] = metadata
        hits[host["host_id"]] = table_hits(directory / ("actino_markers.tbl" if args.actino else "markers.tbl"))
        for name, targets in hits[host["host_id"]].items():
            lengths[name].extend(int(metadata[key]["protein_aa"]) for key in targets)
    medians = {name: statistics.median(values) for name, values in lengths.items()}
    selected, audit, quality = defaultdict(list), [], []
    for host in hosts:
        host_id = host["host_id"]
        sequences = {r.id: str(r.seq) for r in SeqIO.parse(OUT / "native" / host_id / "proteins.faa", "fasta")}
        states = Counter()
        for model in models:
            targets = hits[host_id].get(model["name"], {})
            target = next(iter(targets), "")
            status = "missing" if not targets else "multicopy_masked" if len(targets) > 1 else \
                "BGC_window_overlap_masked" if native[host_id][target]["overlaps_mapped_export_window"] == "True" else \
                "length_outlier_masked" if not 0.8 * medians[model["name"]] <= int(native[host_id][target]["protein_aa"]) <= 1.2 * medians[model["name"]] else \
                "partial_annotation_masked" if native[host_id][target]["partial_location"] == "True" else "usable_single_copy"
            states[status] += 1
            audit.append({"host_id": host_id, "marker_id": model["marker_id"], "model": model["name"], "raw_hits": len(targets),
                          "status": status, "protein_ids": ";".join(targets), "median_length_aa": medians.get(model["name"], ""),
                          "native_length_aa": native[host_id][target]["protein_aa"] if target else ""})
            if status == "usable_single_copy":
                selected[model["marker_id"]].append(SeqRecord(Seq(sequences[target]), id=host_id, description=""))
        quality.append({**host, **tax.get(host["taxid"].split(";")[0], {"taxonomy_name": host["organism"], "phylum": "unresolved", "class": "unresolved", "order": "unresolved", "family": "unresolved", "genus": "unresolved"}),
                        "usable_markers": states["usable_single_copy"], "usable_fraction": states["usable_single_copy"] / len(models),
                        "multicopy_markers": states["multicopy_masked"], "BGC_overlap_markers": states["BGC_window_overlap_masked"],
                        "length_outlier_markers": states["length_outlier_masked"], "partial_markers": states["partial_annotation_masked"],
                        "missing_markers": states["missing"], "host_tree_admission": states["usable_single_copy"] >= cutoff})
    write(destination / "marker_hit_audit.tsv", audit)
    write(destination / "marker_genome_QC.tsv", quality)
    marker_directory = destination / "markers"
    marker_directory.mkdir(exist_ok=True)
    for model in models:
        model["median_protein_length"] = medians.get(model["name"], "")
        model["usable_genomes"] = len(selected[model["marker_id"]])
        SeqIO.write(selected[model["marker_id"]], marker_directory / (model["marker_id"] + ".faa"), "fasta")
    write(destination / "marker_models.tsv", models)
    aliases = {r["genome_id"]: r["host_id"] for r in read(OUT / "host_aliases.tsv")}
    catalogue = {r["locus_id"]: r for r in read(H1 / "H1_carrier_catalogue.tsv")}
    targets = defaultdict(set)
    for row in read(H1 / "H1_locus_genome_correspondence.tsv"):
        if row["carrier_assignment"] == "source_provenance_supported" or row["sequence_evidence"] == "whole_export_exact":
            targets[catalogue[row["locus_id"]]["analytical_lineage"]].add(aliases[row["genome_id"]])
    wanted = {lineage: set(targets[lineage]) for lineage in ["canonical", "noncanonical"]}
    for row in read(H1 / "H1_relative_panel.tsv"):
        for lineage in wanted:
            if any(aliases.get(x) in targets[lineage] for x in row["focal_genome_inputs"].split(";")):
                if row["downloaded_genome_input"] in aliases:
                    wanted[lineage].add(aliases[row["downloaded_genome_input"]])
    admitted = {r["host_id"] for r in quality if r["host_tree_admission"]}
    panels = []
    mags = set(json.loads((ROOT / "config/H2_MAG_decision.json").read_text())["hosts"])
    for lineage, host_ids in wanted.items():
        host_ids &= {r["host_id"] for r in hosts}
        for host_id in sorted(host_ids):
            panels.append({"lineage": lineage, "host_id": host_id, "role": "observed_locus_carrier" if host_id in targets[lineage] else "taxonomic_relative",
                           "marker_QC_pass": host_id in admitted, "MAG": host_id in mags,
                           "primary_admission": host_id in admitted, "no_MAG_admission": host_id in admitted and host_id not in mags})
    write(destination / "host_panels.tsv", panels)
    summary = {"sequence_representatives": len(hosts), "QC_pass": len(admitted), "QC_fail": len(hosts) - len(admitted),
               "marker_status_counts": dict(Counter(r["status"] for r in audit)),
               "panels": {lineage: {"wanted": len(ids), "admitted": len(ids & admitted)} for lineage, ids in wanted.items()},
               "MAG_QC": {r["host_id"]: r["usable_markers"] for r in quality if r["host_id"] in mags},
               "panel_models": len(models), "admission_cutoff": cutoff, "panel_role": "Actinomycetota_sensitivity" if args.actino else "Bacteria74_primary"}
    (destination / "marker_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
