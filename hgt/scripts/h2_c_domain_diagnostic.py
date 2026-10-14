#!/usr/bin/env python3
"""Domain-family diagnostic, separate from the two primary BGC histories."""
from collections import defaultdict
import json

from Bio import SeqIO

from prepare_h2 import ROOT, H1, OUT, read, write
from h2_components import domains
from h2_prepare_trees import save_fasta

pub = {r.id: str(r.seq) for r in SeqIO.parse(H1 / "proteins.faa", "fasta")}
pub_domains = domains(OUT / "publication_domains.domtbl")
selected = {(r["family"], r["tip_id"]): r["publication_gene_id"] for r in read(OUT / "C_annotation_variants.tsv") if r["selected_model"] == "True"}
qc = {r["host_id"]: r for r in read(OUT / "marker_genome_QC.tsv")}
records, aliases, skipped = {}, [], []
by_sequence = {}
for leaf in read(OUT / "component_leaf_map.tsv"):
    if not leaf["family"].endswith("_C"):
        continue
    gene = selected[(leaf["family"], leaf["tip_id"])]
    hits = [r for r in pub_domains[gene] if r["accession"] == "PF00733"]
    if len(hits) != 1:
        skipped.append({"tip_id": leaf["tip_id"], "reason": "PF00733_domain_not_single; no_arbitrary_domain_choice"})
        continue
    hit = hits[0]
    seq = pub[gene][hit["start_aa"]:hit["end_aa"]]
    key = by_sequence.setdefault(seq, leaf["tip_id"])
    records[key] = seq
    aliases.append({"diagnostic_tip": key, "original_tip": leaf["tip_id"], "source": leaf["family"], "host_id": leaf["host_id"],
                    "genus": qc[leaf["host_id"]]["genus"], "protein_id": gene, "domain_start_aa": hit["start_aa"], "domain_end_aa": hit["end_aa"],
                    "basis": "DNA_verified_published_C; exact_domain_sequence_aliases_only_for_this_diagnostic"})
background = defaultdict(list)
for row in read(OUT / "C_metabolic_background_candidates.tsv"):
    background[row["host_id"]].append(row)
genus_best = {}
for host, candidates in background.items():
    directory = OUT / "native" / host
    hits = domains(directory / "domains.domtbl")
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(directory / "proteins.faa", "fasta")}
    genus = qc[host]["genus"]
    if genus == "unresolved":
        genus = host
    for row in candidates:
        match = [r for r in hits[row["native_protein_id"]] if r["accession"] == "PF00733"]
        if len(match) != 1:
            continue
        hit = match[0]
        ranking = (hit["bitscore"], host, row["native_protein_id"])
        if genus not in genus_best or ranking > genus_best[genus][0]:
            genus_best[genus] = (ranking, row, hit, seqs[row["native_protein_id"]][hit["start_aa"]:hit["end_aa"]])
for genus, (_, row, hit, seq) in sorted(genus_best.items()):
    key = row["host_id"] + "_BG" + row["native_protein_id"]
    records[key] = seq
    aliases.append({"diagnostic_tip": key, "original_tip": key, "source": "outside_mapped_loci_background", "host_id": row["host_id"],
                    "genus": genus, "protein_id": row["native_protein_id"], "domain_start_aa": hit["start_aa"], "domain_end_aa": hit["end_aa"],
                    "basis": "one_highest_domain_score_per_genus; background_is_not_a_proven_metabolic_negative_control"})
name = "C_domain_background_diagnostic"
output = OUT / "tree_inputs" / (name + ".faa")
save_fasta(output, records)
write(OUT / "C_domain_diagnostic_aliases.tsv", aliases)
if skipped:
    write(OUT / "C_domain_diagnostic_skipped.tsv", skipped)
jobs = read(OUT / "tree_jobs.tsv")
assert all(r["name"] != name for r in jobs), "Do not overwrite a submitted diagnostic"
index = sum(r["kind"] == "component" for r in jobs)
jobs.append({"kind": "component", "name": name, "input_path": str(output.relative_to(ROOT)), "partition_path": "", "filter": "gappyout", "dataset_role": "PF00733_domain_family_diagnostic_only; no_pooled_primary_BGC_history"})
write(OUT / "tree_jobs.tsv", jobs)
print(json.dumps({"component_index": index, "diagnostic_tips": len(records), "C_domain_sequences": len(by_sequence), "background_genera": len(genus_best), "skipped_C_copies": len(skipped)}))
