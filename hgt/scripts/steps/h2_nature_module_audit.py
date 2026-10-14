#!/usr/bin/env python3
"""Inspect the four held Nature systems using exact alleles and native context."""
from collections import defaultdict
import json

from Bio import SeqIO

from prepare_h2 import H1, OUT, read, write
from h2_components import domains

genes = {r["gene_id"]: r for r in read(H1 / "H1_publication_genes.tsv")}
loci = {r["locus_id"]: r for r in read(H1 / "H1_carrier_catalogue.tsv")}
aa = {r.id: str(r.seq) for r in SeqIO.parse(H1 / "proteins.faa", "fasta")}
dna = {r.id: str(r.seq) for r in SeqIO.parse(H1 / "coding_sequences.fna", "fasta")}
references = {r["gene_id"]: r["family"] for r in read(OUT / "role_references.tsv")}
precursors = [g for g in genes.values() if g["locus_id"].startswith("Merrild_") and g["author_role"] in ["A", "A1"]]
verified = {r["gene_id"] for r in read(H1 / "H1_gene_genome_correspondence.tsv") if r["DNA_status"] == "exact_coding_DNA" and r["translation_status"] in ["exact_translation", "exact_with_bacterial_initiator_M"]}
pub_domains = domains(OUT / "publication_domains.domtbl")
module_rows, short_rows, context_domains, neighborhood = [], [], [], []
for leaf in read(OUT / "component_leaf_map.tsv"):
    if leaf["family"] != "uncertain_Nature_system_C":
        continue
    locus_id = leaf["publication_locus_ids"].split(";")[0]
    cgene = next(genes[g] for g in leaf["publication_gene_ids"].split(";") if genes[g]["locus_id"] == locus_id)
    directory = OUT / "native" / leaf["host_id"]
    native = {r["protein_id"]: r for r in read(directory / "proteins.tsv")}
    local = {key: r for key, r in native.items() if r["contig"] == leaf["contig"] and max(0, int(r["start_0based"]) - int(leaf["end_0based"]), int(leaf["start_0based"]) - int(r["end_0based"])) <= 5000}
    homologues = defaultdict(set)
    with (directory / "role_homology.tsv").open() as handle:
        for line in handle:
            f = line.strip().split("\t")
            if f[0] in local and min((int(f[5])-int(f[4])+1)/int(f[10]), (int(f[7])-int(f[6])+1)/int(f[11])) >= 0.7:
                homologues[references[f[1]]].add(f[0])
    hits = domains(directory / "domains.domtbl")
    for key, protein in native.items():
        distance = max(0, int(protein["start_0based"]) - int(leaf["end_0based"]), int(leaf["start_0based"]) - int(protein["end_0based"]))
        if protein["contig"] == leaf["contig"] and distance <= 12000:
            neighborhood.append({"Nature_locus": locus_id, "host_id": leaf["host_id"], "native_protein_id": key,
                "product": protein["product"], "protein_aa": protein["protein_aa"], "gene_start_0based": protein["start_0based"],
                "gene_end_0based": protein["end_0based"], "strand": protein["strand"], "distance_from_target_C_bp": distance,
                "C_start_0based": leaf["start_0based"], "C_end_0based": leaf["end_0based"],
                "domain_accessions": ";".join(sorted({x["accession"] for x in hits[key]})), "partial_location": protein["partial_location"]})
    for key, protein in local.items():
        for hit in hits[key]:
            context_domains.append({"Nature_locus": locus_id, "host_id": leaf["host_id"], "native_protein_id": key,
                "product": protein["product"], "protein_aa": protein["protein_aa"], "contig": protein["contig"],
                "gene_start_0based": protein["start_0based"], "gene_end_0based": protein["end_0based"],
                "C_start_0based": leaf["start_0based"], "C_end_0based": leaf["end_0based"],
                "reference_families": ";".join(k for k, v in homologues.items() if key in v), **hit})
    rre_c39 = sorted({x["accession"] for key in local for x in hits[key] if x["accession"] in ["PF05402", "PF13471", "PF03412"]})
    matches, matching_precursor_loci = [], set()
    for gene in genes.values():
        if gene["locus_id"] != locus_id or not 0 < int(gene["protein_aa"] or 0) <= 120:
            continue
        distance = max(0, int(gene["start_0based"]) - int(cgene["end_0based"]), int(cgene["start_0based"]) - int(gene["end_0based"]))
        if distance > 2000:
            continue
        same_aa = [p["gene_id"] for p in precursors if aa[p["gene_id"]] == aa[gene["gene_id"]]]
        same_dna = [p["gene_id"] for p in precursors if dna[p["gene_id"]] == dna[gene["gene_id"]]]
        matching_precursor_loci.update(genes[key]["locus_id"] for key in same_dna)
        if same_aa or same_dna:
            matches.append(gene["gene_id"])
        short_rows.append({"Nature_locus": locus_id, "gene_id": gene["gene_id"], "protein_aa": gene["protein_aa"],
            "distance_from_C_bp": distance, "amino_acid_sequence": aa[gene["gene_id"]], "coding_DNA": dna[gene["gene_id"]],
            "native_DNA_translation_confirmed": gene["gene_id"] in verified,
            "exact_Merrild_precursor_AA_matches": ";".join(same_aa), "exact_Merrild_precursor_DNA_matches": ";".join(same_dna),
            "role_status": "exact_published_precursor_allele" if same_dna else "short_ORF_candidate; function_not_assigned"})
    matching_C_loci = {key for key, r in loci.items() if key.startswith("Merrild_") and
                       dna[r["C_gene_id"]] == dna[cgene["gene_id"]] and aa[r["C_gene_id"]] == aa[cgene["gene_id"]]}
    paired = matching_C_loci & matching_precursor_loci
    module_rows.append({"Nature_locus": locus_id, "host_id": leaf["host_id"], "C_gene_id": cgene["gene_id"], "C_aa": cgene["protein_aa"],
        "exact_Merrild_C_DNA_AA_loci": ";".join(sorted(matching_C_loci)),
        "exact_precursor_and_C_DNA_in_same_Merrild_locus": ";".join(sorted(paired)),
        "paired_Merrild_analytical_lineages": ";".join(sorted({loci[key]["analytical_lineage"] for key in paired})),
        "C_domains": ";".join(sorted({x["accession"] for x in pub_domains[cgene["gene_id"]]})),
        "near_C_reference_homologues": json.dumps({k: sorted(v) for k, v in homologues.items()}, sort_keys=True),
        "near_C_RRE_C39_profile_hits": ";".join(rre_c39), "short_ORFs_matching_Merrild_precursors": ";".join(matches),
        "limit": "reference_and_profile_nondetection_is_not_proven_absence; GNAT_not_confirmed_triculamin_core; genomic_family_does_not_prove_product"})
write(OUT / "Nature_module_audit.tsv", module_rows)
write(OUT / "Nature_short_ORF_alleles.tsv", short_rows)
write(OUT / "Nature_native_context_domains.tsv", context_domains)
write(OUT / "Nature_native_neighborhood.tsv", neighborhood)
print(json.dumps({"Nature_systems": len(module_rows), "short_ORFs": len(short_rows), "exact_precursor_alleles": sum(bool(r["exact_Merrild_precursor_DNA_matches"]) for r in short_rows)}))
