#!/usr/bin/env python3
"""Export every direct hit's +/-10 kb context and explicit small-ORF candidates.

ORF enumeration is permissive and does not call precursors by itself. Coordinates
are always 0-based, half-open on the original MGnify contig.
"""
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from maize_genomes import read_genome
from profile_search import ROOT


def small_orfs(dna, offset, candidate):
    records, rows = [], []
    for strand, oriented in ((1, dna), (-1, dna.reverse_complement())):
        text = str(oriented)
        for frame in range(3):
            last_stop = frame - 3
            for stop in range(frame, len(text) - 2, 3):
                if text[stop:stop+3] not in ("TAA", "TGA", "TAG"):
                    continue
                first = max(last_stop + 3, stop - 300)
                for start in range(first, stop - 59, 3):
                    if text[start:start+3] not in ("ATG", "GTG", "TTG"):
                        continue
                    protein = "M" + str(Seq(text[start:stop]).translate(table=11))[1:]
                    low, high = (start, stop + 3) if strand == 1 else (len(text) - stop - 3, len(text) - start)
                    identifier = f"{candidate}__orf_{offset+low}_{offset+high}_{'p' if strand==1 else 'm'}"
                    records.append(SeqRecord(Seq(protein), id=identifier, description=""))
                    rows.append({"orf_id": identifier, "candidate_id": candidate,
                                 "start_0based": offset + low, "end_0based_exclusive": offset + high,
                                 "strand": strand, "aa_length": len(protein), "start_codon": text[start:start+3]})
                last_stop = stop
    return records, rows


def main():
    tables, processed = ROOT / "results/tables", ROOT / "data/processed"
    scores = pd.read_csv(tables / "maize_direct_profile_hits.tsv", sep="\t")
    hit_ids = set(scores.target_id)
    by_genome = defaultdict(set)
    for identifier in hit_ids:
        by_genome[identifier.split("_")[0]].add(identifier)
    loci, neighbors, all_orfs, orf_rows, genbank = [], [], [], [], []
    proteins, nucleotide_windows = {}, []
    for genome, genome_hits in sorted(by_genome.items()):
        features, sequences = read_genome(genome)
        lookup = {f["target_id"]: f for f in features}
        assert genome_hits <= set(lookup)
        for candidate in sorted(genome_hits):
            focal = lookup[candidate]
            start = max(0, focal["start_0based"] - 10000)
            end = min(focal["contig_bp"], focal["end_0based_exclusive"] + 10000)
            dna = sequences[focal["contig"]].seq[start:end]
            loci.append({k: v for k, v in focal.items() if k not in ("attributes", "protein")}
                        | {"window_start_0based": start, "window_end_0based_exclusive": end,
                           "C_left_edge_bp": focal["start_0based"],
                           "C_right_edge_bp": focal["contig_bp"] - focal["end_0based_exclusive"],
                           "aa_length": len(focal["protein"])})
            nucleotide_windows.append(SeqRecord(dna, id=candidate, description=f"{focal['contig']}:{start}:{end}"))
            record = SeqRecord(dna, id=candidate, name=candidate, description="Uncurated candidate window")
            record.annotations["molecule_type"] = "DNA"
            for f in features:
                if f["contig"] != focal["contig"] or f["start_0based"] >= end or f["end_0based_exclusive"] <= start:
                    continue
                neighbors.append({k: v for k, v in f.items() if k not in ("attributes", "protein")}
                                 | {"candidate_id": candidate, "aa_length": len(f["protein"]),
                                    "relative_start_to_C": f["start_0based"] - focal["start_0based"]})
                proteins[f["target_id"]] = SeqRecord(Seq(f["protein"]), id=f["target_id"], description=f["product"])
                # Only complete CDS inside the exported window receive translation.
                if f["start_0based"] >= start and f["end_0based_exclusive"] <= end:
                    record.features.append(SeqFeature(FeatureLocation(f["start_0based"] - start,
                        f["end_0based_exclusive"] - start, strand=1 if f["strand"] == "+" else -1),
                        type="CDS", qualifiers={"locus_tag": [f["target_id"]], "product": [f["product"]],
                                                "translation": [f["protein"]]}))
            genbank.append(record)
            orfs, rows = small_orfs(dna, start, candidate)
            all_orfs.extend(orfs)
            orf_rows.extend(rows)
    pd.DataFrame(loci).to_csv(tables / "maize_candidate_loci.tsv", sep="\t", index=False)
    pd.DataFrame(neighbors).to_csv(tables / "maize_candidate_neighbors.tsv", sep="\t", index=False)
    pd.DataFrame(orf_rows).to_csv(tables / "maize_candidate_small_orfs.tsv", sep="\t", index=False)
    SeqIO.write(proteins.values(), processed / "maize_candidate_neighbors.faa", "fasta")
    SeqIO.write(all_orfs, processed / "maize_candidate_small_orfs.faa", "fasta")
    SeqIO.write(nucleotide_windows, processed / "maize_candidate_windows.fna", "fasta")
    SeqIO.write(genbank, processed / "maize_candidate_windows.gbk", "genbank")
    print({"candidate_loci": len(loci), "neighbor_CDS": len(proteins), "small_ORF_candidates": len(all_orfs)})


if __name__ == "__main__":
    main()
