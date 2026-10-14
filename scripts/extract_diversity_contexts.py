#!/usr/bin/env python3
"""Addressed annotated-genome context extraction for a completed screen shard."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import ExactPosition

from extract_maize_contexts import small_orfs
from profile_search import ROOT
from screen_diversity_proteomes import fetch_batch


def translated(dna, table=11, phase=0):
    dna = dna[phase:]
    # Partial annotation at contig boundaries may end in an incomplete codon.
    dna = dna[:len(dna)//3*3]
    protein = str(dna.translate(table=table)).removesuffix("*")
    if protein and str(dna[:3]) in CodonTable.unambiguous_dna_by_id[table].start_codons:
        protein = "M" + protein[1:]
    return protein


def read_genome(path, source):
    features, dna = [], {}
    if source == "RefSeq":
        with gzip.open(path, "rt") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                dna[record.id] = record
                for ordinal, feature in enumerate(record.features):
                    q = feature.qualifiers
                    if feature.type != "CDS" or "translation" not in q:
                        continue
                    features.append({"feature_id": q.get("locus_tag", [f"{record.id}_cds{ordinal}"])[0],
                                     "source_protein_id": q.get("protein_id", [""])[0], "contig": record.id,
                                     "start": int(feature.location.start), "end": int(feature.location.end),
                                     "strand": int(feature.location.strand), "product": q.get("product", [""])[0],
                                     "table": int(q.get("transl_table", [11])[0]),
                                     "phase": int(q.get("codon_start", [1])[0])-1,
                                     "partial_CDS": "partial" in q or any(not isinstance(p, ExactPosition)
                                          for part in feature.location.parts for p in (part.start, part.end)),
                                     "protein": q["translation"][0], "coding_dna": str(feature.extract(record.seq))})
    else:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip().split("\t")
                if f[2] != "CDS":
                    continue
                attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                features.append({"feature_id": attrs["ID"], "source_protein_id": attrs["ID"],
                                 "contig": f[0], "start": int(f[3])-1, "end": int(f[4]),
                                 "strand": 1 if f[6] == "+" else -1,
                                 "product": unquote(attrs.get("product", "")),
                                 "table": int(attrs.get("transl_table", 11)), "phase": int(f[7]),
                                 "partial_CDS": attrs.get("partial", "00") != "00"})
            dna = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        for feature in features:
            sequence = dna[feature["contig"]].seq[feature["start"]:feature["end"]]
            if feature["strand"] == -1:
                sequence = sequence.reverse_complement()
            feature["coding_dna"] = str(sequence)
            feature["protein"] = translated(sequence, feature["table"], feature["phase"])
    if not dna:
        raise ValueError(f"No assembly DNA in {path}")
    return features, dna


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen-name", required=True)
    args = parser.parse_args()
    tables = ROOT / "results/tables"
    hits = pd.read_csv(tables / f"{args.screen_name}_profile_hits.tsv", sep="\t")
    best = hits.sort_values("sequence_score", ascending=False).drop_duplicates("target_id")
    selected = pd.read_csv(ROOT / "config" / f"{args.screen_name}_selection.tsv", sep="\t", keep_default_na=False).set_index("genome_id")
    soil = pd.read_csv(ROOT / "config/diversity_soil_species.tsv", sep="\t", keep_default_na=False).set_index("Genome")
    candidates = {}
    for path in sorted((ROOT / "results/hmmer" / args.screen_name).glob("batch_*/candidates.faa")):
        candidates.update({r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")})
    if set(best.target_id) != set(candidates):
        raise ValueError("Nominated sequence set differs from HMM hits")
    protocol = json.loads((ROOT / "config/diversity_protocol.json").read_text())
    out = ROOT / "data/processed" / f"{args.screen_name}_contexts"
    out.mkdir(parents=True, exist_ok=True)
    raw = ROOT / "data/raw/diversity_contexts"
    raw.mkdir(parents=True, exist_ok=True)
    manifest = ROOT / "config" / f"{args.screen_name}_context_downloads.jsonl"
    pinned = {r["genome_id"]: r for r in (json.loads(x) for x in manifest.read_text().splitlines())} if manifest.exists() else {}
    loci, neighbors, validations = [], [], []
    groups = list(best.groupby("genome_id"))
    downloaded = {}
    with (out / "neighbors.faa").open("w") as neighbor_fasta, (out / "precursor_queries.faa").open("w") as precursor_fasta, (out / "windows.fna").open("w") as windows, manifest.open("a") as ledger:
        for index, (genome, genome_hits) in enumerate(groups, 1):
            if (index-1) % 50 == 0:
                requests = []
                for identifier, _ in groups[index-1:index+49]:
                    source = selected.loc[identifier, "source"]
                    url = selected.loc[identifier].get("context_url") or (selected.loc[identifier, "url"].replace("_protein.faa.gz", "_genomic.gbff.gz") if source == "RefSeq"
                           else soil.loc[identifier, "FTP_download"].replace("ftp://", "https://", 1))
                    requests.append({"genome_id": identifier, "source": source, "url": url, "compression": "gzip",
                                     "file_suffix": ".gbff.gz" if source == "RefSeq" else ".gff.gz"})
                downloaded = {r["genome_id"]: r for r in fetch_batch(requests, raw, workers=1)}
            record = downloaded[genome]
            source, url, digest = record["source"], record["url"], record["sha256"]
            path = ROOT / record["local_path"]
            if genome in pinned and (pinned[genome]["sha256"] != digest or pinned[genome]["url"] != url):
                raise ValueError(f"Changed context source: {genome}")
            if genome not in pinned:
                ledger.write(json.dumps({"genome_id": genome, "url": url, "sha256": digest,
                                         "local_path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                                         "retrieved_utc": datetime.now(timezone.utc).isoformat()}) + "\n")
                ledger.flush()
            features, contigs = read_genome(path, source)
            contig_hashes = sorted(hashlib.sha256(min(str(r.seq).upper(), str(r.seq.reverse_complement()).upper()).encode()).hexdigest()
                                   for r in contigs.values())
            assembly_hash = hashlib.sha256("\n".join(contig_hashes).encode()).hexdigest()
            seen = set()
            for hit in genome_hits.to_dict("records"):
                matches = [f for f in features if f["protein"] == candidates[hit["target_id"]]]
                validations.append({"target_id": hit["target_id"], "genome_id": genome,
                                    "annotated_exact_matches": len(matches),
                                    "status": "exact_C_annotation_match" if matches else "C_context_annotation_unresolved"})
                for focal in matches:
                    key = (focal["contig"], focal["start"], focal["end"])
                    if key in seen:
                        continue
                    seen.add(key)
                    locus = f"{genome}__{focal['contig']}__C{focal['start']}_{focal['end']}"
                    contig = contigs[focal["contig"]]
                    low = max(0, focal["start"]-protocol["context_flank_bp"])
                    high = min(len(contig), focal["end"]+protocol["context_flank_bp"])
                    c_dna_match = translated(Seq(focal["coding_dna"]), focal["table"], focal["phase"]) == focal["protein"]
                    loci.append({**hit, "locus_id": locus, "contig": focal["contig"], "C_start": focal["start"],
                                 "C_end": focal["end"], "C_strand": focal["strand"], "C_source_feature": focal["feature_id"],
                                 "C_DNA_translation_match": c_dna_match, "contig_bp": len(contig),
                                 "C_partial_CDS": focal["partial_CDS"],
                                 "window_start": low, "window_end": high,
                                 "left_edge_bp": focal["start"], "right_edge_bp": len(contig)-focal["end"],
                                 "source": source, "assembly": selected.loc[genome,"assembly"],
                                 "biosample": selected.loc[genome,"biosample"],
                                 "assembly_contig_multiset_sha256": assembly_hash})
                    SeqIO.write(SeqRecord(contig.seq[low:high], id=locus, description=""), windows, "fasta")
                    annotated_coordinates = set()
                    for ordinal, feature in enumerate(features):
                        if feature["contig"] != focal["contig"] or feature["start"] < low or feature["end"] > high:
                            continue
                        identifier = f"{locus}__cds{ordinal}"
                        row = {"locus_id": locus, "query_id": identifier, "source_feature_id": feature["feature_id"],
                               "source_protein_id": feature["source_protein_id"], "start": feature["start"],
                               "end": feature["end"], "strand": feature["strand"], "product": feature["product"],
                               "aa_length": len(feature["protein"]), "sequence_sha256": hashlib.sha256(feature["protein"].encode()).hexdigest(),
                               "origin": "annotated_CDS", "is_C": key == (feature["contig"], feature["start"], feature["end"])}
                        neighbors.append(row)
                        annotated_coordinates.add((feature["start"], feature["end"], feature["strand"]))
                        neighbor_fasta.write(f">{identifier}\n{feature['protein']}\n")
                        if len(feature["protein"]) <= 150:
                            precursor_fasta.write(f">{identifier}\n{feature['protein']}\n")
                    if hit["sequence_score"] >= protocol["strong_C_sequence_bits"] and hit["model_coverage"] >= protocol["strong_C_model_coverage"] and focal["table"] == 11:
                        records, orfs = small_orfs(contig.seq[low:high], low, locus)
                        for record, orf in zip(records, orfs):
                            strand = 1 if orf["strand"] == 1 else -1
                            if (orf["start_0based"], orf["end_0based_exclusive"], strand) in annotated_coordinates:
                                continue
                            neighbors.append({"locus_id": locus, "query_id": record.id, "source_feature_id": "",
                                              "source_protein_id": "", "start": orf["start_0based"], "end": orf["end_0based_exclusive"],
                                              "strand": strand, "product": "addressed_precursor_ORF_candidate", "aa_length": len(record),
                                              "sequence_sha256": hashlib.sha256(str(record.seq).encode()).hexdigest(),
                                              "origin": "addressed_sORF_prediction", "is_C": False})
                            SeqIO.write(record, precursor_fasta, "fasta")
            if index % 25 == 0:
                print(f"context genomes {index}/{best.genome_id.nunique()}; loci {len(loci)}", flush=True)
    pd.DataFrame(loci).to_csv(out / "loci.tsv", sep="\t", index=False)
    pd.DataFrame(neighbors).to_csv(out / "features.tsv", sep="\t", index=False)
    pd.DataFrame(validations).to_csv(out / "validation.tsv", sep="\t", index=False)
    print(json.dumps({"nominations": len(best), "context_loci": len(loci), "features": len(neighbors),
                      "unmapped_nominations": sum(not r["annotated_exact_matches"] for r in validations)}))


if __name__ == "__main__":
    main()
