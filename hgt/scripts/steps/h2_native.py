#!/usr/bin/env python3
"""Extract frozen annotated proteins and run standard HMM/BLAST searches."""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_h2 import ROOT, SNAP, H1, OUT, read, write


def call(arguments, log):
    with log.open("w") as handle:
        subprocess.run([str(x) for x in arguments], stdout=handle, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--shards", type=int, default=16)
    args = parser.parse_args()
    native_root = Path(os.environ.get("HGT_GENOME_ROOT", str(ROOT)))
    bin_root = Path(os.environ["HGT_BIN"])
    small_bin = Path(os.environ["HGT_SMALL_BIN"])
    expected = {r["genome_id"]: r for r in json.loads((ROOT / "config/H2/H2_input_manifest.json").read_text())["genomes"]}
    aliases = {r["genome_id"]: r["host_id"] for r in read(OUT / "host_aliases.tsv")}
    windows = defaultdict(list)
    for row in read(H1 / "H1_locus_genome_correspondence.tsv"):
        if row["carrier_assignment"] == "source_provenance_supported" or row["sequence_evidence"] == "whole_export_exact":
            windows[row["genome_id"]].append(row)
    for row in read(OUT / "host_inputs.tsv")[args.shard::args.shards]:
        host_id = row["host_id"]
        directory = OUT / "native" / host_id
        directory.mkdir(parents=True, exist_ok=True)
        path = native_root / row["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected[row["genome_id"]]["sha256"]:
            raise ValueError(f"Frozen genome changed: {host_id}")
        opening = gzip.open if path.suffix == ".gz" else open
        with opening(path, "rt") as handle, warnings.catch_warnings():
            warnings.simplefilter("ignore")
            records = list(SeqIO.parse(handle, "genbank"))
        metadata, proteins = [], []
        sequence_info = []
        genome_windows = windows[row["genome_id"]]
        for record in records:
            sequence = str(record.seq).upper()
            reverse = str(Seq(sequence).reverse_complement())
            canonical = min(sequence, reverse)
            sequence_info.append({"contig": record.id, "bp": len(sequence), "canonical_contig_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
                                  "native_to_canonical_reverse": sequence != canonical})
            for feature_index, feature in enumerate(record.features):
                if feature.type != "CDS" or not feature.qualifiers.get("translation"):
                    continue
                protein = "".join(feature.qualifiers["translation"]).rstrip("*")
                if not protein or "*" in protein:
                    continue
                protein_id = f"{host_id}_P{len(proteins) + 1:05}"
                start, end = int(feature.location.start), int(feature.location.end)
                overlaps = any(w["contig"] == record.id and start < int(w["export_end_0based"]) and end > int(w["export_start_0based"]) for w in genome_windows)
                metadata.append({"protein_id": protein_id, "host_id": host_id, "source": "native_GenBank_annotation",
                                 "contig": record.id, "start_0based": start, "end_0based": end, "strand": feature.location.strand,
                                 "native_accession": ";".join(feature.qualifiers.get("protein_id", [])),
                                 "product": ";".join(feature.qualifiers.get("product", [])), "protein_aa": len(protein),
                                 "overlaps_mapped_export_window": overlaps,
                                 "partial_location": any(type(p).__name__ != "ExactPosition" for p in [feature.location.start, feature.location.end])})
                proteins.append(SeqRecord(Seq(protein), id=protein_id, description=""))
        SeqIO.write((SeqRecord(r.seq, id=r.id, description="") for r in records), directory / "genome.fna", "fasta")
        if not proteins:
            call([small_bin / "prodigal", "-i", directory / "genome.fna", "-a", directory / "prodigal.faa", "-o", directory / "prodigal.gff", "-f", "gff", "-p", "meta", "-q"], directory / "prodigal.log")
            for record in SeqIO.parse(directory / "prodigal.faa", "fasta"):
                fields = record.description.split(" # ")
                contig = record.id.rsplit("_", 1)[0]
                start, end, strand = int(fields[1]) - 1, int(fields[2]), int(fields[3])
                protein = str(record.seq).rstrip("*")
                protein_id = f"{host_id}_P{len(proteins) + 1:05}"
                overlaps = any(w["contig"] == contig and start < int(w["export_end_0based"]) and end > int(w["export_start_0based"]) for w in genome_windows)
                metadata.append({"protein_id": protein_id, "host_id": host_id, "source": "Prodigal_meta_for_unannotated_input",
                                 "contig": contig, "start_0based": start, "end_0based": end, "strand": strand,
                                 "native_accession": "", "product": "not_functionally_annotated", "protein_aa": len(protein),
                                 "overlaps_mapped_export_window": overlaps, "partial_location": "partial=00" not in fields[-1]})
                proteins.append(SeqRecord(Seq(protein), id=protein_id, description=""))
        write(directory / "proteins.tsv", metadata)
        write(directory / "contigs.tsv", sequence_info)
        SeqIO.write(proteins, directory / "proteins.faa", "fasta")
        for name, hmm in [("markers", "Bacteria.hmm"), ("domains", "BGC_domains.hmm")]:
            call([bin_root / "hmmsearch", "--cpu", "2", "--cut_ga", "--tblout", directory / (name + ".tbl"),
                  "--domtblout", directory / (name + ".domtbl"), SNAP / "references" / hmm, directory / "proteins.faa"], directory / (name + ".log"))
        call([bin_root / "blastp", "-query", directory / "proteins.faa", "-subject", OUT / "role_references.faa",
              "-evalue", "1e-5", "-qcov_hsp_perc", "70", "-num_threads", "2", "-max_target_seqs", "20",
              "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen",
              "-out", directory / "role_homology.tsv"], directory / "blast.log")
        (directory / "complete.json").write_text(json.dumps({"host_id": host_id, "native_proteins": len(proteins), "genome_id": row["genome_id"], "input_sha256": expected[row["genome_id"]]["sha256"]}, indent=2) + "\n")
        print(f"{host_id}: {len(proteins)} proteins; HMM and reference searches complete", flush=True)


if __name__ == "__main__":
    main()
