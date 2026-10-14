#!/usr/bin/env python3
"""Screen annotated assembly proteomes in bounded batches with existing HMMs.

The default small deterministic sample measures resource costs. Set --limit 0
only for the full fixed species selection. No translated-genome search is used.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import resource
import shutil
import subprocess
import time
from pathlib import Path
from datetime import datetime, timezone

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

from profile_search import ROOT, search


def input_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def selection(limit, manifest=None):
    if manifest:
        inputs = input_rows(manifest)
        inputs.sort(key=lambda r: hashlib.sha256(r["genome_id"].encode()).hexdigest())
        if limit:
            inputs = [r for source in sorted({x["source"] for x in inputs})
                      for r in [x for x in inputs if x["source"] == source][:limit]]
        if len({r["genome_id"] for r in inputs}) != len(inputs):
            raise ValueError("Selection contains repeated genome IDs")
        return inputs
    inputs = []
    for r in input_rows(ROOT / "inputs/selections/diversity_refseq_species.tsv"):
        acc = r["assembly"]
        base = r["ftp_path"]
        folder = base.rsplit("/", 1)[1]
        url = f"{base}/{folder}_protein.faa.gz"
        inputs.append({"source": "RefSeq", "genome_id": acc, "assembly": acc, "biosample": r["biosample"],
                       "url": url, "compression": "gzip"})
    for r in input_rows(ROOT / "inputs/selections/diversity_soil_species.tsv"):
        acc = r["Genome"]
        url = f"https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/soil/v1.0/species_catalogue/{acc[:-2]}/{acc}/genome/{acc}.faa"
        inputs.append({"source": "Soil", "genome_id": acc, "assembly": r["Genome_accession"],
                       "biosample": r["Sample_accession"], "url": url, "compression": "plain"})
    inputs.sort(key=lambda r: hashlib.sha256(r["genome_id"].encode()).hexdigest())
    if limit:
        inputs = [r for source in ("RefSeq", "Soil") for r in
                  [x for x in inputs if x["source"] == source][:limit]]
    return inputs


def fetch_batch(rows, raw, workers):
    paths = {row["genome_id"]: raw / (row["genome_id"] + row.get("file_suffix", ".faa.gz" if row["compression"] == "gzip" else ".faa")) for row in rows}
    missing = [row for row in rows if not paths[row["genome_id"]].exists()]
    durations = {}
    if missing:
        # One standard curl multi-transfer reuses connections across a batch.
        # Bound concurrency; validate individual HTTP responses before rename.
        command = ["curl", "--parallel", "--parallel-max", str(workers), "--fail",
                   "--location", "--retry", "3", "--retry-delay", "5", "--retry-all-errors",
                   "--silent", "--show-error", "--max-time", "1200",
                   "--write-out", "%{filename_effective}\t%{response_code}\t%{exitcode}\t%{time_total}\n"]
        temporary = {}
        for row in missing:
            path = paths[row["genome_id"]]
            temp = path.with_suffix(path.suffix + ".download")
            temporary[str(temp)] = (row["genome_id"], path)
            command += ["--url", row["url"], "--output", str(temp)]
        result = subprocess.run(command, capture_output=True, text=True)
        failed = []
        for line in result.stdout.splitlines():
            filename, code, exitcode, duration = line.split("\t")
            if code == "200" and exitcode == "0":
                identifier, path = temporary[filename]
                Path(filename).rename(path)
                durations[identifier] = float(duration)
            elif code != "000":
                failed.append((temporary[filename][0], code))
        if result.returncode or len(durations) != len(missing):
            raise RuntimeError(f"Incomplete proteome batch download, HTTP failures {failed[:5]}: {result.stderr[-1500:]}")
    records = []
    for row in rows:
        path = paths[row["genome_id"]]
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        records.append({**row, "local_path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                        "sha256": digest, "fetch_seconds": durations.get(row["genome_id"], 0),
                        "retrieved_utc": datetime.now(timezone.utc).isoformat()})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=12, help="Genomes per source; 0 means all")
    parser.add_argument("--selection", type=Path, help="Frozen annotated-proteome selection TSV")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--output-name", default="diversity_pilot")
    args = parser.parse_args()
    if (args.limit < 0 or args.batch_size < 1 or args.workers < 1 or args.cpu < 1
            or not 0 <= args.shard_index < args.shard_count):
        raise ValueError("Invalid screen size")
    selected = selection(args.limit, args.selection)[args.shard_index::args.shard_count]
    if not selected:
        raise ValueError("Empty species selection")
    if args.limit == 0 and args.output_name == "diversity_pilot":
        raise ValueError("Use a separate --output-name for the full screen to preserve the pilot")
    tables = ROOT / "results/tables"
    with (ROOT / "inputs/selections" / f"{args.output_name}_selection.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)
    raw = ROOT / "data/raw/diversity_proteomes"
    raw.mkdir(parents=True, exist_ok=True)
    interim = ROOT / "data/interim" / args.output_name
    interim.mkdir(parents=True, exist_ok=True)
    summaries, scores = [], []
    started = time.monotonic()
    # Append per-batch records: avoid repeatedly rewriting a genome-sized JSON
    # manifest. Pin URLs and hashes before deriving any scientific result.
    ledger = ROOT / "inputs/provenance" / f"{args.output_name}_downloads.jsonl"
    pinned = {r["genome_id"]: r for r in (json.loads(line) for line in ledger.read_text().splitlines())} if ledger.exists() else {}
    for offset in range(0, len(selected), args.batch_size):
        batch = selected[offset:offset+args.batch_size]
        result_dir = ROOT / "results/hmmer" / args.output_name / f"batch_{offset:05d}"
        saved_summary, saved_hits, saved_proteomes = (result_dir / name for name in
                                                    ("summary.json", "hits.tsv", "proteomes.tsv"))
        if saved_summary.exists() and saved_hits.exists() and saved_proteomes.exists():
            previous_rows = pd.read_csv(saved_proteomes, sep="\t", keep_default_na=False)
            if list(previous_rows.genome_id) != [r["genome_id"] for r in batch]:
                raise ValueError("Saved batch does not match fixed selection")
            summaries.extend(previous_rows.to_dict("records"))
            scores.append(pd.read_csv(saved_hits, sep="\t", keep_default_na=False))
            continue
        available = shutil.disk_usage(raw).free
        if available < 3 * 1024**3:
            raise OSError("Less than 3 GiB free; stop before downloading another batch")
        downloaded = fetch_batch(batch, raw, args.workers)
        with ledger.open("a") as handle:
            for r in downloaded:
                previous = pinned.get(r["genome_id"])
                if previous and any(previous[k] != r[k] for k in ["url", "bytes", "sha256"]):
                    raise ValueError(f"Pinned proteome changed: {r['genome_id']}")
                if not previous:
                    handle.write(json.dumps(r) + "\n")
        fasta = interim / "batch.faa"
        mapping = []
        with fasta.open("w") as output:
            for r in downloaded:
                opener = gzip.open if r["compression"] == "gzip" else open
                count = 0
                with opener(ROOT / r["local_path"], "rt") as source:
                    for count, (title, sequence) in enumerate(SimpleFastaParser(source), start=1):
                        identifier = f"{r['genome_id']}__protein{count:06d}"
                        output.write(f">{identifier}\n{sequence}\n")
                        mapping.append({"target_id": identifier, "genome_id": r["genome_id"], "source_protein_id": title.split()[0]})
                if not count:
                    raise ValueError(f"Empty annotated proteome: {r['genome_id']}")
                summaries.append({**r, "protein_count": count})
        screening_start = time.monotonic()
        result = search(fasta, result_dir, cpu=args.cpu)
        result = result.merge(pd.DataFrame(mapping), on="target_id", validate="many_to_one")
        scores.append(result)
        # Persist the exact nominated sequences with original IDs before the
        # batch FASTA is reused. This also makes addressed context checks cheap.
        nominated = set(result.target_id)
        with fasta.open() as source, (result_dir / "candidates.faa").open("w") as output:
            for title, sequence in SimpleFastaParser(source):
                if title.split()[0] in nominated:
                    output.write(f">{title}\n{sequence}\n")
        result.to_csv(saved_hits, sep="\t", index=False)
        pd.DataFrame(summaries[-len(batch):]).to_csv(saved_proteomes, sep="\t", index=False)
        batch_summary = {"batch_start": offset, "genomes": len(batch), "proteins": len(mapping),
                         "fasta_bytes": fasta.stat().st_size, "hmm_seconds": time.monotonic()-screening_start}
        saved_summary.write_text(json.dumps(batch_summary, indent=2) + "\n")
        print(f"{args.output_name}: genomes {offset+len(batch)}/{len(selected)}, candidates {len(result)}", flush=True)
    pd.concat(scores, ignore_index=True).to_csv(tables / f"{args.output_name}_profile_hits.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(tables / f"{args.output_name}_proteomes.tsv", sep="\t", index=False)
    summary = {"status": "resource_pilot_not_full_species_search" if args.limit else "completed_fixed_species_protein_screen_not_BGC_confirmation",
               "genomes": len(summaries), "proteins": sum(r["protein_count"] for r in summaries),
               "download_bytes": sum(r["bytes"] for r in summaries), "wall_seconds": time.monotonic()-started,
               "profile_hit_rows": sum(len(r) for r in scores)}
    summary["parent_max_rss_platform_units"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    summary["largest_child_max_rss_platform_units"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    summary["rss_units"] = "bytes_on_macOS_KiB_on_Linux"
    (tables / f"{args.output_name}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
