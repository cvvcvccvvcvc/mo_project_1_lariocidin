#!/usr/bin/env python3
"""One standard alignment/ML task with explicit input and site diagnostics."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess

from Bio import SeqIO

from prepare_h2 import ROOT, OUT, read

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def alignment_stats(path):
    seqs = [str(r.seq).upper() for r in SeqIO.parse(path, "fasta")]
    if not seqs or len({len(s) for s in seqs}) != 1:
        raise ValueError(f"Empty or unequal alignment: {path}")
    variable = informative = 0
    for column in zip(*seqs):
        counts = Counter(x for x in column if x in AMINO_ACIDS)
        variable += len(counts) > 1
        informative += sum(n >= 2 for n in counts.values()) >= 2
    return {"tips": len(seqs), "unique_sequences": len(set(seqs)), "columns": len(seqs[0]),
            "variable_sites": variable, "parsimony_informative_sites": informative,
            "missing_fraction": sum(x not in AMINO_ACIDS for s in seqs for x in s) / (len(seqs) * len(seqs[0]))}


def execute(arguments, output):
    with output.open("w") as handle:
        subprocess.run([str(x) for x in arguments], stdout=handle, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["marker", "component", "host"], required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--actino", action="store_true")
    parser.add_argument("--fixed-host-models", action="store_true")
    parser.add_argument("--diagnostic-fast", action="store_true")
    args = parser.parse_args()
    destination = OUT / "actino" if args.actino else OUT
    bin_root = Path(os.environ["HGT_BIN"])
    small_bin = Path(os.environ["HGT_SMALL_BIN"])
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
    if args.mode == "marker":
        row = read(destination / "marker_models.tsv")[args.index]
        name = row["marker_id"]
        input_file = destination / "markers" / (name + ".faa")
        directory = destination / "marker_trees" / name
    else:
        rows = [r for r in read(destination / "tree_jobs.tsv") if r["kind"] == args.mode]
        row = rows[args.index]
        name = row["name"]
        input_file = ROOT / row["input_path"]
        directory = destination / ("host_trees" if args.mode == "host" else "component_trees") / name
    directory.mkdir(parents=True, exist_ok=True)
    if args.mode != "host":
        with (directory / "raw_alignment.faa").open("w") as alignment, (directory / "mafft.log").open("w") as error:
            subprocess.run([str(bin_root / "mafft"), "--auto", "--thread", str(threads), str(input_file)], stdout=alignment, stderr=error, check=True)
        if args.mode == "component" and row.get("filter") == "untrimmed":
            alignment_file = directory / "raw_alignment.faa"
        else:
            alignment_file = directory / "trimmed_alignment.faa"
            execute([small_bin / "trimal", "-in", directory / "raw_alignment.faa", "-out", alignment_file, "-gappyout"], directory / "trimal.log")
        raw = alignment_stats(directory / "raw_alignment.faa")
        filtered = alignment_stats(alignment_file)
        stats = {"name": name, "kind": args.mode, "raw": raw, "filtered": filtered,
                 "retained_fraction": filtered["columns"] / raw["columns"]}
    else:
        alignment_file = input_file
        stats = {"name": name, "kind": args.mode, "filtered": alignment_stats(input_file)}
    if stats["filtered"]["unique_sequences"] < 4 or stats["filtered"]["parsimony_informative_sites"] < 10:
        stats["tree_status"] = "insufficient_distinct_sequences_or_informative_sites; no_supported_tree_claim"
        (directory / "task_summary.json").write_text(json.dumps(stats, indent=2) + "\n")
        print(json.dumps(stats))
        return
    prefix = "ML_fastdiagnostic" if args.mode == "marker" else "ML"
    command = [bin_root / "iqtree", "-s", alignment_file, "-nt", str(threads), "-seed", str(131013 + args.index),
               "-alrt", "1000", "--keep-ident", "--prefix", directory / prefix]
    if not (args.mode == "host" and args.fixed_host_models):
        command += ["-m", "MFP+MERGE" if args.mode == "host" else "MFP", "-mset", "LG,WAG,JTT", "-mrate", "G"]
    if args.mode == "marker" or args.diagnostic_fast:
        command += ["-fast"]
    else:
        command += ["-B", "1000", "--bnni"]
    if args.mode == "host":
        command += ["-p", ROOT / row["partition_path"]]
        if not args.fixed_host_models:
            command += ["-rcluster", "10"]
    execute(command, directory / "iqtree.log")
    if not (directory / (prefix + ".treefile")).exists():
        raise ValueError("IQ-TREE completed without ML tree")
    stats["tree_status"] = "approximate_ML_diagnostic; SH_aLRT_support_not_UFBoot" if args.mode == "marker" or args.diagnostic_fast else "ML_complete; support_is_not_HGT_probability"
    (directory / "task_summary.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
