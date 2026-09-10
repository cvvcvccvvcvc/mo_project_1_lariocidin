#!/usr/bin/env python3
"""Shared, bounded HMMER search and explicit domain-coverage accounting."""
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROFILES = ("canonical", "noncanonical")
SCORE_COLUMNS = ["target_id", "profile", "sequence_score", "evalue", "bias",
                 "target_length", "model_length", "model_coverage", "target_coverage",
                 "best_domain_score", "best_domain_ievalue", "domain_count", "priority"]


def union_length(intervals):
    total, last = 0, 0
    for start, end in sorted(intervals):
        total += max(0, end - max(start - 1, last))
        last = max(last, end)
    return total


def read_scores(tbl, domtbl, profile):
    domains = defaultdict(list)
    for line in domtbl.read_text().splitlines():
        if line and not line.startswith("#"):
            fields = line.split()
            domains[fields[0]].append(fields)
    rows = []
    for line in tbl.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        ds = domains[f[0]]
        model_length = int(ds[0][5]) if ds else 0
        target_length = int(ds[0][2]) if ds else 0
        mcov = union_length([(int(d[15]), int(d[16])) for d in ds]) / model_length if ds else 0
        tcov = union_length([(int(d[17]), int(d[18])) for d in ds]) / target_length if ds else 0
        rows.append(dict(zip(SCORE_COLUMNS, [
            f[0], profile, float(f[5]), float(f[4]), float(f[6]), target_length,
            model_length, mcov, tcov, max((float(d[13]) for d in ds), default=None),
            min((float(d[12]) for d in ds), default=None), len(ds),
            "high" if float(f[5]) >= 250 and mcov >= .70 else "review",
        ])))
    return rows


def search(fasta, output_dir, exhaustive=False):
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, commands = [], []
    for profile in PROFILES:
        prefix = output_dir / profile
        tbl, dom = prefix.with_suffix(".tbl"), prefix.with_suffix(".domtbl")
        command = [str(ROOT / ".env/bin/hmmsearch"), "--cpu", "4", "--noali",
                   "--seed", "42", "-T", "80", "--domT", "20",
                   "--tblout", str(tbl), "--domtblout", str(dom)]
        if exhaustive:
            command.append("--max")
        command += [str(ROOT / f"data/processed/literature_profiles/final.{profile}.hmm"), str(fasta)]
        with prefix.with_suffix(".log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
        commands.append(command)
        rows.extend(read_scores(tbl, dom, profile))
    (output_dir / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    return pd.DataFrame(rows, columns=SCORE_COLUMNS)
