#!/usr/bin/env python3
"""Cluster-aware cross-validation and final HMMs for published cyclases.

The cross-validation unit is a 70%-identity connected component with at least
80% reciprocal query coverage.  Close homologues therefore stay in the same fold.
The Barrett-review set is used only as a lasso-cyclase comparator, not as proof
that every possible non-family protein is represented.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


BLAST_FIELDS = [
    "query_id",
    "subject_id",
    "pident",
    "alignment_length",
    "query_length",
    "subject_length",
    "evalue",
    "bitscore",
    "query_coverage_all_hsps",
    "query_coverage_hsp",
]
POSITIVE_GROUPS = ("canonical", "noncanonical")
COMPARATOR_GROUP = "Barrett et al review"


class UnionFind:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: str) -> str:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def run(command: list[str], stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(command, check=True, capture_output=True, text=True)
    else:
        with stdout.open("w", encoding="utf-8") as handle:
            subprocess.run(command, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)


def connected_components(
    sequence_ids: list[str],
    pair_hits: dict[tuple[str, str], tuple[float, float]],
    identity: float,
    reciprocal_coverage: float,
) -> list[list[str]]:
    uf = UnionFind(sequence_ids)
    sequence_set = set(sequence_ids)
    for (left, right), (left_identity, left_coverage) in pair_hits.items():
        if left >= right or left not in sequence_set or right not in sequence_set:
            continue
        reverse = pair_hits.get((right, left))
        if (
            reverse
            and left_identity >= identity
            and reverse[0] >= identity
            and left_coverage >= reciprocal_coverage
            and reverse[1] >= reciprocal_coverage
        ):
            uf.union(left, right)
    components: dict[str, list[str]] = defaultdict(list)
    for sequence_id in sequence_ids:
        components[uf.find(sequence_id)].append(sequence_id)
    return sorted(components.values(), key=lambda values: (-len(values), values[0]))


def balanced_cluster_folds(clusters: list[list[str]], fold_count: int) -> dict[str, int]:
    fold_sizes = [0] * fold_count
    assignment: dict[str, int] = {}
    for cluster in clusters:
        fold = min(range(fold_count), key=lambda index: (fold_sizes[index], index))
        for sequence_id in cluster:
            assignment[sequence_id] = fold
        fold_sizes[fold] += len(cluster)
    return assignment


def write_fasta(records: dict[str, object], ids: list[str], path: Path) -> None:
    SeqIO.write([records[sequence_id] for sequence_id in ids], path, "fasta")


def build_profile(
    group: str,
    sequence_ids: list[str],
    records: dict[str, object],
    prefix: Path,
    mafft: Path,
    hmmbuild: Path,
) -> Path:
    safe_group = group.replace(" ", "_")
    fasta = prefix.with_name(f"{prefix.name}.{safe_group}.faa")
    alignment = prefix.with_name(f"{prefix.name}.{safe_group}.afa")
    hmm = prefix.with_name(f"{prefix.name}.{safe_group}.hmm")
    write_fasta(records, sequence_ids, fasta)
    run([str(mafft), "--auto", "--quiet", str(fasta)], stdout=alignment)
    run(
        [
            str(hmmbuild),
            "--amino",
            "--informat",
            "afa",
            "-n",
            f"lariocidin_{safe_group}",
            str(hmm),
            str(alignment),
        ]
    )
    return hmm


def search_profile(hmmsearch: Path, hmm: Path, database: Path, tblout: Path) -> dict[str, float]:
    run(
        [
            str(hmmsearch),
            "--max",
            "--noali",
            "-T",
            "-999",
            "--tblout",
            str(tblout),
            str(hmm),
            str(database),
        ]
    )
    scores: dict[str, float] = {}
    with tblout.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split(maxsplit=18)
            scores[fields[0]] = float(fields[5])
    return scores


def metrics(scores: pd.DataFrame, threshold: float) -> dict[str, float | int]:
    positive = scores["label"] == "positive"
    predicted = scores["best_bitscore"] >= threshold
    tp = int((positive & predicted).sum())
    fn = int((positive & ~predicted).sum())
    fp = int((~positive & predicted).sum())
    tn = int((~positive & ~predicted).sum())
    return {
        "threshold": float(threshold),
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "sensitivity": tp / (tp + fn),
        "comparator_specificity": tn / (tn + fp),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cyclase-fasta", type=Path, required=True)
    parser.add_argument("--cyclase-table", type=Path, required=True)
    parser.add_argument("--all-v-all", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--hmmbuild", type=Path, required=True)
    parser.add_argument("--hmmsearch", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--scores-out", type=Path, required=True)
    parser.add_argument("--folds-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--cluster-identity", type=float, default=70.0)
    parser.add_argument("--cluster-coverage", type=float, default=80.0)
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(args.cyclase_table, sep="\t", keep_default_na=False)
    group_by_id = table.set_index("cyclase_id")["source_groups"].to_dict()
    records = {record.id: record for record in SeqIO.parse(args.cyclase_fasta, "fasta")}
    if set(records) != set(group_by_id):
        raise ValueError("Cyclase FASTA and metadata IDs differ")

    blast = pd.read_csv(args.all_v_all, sep="\t", names=BLAST_FIELDS)
    pair_best = (
        blast.sort_values(
            ["query_id", "subject_id", "bitscore"], ascending=[True, True, False]
        )
        .drop_duplicates(["query_id", "subject_id"], keep="first")
    )
    pair_hits = {
        (row.query_id, row.subject_id): (row.pident, row.query_coverage_all_hsps)
        for row in pair_best.itertuples()
    }

    cluster_count_sensitivity: dict[str, dict[str, int]] = {}
    for identity in (50.0, 70.0, 80.0):
        cluster_count_sensitivity[str(int(identity))] = {}
        for group in POSITIVE_GROUPS:
            ids = table.loc[table["source_groups"] == group, "cyclase_id"].tolist()
            clusters = connected_components(ids, pair_hits, identity, args.cluster_coverage)
            cluster_count_sensitivity[str(int(identity))][group] = len(clusters)

    fold_rows: list[dict[str, object]] = []
    fold_by_group: dict[str, dict[str, int]] = {}
    for group in POSITIVE_GROUPS:
        ids = table.loc[table["source_groups"] == group, "cyclase_id"].tolist()
        clusters = connected_components(
            ids, pair_hits, args.cluster_identity, args.cluster_coverage
        )
        assignment = balanced_cluster_folds(clusters, args.folds)
        fold_by_group[group] = assignment
        cluster_id_by_sequence = {
            sequence_id: cluster_index
            for cluster_index, cluster in enumerate(clusters, start=1)
            for sequence_id in cluster
        }
        for sequence_id in ids:
            fold_rows.append(
                {
                    "cyclase_id": sequence_id,
                    "source_group": group,
                    "cluster_id": f"{group}_{cluster_id_by_sequence[sequence_id]:03d}",
                    "fold": assignment[sequence_id] + 1,
                }
            )
    fold_table = pd.DataFrame(fold_rows).sort_values(["source_group", "fold", "cluster_id"])
    args.folds_out.parent.mkdir(parents=True, exist_ok=True)
    fold_table.to_csv(args.folds_out, sep="\t", index=False)

    all_scores: dict[tuple[int, str], dict[str, float]] = {}
    for fold in range(args.folds):
        for group in POSITIVE_GROUPS:
            ids = [
                sequence_id
                for sequence_id, source_group in group_by_id.items()
                if source_group == group and fold_by_group[group][sequence_id] != fold
            ]
            prefix = args.work_dir / f"cv_fold_{fold + 1}"
            hmm = build_profile(group, ids, records, prefix, args.mafft, args.hmmbuild)
            tblout = args.work_dir / f"cv_fold_{fold + 1}.{group}.tblout"
            all_scores[(fold, group)] = search_profile(
                args.hmmsearch, hmm, args.cyclase_fasta, tblout
            )

    score_rows: list[dict[str, object]] = []
    for sequence_id, source_group in group_by_id.items():
        if source_group in POSITIVE_GROUPS:
            fold = fold_by_group[source_group][sequence_id]
            canonical_score = all_scores[(fold, "canonical")][sequence_id]
            noncanonical_score = all_scores[(fold, "noncanonical")][sequence_id]
            score_rows.append(
                {
                    "cyclase_id": sequence_id,
                    "source_group": source_group,
                    "label": "positive",
                    "evaluation_fold": fold + 1,
                    "canonical_bitscore": canonical_score,
                    "noncanonical_bitscore": noncanonical_score,
                    "best_bitscore": max(canonical_score, noncanonical_score),
                }
            )
        else:
            fold_scores = []
            for fold in range(args.folds):
                canonical_score = all_scores[(fold, "canonical")][sequence_id]
                noncanonical_score = all_scores[(fold, "noncanonical")][sequence_id]
                fold_scores.append((max(canonical_score, noncanonical_score), fold, canonical_score, noncanonical_score))
            best, fold, canonical_score, noncanonical_score = max(fold_scores)
            score_rows.append(
                {
                    "cyclase_id": sequence_id,
                    "source_group": source_group,
                    "label": "comparator",
                    "evaluation_fold": f"worst:{fold + 1}",
                    "canonical_bitscore": canonical_score,
                    "noncanonical_bitscore": noncanonical_score,
                    "best_bitscore": best,
                }
            )
    scores = pd.DataFrame(score_rows).sort_values(["label", "source_group", "cyclase_id"])
    args.scores_out.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.scores_out, sep="\t", index=False)

    positive_scores = scores.loc[scores["label"] == "positive", "best_bitscore"]
    comparator_scores = scores.loc[scores["label"] == "comparator", "best_bitscore"]
    zero_comparator_fp_threshold = float(comparator_scores.max() + 0.01)
    positive_recall_95_threshold = float(np.quantile(positive_scores, 0.05))
    tested_thresholds = sorted(
        {
            round(zero_comparator_fp_threshold, 2),
            round(positive_recall_95_threshold, 2),
            round(zero_comparator_fp_threshold * 0.8, 2),
            round(zero_comparator_fp_threshold * 1.2, 2),
        }
    )

    final_profiles = {}
    for group in POSITIVE_GROUPS:
        ids = table.loc[table["source_groups"] == group, "cyclase_id"].tolist()
        hmm = build_profile(
            group,
            ids,
            records,
            args.work_dir / "final",
            args.mafft,
            args.hmmbuild,
        )
        final_profiles[group] = str(hmm)

    fold_sizes = {
        group: {
            str(fold + 1): int(sum(value == fold for value in assignment.values()))
            for fold in range(args.folds)
        }
        for group, assignment in fold_by_group.items()
    }
    summary = {
        "status": "preliminary_internal_calibration_not_external_validation",
        "positive_definition": "canonical or noncanonical cyclase explicitly labelled in the published Merrild archive",
        "comparator_definition": "other lasso cyclases from the Barrett-review archive subset",
        "cluster_rule": {
            "identity_percent": args.cluster_identity,
            "reciprocal_query_coverage_percent": args.cluster_coverage,
            "linkage": "connected_components_single_linkage",
        },
        "cluster_count_sensitivity": cluster_count_sensitivity,
        "fold_count": args.folds,
        "fold_sizes": fold_sizes,
        "positive_score_quantiles": {
            str(q): float(positive_scores.quantile(q)) for q in [0, 0.01, 0.05, 0.5, 0.95, 1]
        },
        "comparator_score_quantiles": {
            str(q): float(comparator_scores.quantile(q)) for q in [0, 0.5, 0.95, 0.99, 1]
        },
        "empirical_zero_comparator_fp_threshold": zero_comparator_fp_threshold,
        "empirical_95pct_positive_recall_threshold": positive_recall_95_threshold,
        "threshold_metrics": [metrics(scores, value) for value in tested_thresholds],
        "final_profiles": final_profiles,
        "limitations": [
            "The source groups were sequence-informed, so this is not independent validation.",
            "The comparator set is small and does not represent all non-family proteins in MAGs.",
            "A profile hit remains a candidate until precursor and neighbourhood context are checked.",
        ],
    }
    args.summary_out.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
