#!/usr/bin/env python3
"""Collect the frozen H3 AU screen and verify each constrained host split."""

import csv
import hashlib
from decimal import Decimal
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "results/H3/global_screen"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def au_scores(path):
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "USER TREES")
    scores = {}
    for line in lines[start + 1:]:
        parts = line.split()
        if parts and parts[0] in ("1", "2") and len(parts) >= 8:
            scores[int(parts[0])] = (
                float(parts[1]), float(parts[2]), float(parts[-2]),
                10 ** Decimal(parts[2]).as_tuple().exponent,
            )
        if len(scores) == 2:
            break
    if len(scores) != 2 or not any("All tests performed 10000 resamplings" in line for line in lines):
        raise ValueError(f"Incomplete two-tree AU report: {path}")
    return scores


def tip_names(tree):
    names = [tip.name for tip in tree.get_terminals()]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate leaf name")
    return set(names)


def has_projected_split(tree, side, shared):
    for node in tree.get_nonterminals():
        projected = {tip.name for tip in node.get_terminals()} & shared
        if projected == side or shared - projected == side:
            return True
    return False


def main():
    inventory = read_tsv(SCREEN / "constraint_inventory.tsv")
    previous = {(r["scope"], r["family"], r["host_split"]): r for r in read_tsv(SCREEN / "AU_summary.tsv")}
    quality = {
        (r["scope"], r["dataset"]): r
        for r in read_tsv(ROOT / "results/H2/comparisons/tree_quality.tsv") if r["kind"] == "component"
    }
    output = []
    for rank, r in enumerate(inventory, 1):
        scope, family, split = r["scope"], r["family"], r["host_split"]
        key = scope, family, split
        case = "_".join(key)
        constraint_path = SCREEN / "constraints" / f"{case}.nwk"
        row = {
            "gCF_rank": rank, "case": case, "scope": scope, "family": family, "host_split": split,
            "host_side": r["host_side"], "host_gCF_percent": r["host_gCF_percent"],
            "host_gCF_N_over_gN": r["host_gCF_N_over_gN"], "host_SH_aLRT": r["host_SH_aLRT"],
            "host_UFBoot": r["host_UFBoot"], "conflicting_H2_pairs": r["conflicting_pair_rows"],
            "component_tips": quality[(scope, family)]["tips"],
            "alignment_columns": quality[(scope, family)]["columns"],
            "composition_chi2_failures": quality[(scope, family)]["composition_chi2_failures"],
            "model": quality[(scope, family)]["best_model"],
            "free_logL": "", "constrained_logL": "", "delta_logL": "",
            "constrained_p_AU": "", "AU_excludes_host_at_0.05": "",
            "bonferroni_p_174": "", "bonferroni_pass_familywise_0.05": "",
            "constrained_split_present": "",
            "same_tip_set": "", "source": "", "status": "pending",
            "constraint_sha256": hashlib.sha256(constraint_path.read_bytes()).hexdigest(),
            "constrained_tree_sha256": "",
            "p_AU_rank_of_174": "",
        }
        if key in previous:
            old = previous[key]
            row.update({k: old[k] for k in ("free_logL", "constrained_logL", "delta_logL", "constrained_p_AU")})
            row["AU_excludes_host_at_0.05"] = str(float(old["constrained_p_AU"]) < 0.05).lower()
            row["constrained_split_present"] = "true"
            row["same_tip_set"] = "true"
            row["source"] = old["run"]
            row["status"] = "complete_prior"
            prior_tree = (ROOT / "results/H3/topology/BURK_C_LG_G4/constrained.treefile"
                          if old["case"] == "H2-CAN-BURK" else SCREEN / "AU" / case / "constrained.treefile")
            row["constrained_tree_sha256"] = hashlib.sha256(prior_tree.read_bytes()).hexdigest()
        else:
            folder = SCREEN / "AU" / case
            report, tree_path = folder / "AU.iqtree", folder / "constrained.treefile"
            if report.is_file() and tree_path.is_file():
                scores = au_scores(report)
                panel = "actino/" if scope == "Actinobacteria138" else ""
                free_path = ROOT / "results/H2" / panel / "component_trees" / family / "ML.treefile"
                free, constrained = Phylo.read(free_path, "newick"), Phylo.read(tree_path, "newick")
                side_tree = Phylo.read(constraint_path, "newick")
                side = tip_names(side_tree.clade.clades[0])
                shared = tip_names(side_tree)
                same = tip_names(free) == tip_names(constrained)
                present = has_projected_split(constrained, side, shared)
                if not same or not present or not shared <= tip_names(free):
                    raise ValueError(f"Tree/constraint mismatch: {case}")
                best_logl = max(scores[1][0], scores[2][0])
                for number in (1, 2):
                    logl, delta, _, resolution = scores[number]
                    if abs(best_logl - logl - delta) > 0.5 * float(resolution) + 0.0001:
                        raise ValueError(f"Inconsistent AU likelihoods: {case}")
                row.update({
                    "free_logL": scores[1][0], "constrained_logL": scores[2][0],
                    "delta_logL": round(scores[1][0] - scores[2][0], 5),
                    "constrained_p_AU": scores[2][2],
                    "AU_excludes_host_at_0.05": str(scores[2][2] < 0.05).lower(),
                    "constrained_split_present": "true", "same_tip_set": "true",
                    "source": "H3_AU_array_2026-10-20", "status": "complete",
                    "constrained_tree_sha256": hashlib.sha256(tree_path.read_bytes()).hexdigest(),
                })
        if row["constrained_p_AU"] != "":
            adjusted = min(1.0, 174 * float(row["constrained_p_AU"]))
            row["bonferroni_p_174"] = f"{adjusted:.6g}"
            row["bonferroni_pass_familywise_0.05"] = str(adjusted < 0.05).lower()
        output.append(row)
    complete = sum(r["status"].startswith("complete") for r in output)
    if complete == len(output):
        ranked = sorted(output, key=lambda r: (float(r["constrained_p_AU"]), r["gCF_rank"]))
        for rank, row in enumerate(ranked, 1):
            row["p_AU_rank_of_174"] = rank
    target = SCREEN / "AU_full_screen.tsv"
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"AU complete {complete}/{len(output)}; pending {len(output) - complete}")


if __name__ == "__main__":
    main()
