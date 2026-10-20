#!/usr/bin/env python3
"""Inventory distinct host constraints behind all supported H2 split conflicts."""

from collections import defaultdict
import csv
import json

from Bio import Phylo, SeqIO

from prepare_h2 import OUT, ROOT, read
from h2_review_relations import edges


def main():
    conflicts = read(OUT / "comparisons/supported_split_conflicts.tsv")
    queries = json.loads((ROOT / "config/H2/H2_review_relation_queries.json").read_text())["candidates"]
    query_ids = defaultdict(list)
    for query in queries:
        query_ids[(query["scope"], query["family"], query["host_split"])].append(query["candidate_id"])
    grouped = defaultdict(list)
    for row in conflicts:
        grouped[(row["scope"], row["family"], row["host_split"])].append(row)

    output = []
    destination = OUT.parent / "H3/global_screen"
    constraints_dir = destination / "constraints"
    constraints_dir.mkdir(parents=True, exist_ok=True)
    for (scope_name, family, split_id), rows in sorted(grouped.items()):
        row = rows[0]
        scope = OUT / "actino" if scope_name == "Actinobacteria138" else OUT
        host_name = "uncertain_Nature" if family == "uncertain_Nature_system_C" else family.split("_")[0]
        host_dir = scope / "host_trees" / host_name
        tree = Phylo.read(host_dir / "concordance.cf.branch", "newick")
        with (host_dir / "concordance.cf.stat").open() as handle:
            stats = {r["ID"]: r for r in csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")}

        shared = set(row["shared_hosts"].split(";"))
        host_side = set(row["host_side"].split(";"))
        matches = [stats[str(int(node.confidence))] for node in edges(tree, host_side, shared)]
        if not matches:
            raise ValueError(f"No gCF branch for {scope_name}/{family}/{split_id}")
        if any(r["shared_hosts"] != row["shared_hosts"] or r["host_side"] != row["host_side"] for r in rows):
            raise ValueError(f"Inconsistent split identity for {scope_name}/{family}/{split_id}")

        alignment = scope / "component_trees" / family / "trimmed_alignment.faa"
        tips = [record.id for record in SeqIO.parse(alignment, "fasta")]
        by_host = defaultdict(list)
        for tip in tips:
            by_host[tip.split("_", 1)[0]].append(tip)
        if any(len(by_host[host]) != 1 for host in shared):
            raise ValueError(f"Constraint has missing/multiple copies for {scope_name}/{family}/{split_id}")
        side_tips = sorted(by_host[host][0] for host in host_side)
        other_tips = sorted(by_host[host][0] for host in shared - host_side)
        if min(len(side_tips), len(other_tips)) < 2:
            raise ValueError(f"Constraint has trivial side for {scope_name}/{family}/{split_id}")
        constraint = constraints_dir / f"{scope_name}_{family}_{split_id}.nwk"
        constraint.write_text("((" + ",".join(side_tips) + ")," + ",".join(other_tips) + ");\n")
        selected = query_ids[(scope_name, family, split_id)]

        output.append({
            "scope": scope_name,
            "family": family,
            "host_split": split_id,
            "shared_single_copy_hosts": len(shared),
            "host_side": row["host_side"],
            "host_SH_aLRT": row["host_SH_aLRT"],
            "host_UFBoot": row["host_UFBoot"],
            "host_gCF_percent": ";".join(r["gCF"] for r in matches),
            "host_gCF_N_over_gN": ";".join(f'{r["gCF_N"]}/{r["gN"]}' for r in matches),
            "host_gCF_branch_ids": ";".join(r["ID"] for r in matches),
            "conflicting_gene_splits": len({r["gene_split"] for r in rows}),
            "conflicting_pair_rows": len(rows),
            "gene_split_ids": ";".join(sorted({r["gene_split"] for r in rows})),
            "selected_H2_relation_ids": ";".join(selected),
            "alignment": str(alignment.relative_to(ROOT)),
            "interpretation": "host_marker_concordance_only; not_HGT_probability_or_event",
        })

    output.sort(key=lambda r: (-min(float(x) for x in r["host_gCF_percent"].split(";")),
                               -min(int(x.split("/")[1]) for x in r["host_gCF_N_over_gN"].split(";")),
                               r["scope"], r["family"], r["host_split"]))
    with (destination / "constraint_inventory.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"{len(conflicts)} supported conflict pairs -> {len(output)} distinct family/panel/host constraints; "
          "gCF mapped for all; AU reported separately")


if __name__ == "__main__":
    main()
