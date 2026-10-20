#!/usr/bin/env python3
"""Queue untested H2 host constraints in the frozen gCF inventory order."""

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCREEN = ROOT / "results/H3/global_screen"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    inventory = rows(SCREEN / "constraint_inventory.tsv")
    done = {(r["scope"], r["family"], r["host_split"]) for r in rows(SCREEN / "AU_summary.tsv")}
    models = {
        (r["scope"], r["dataset"]): r["best_model"]
        for r in rows(ROOT / "results/H2/comparisons/tree_quality.tsv")
        if r["kind"] == "component"
    }
    if len(inventory) != 174 or len(done) != 4:
        raise ValueError(f"Unexpected frozen screen size: {len(inventory)} constraints, {len(done)} AU tests")
    queue = []
    for rank, row in enumerate(inventory, start=1):
        key = row["scope"], row["family"], row["host_split"]
        if key in done:
            continue
        if key[:2] not in models:
            raise ValueError(f"No H2 model for {key}")
        constraint = SCREEN / "constraints" / "_".join(key)
        if not constraint.with_suffix(".nwk").is_file():
            raise FileNotFoundError(constraint)
        queue.append((len(queue) + 1, rank, *key, models[key[:2]]))
    if len(queue) != 170:
        raise ValueError(f"Expected 170 remaining tests, got {len(queue)}")
    with (SCREEN / "AU_queue.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("queue_index", "gCF_rank", "scope", "family", "host_split", "model"))
        writer.writerows(queue)
    print(f"Queued {len(queue)} untested constraints in gCF order")


if __name__ == "__main__":
    main()
