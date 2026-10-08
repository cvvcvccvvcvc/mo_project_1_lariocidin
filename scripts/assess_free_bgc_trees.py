#!/usr/bin/env python3
"""Assess the prespecified order split and extract selected IQ-TREE models."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import Phylo


ROLES = ("B1", "B2", "C", "D")


def validate_taxa(names, expected: set[str], label: str) -> None:
    names = list(names)
    if len(names) != len(set(names)) or set(names) != expected:
        raise ValueError(f"Taxon mismatch or duplicate IDs in {label}")


def has_split(tree_path: Path, group: set[str], expected: set[str] | None = None) -> bool:
    tree = Phylo.read(tree_path, "newick")
    if expected is not None:
        validate_taxa((leaf.name for leaf in tree.get_terminals()), expected, str(tree_path))
    leaves = {leaf.name for leaf in tree.get_terminals()}
    if not group < leaves:
        raise ValueError(f"Invalid or trivial group for {tree_path}")
    complement = leaves - group
    return any(
        {leaf.name for leaf in clade.get_terminals()} in (group, complement)
        for clade in tree.find_clades(order="preorder")
        if clade is not tree.root
    )


def model_from_report(path: Path) -> str:
    match = re.search(r"^Best-fit model according to BIC:\s+(\S+)", path.read_text(), re.MULTILINE)
    if not match:
        raise ValueError(f"Best-fit model missing from {path}")
    return match.group(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--constraint-out", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    selected = selected.loc[selected["phylum"] == "Pseudomonadota"]
    groups = {
        name: set(subset["sample_id"])
        for name, subset in selected.groupby("group")
    }
    if set(groups) != {"Enterobacterales", "Burkholderiales"} or sorted(map(len, groups.values())) != [6, 8]:
        raise ValueError("Expected frozen 6/8 order groups")
    enterobacterales = sorted(groups["Enterobacterales"])
    burkholderiales = sorted(groups["Burkholderiales"])
    args.constraint_out.parent.mkdir(parents=True, exist_ok=True)
    args.constraint_out.write_text(
        f"(({','.join(enterobacterales)}),({','.join(burkholderiales)}));\n"
    )

    rows = []
    for role in ROLES:
        prefix = args.tree_dir / f"{role}.free"
        rows.append({
            "role": role,
            "order_split_present": has_split(
                Path(str(prefix) + ".treefile"), groups["Enterobacterales"], set(selected.sample_id)
            ),
            "model": model_from_report(Path(str(prefix) + ".iqtree")),
        })
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
