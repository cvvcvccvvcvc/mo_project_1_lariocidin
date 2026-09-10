#!/usr/bin/env python3
"""Locate the focal unrooted split and summarize support across tree variants."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import AlignIO, Phylo


FOCAL = {
    "burkholderia_imcc1007",
    "ralstonia_pickettii_dtp0602",
    "paraburkholderia_mimosarum_stm3621",
    "photorhabdus_asymbiotica_dsm15149",
}


def parse_support(clade) -> tuple[float | None, float | None]:
    value = clade.name
    if value is None and clade.confidence is not None:
        value = str(clade.confidence)
    if value is None:
        return None, None
    match = re.fullmatch(r"([0-9.]+)(?:/([0-9.]+))?", str(value))
    if not match:
        return None, None
    return float(match.group(1)), float(match.group(2)) if match.group(2) else None


def focal_support(tree_path: Path) -> tuple[bool, float | None, float | None]:
    tree = Phylo.read(tree_path, "newick")
    all_tips = {tip.name for tip in tree.get_terminals()}
    if not FOCAL <= all_tips:
        raise ValueError(f"Focal taxa absent from {tree_path}")
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        side = {tip.name for tip in clade.get_terminals()}
        if side == FOCAL or all_tips - side == FOCAL:
            sh_alrt, ufboot = parse_support(clade)
            return True, sh_alrt, ufboot
    return False, None, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", action="append", nargs=3, metavar=("NAME", "TREE", "ALIGNMENT"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for name, tree_name, alignment_name in args.variant:
        tree_path, alignment_path = Path(tree_name), Path(alignment_name)
        alignment = AlignIO.read(alignment_path, "fasta")
        found, sh_alrt, ufboot = focal_support(tree_path)
        rows.append(
            {
                "variant": name,
                "sequence_count": len(alignment),
                "alignment_columns": alignment.get_alignment_length(),
                "focal_split_present": found,
                "sh_alrt": "" if sh_alrt is None else sh_alrt,
                "ufboot": "" if ufboot is None else ufboot,
                "meets_stronger_support_guide": bool(
                    found
                    and sh_alrt is not None
                    and ufboot is not None
                    and sh_alrt >= 80
                    and ufboot >= 95
                ),
            }
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
