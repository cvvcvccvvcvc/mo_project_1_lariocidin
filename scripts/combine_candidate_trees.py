#!/usr/bin/env python3
"""Write free and constrained Newick topologies for an IQ-TREE topology test."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--free", type=Path, required=True)
    parser.add_argument("--constrained", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    trees = [args.free.read_text().strip(), args.constrained.read_text().strip()]
    if not all(tree.endswith(";") for tree in trees):
        raise ValueError("Expected two complete Newick trees")
    args.output.write_text("\n".join(trees) + "\n")


if __name__ == "__main__":
    main()
