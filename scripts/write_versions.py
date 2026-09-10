#!/usr/bin/env python3
"""Record exact versions of the top-level tools used in this analysis."""

from __future__ import annotations

import argparse
import importlib.metadata
import subprocess
import sys
from pathlib import Path


def command_version(path: Path, args: list[str]) -> str:
    result = subprocess.run([str(path), *args], text=True, capture_output=True, check=True)
    return (result.stdout or result.stderr).splitlines()[0].strip()


def hmmer_version(path: Path) -> str:
    result = subprocess.run([str(path), "-h"], text=True, capture_output=True, check=True)
    return next(line.removeprefix("# ") for line in result.stdout.splitlines() if line.startswith("# HMMER "))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    bin_dir = Path(sys.executable).parent
    versions = [
        ("python", sys.version.split()[0]),
        ("biopython", importlib.metadata.version("biopython")),
        ("pandas", importlib.metadata.version("pandas")),
        ("matplotlib", importlib.metadata.version("matplotlib")),
        ("mafft", command_version(bin_dir / "mafft", ["--version"])),
        ("iqtree", command_version(bin_dir / "iqtree", ["-version"])),
        ("blastp", command_version(bin_dir / "blastp", ["-version"])),
        ("hmmer", hmmer_version(bin_dir / "hmmsearch")),
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("tool\tversion\n" + "".join(f"{tool}\t{version}\n" for tool, version in versions))
    print(f"recorded {len(versions)} tool versions")


if __name__ == "__main__":
    main()
