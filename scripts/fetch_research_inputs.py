#!/usr/bin/env python3
"""Restore archived public metadata; fetch the pinned literature archive.

The ZIP contains original API/XML/TSV snapshots, checked against the existing
source manifests. Replaying an old selection must not query today's live API.
"""
import csv
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def restore_snapshots(root=ROOT):
    expected = {}
    for name in ("s3_download_manifest.json", "s5_download_manifest.json"):
        for row in json.loads((root / "config" / name).read_text()).values():
            expected[row["local_path"]] = row["sha256"]
    for name, row in json.loads((root / "config/maize_biosample_downloads.json").read_text()).items():
        expected[f"data/raw/s3/{name}"] = row["sha256"]
    with (root / "config/research_sources_2026-09-10.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected[row["local_path"]] = row["sha256"]
    with zipfile.ZipFile(root / "config/source_snapshots.zip") as archive:
        for name in archive.namelist():
            path = root / name
            if name not in expected or not path.resolve().is_relative_to((root / "data/raw").resolve()):
                raise ValueError(f"Unlisted or unsafe snapshot path: {name}")
            data = path.read_bytes() if path.exists() else archive.read(name)
            if hashlib.sha256(data).hexdigest() != expected[name]:
                raise ValueError(f"Changed snapshot: {name}")
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
    return len(archive.namelist())


def main():
    count = restore_snapshots()
    with (ROOT / "config/research_sources_2026-09-10.tsv").open() as handle:
        row = next(row for row in csv.DictReader(handle, delimiter="\t") if row["source_id"] == "triculamin_bgc_archive")
    path = ROOT / row["local_path"]
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".download")
        subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--silent", "--show-error",
                        "--max-time", "300", row["url"], "--output", str(temporary)], check=True)
        if hashlib.sha256(temporary.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError("Published literature archive changed")
        temporary.rename(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise ValueError("Local literature archive differs from pinned source")
    print(f"Verified {count} original metadata snapshots and the literature archive")


if __name__ == "__main__":
    main()
