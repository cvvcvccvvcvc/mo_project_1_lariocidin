#!/usr/bin/env python3
"""Restore archived public metadata; fetch the pinned literature archive.

The ZIP contains original API/XML/TSV snapshots, checked against the existing
source manifests. Replaying an old selection must not query today's live API.
"""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def restore_snapshots(root=ROOT, diversity=False):
    expected = {}
    manifests = sorted((root / "inputs/provenance").glob("diversity_*download_manifest.json")) if diversity else [
        root / "inputs/provenance" / name for name in ("s3_download_manifest.json", "s5_download_manifest.json")]
    for path in manifests:
        for row in json.loads(path.read_text()).values():
            expected[row["local_path"]] = row["sha256"]
    if not diversity:
        for name, row in json.loads((root / "inputs/provenance/maize_biosample_downloads.json").read_text()).items():
            expected[f"data/raw/s3/{name}"] = row["sha256"]
        with (root / "inputs/selections/research_sources_2026-09-10.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                expected[row["local_path"]] = row["sha256"]
    archive_name = "diversity_metadata_snapshots.zip" if diversity else "source_snapshots.zip"
    with zipfile.ZipFile(root / "inputs/archives" / archive_name) as archive:
        for name in archive.namelist():
            path = root / name
            if name not in expected or not path.resolve().is_relative_to((root / "data/raw").resolve()):
                raise ValueError(f"Unlisted or unsafe snapshot path: {name}")
            missing = not path.exists()
            target = path.with_suffix(path.suffix + ".restore") if missing else path
            if missing:
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(name) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output)
            with target.open("rb") as handle:
                if hashlib.file_digest(handle, "sha256").hexdigest() != expected[name]:
                    raise ValueError(f"Changed snapshot: {name}")
            if missing:
                target.rename(path)
    return len(archive.namelist())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--diversity", action="store_true", help="Restore the separate diversity-stage metadata archive")
    args = parser.parse_args()
    count = restore_snapshots(diversity=args.diversity)
    if args.diversity:
        print(f"Verified {count} diversity metadata/model snapshots")
        return
    with (ROOT / "inputs/selections/research_sources_2026-09-10.tsv").open() as handle:
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
