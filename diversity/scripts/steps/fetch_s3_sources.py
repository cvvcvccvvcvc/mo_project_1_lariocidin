#!/usr/bin/env python3
"""Cache the explicitly listed S3 inputs; refuse changes to pinned snapshots."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    manifest = ROOT / "inputs/provenance/s3_download_manifest.json"
    records = json.loads(manifest.read_text()) if manifest.exists() else {}
    with (ROOT / "inputs/selections/baseline/s3_sources.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            path = ROOT / row["local_path"]
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                temporary = path.with_suffix(path.suffix + ".download")
                subprocess.run([
                    "curl", "--fail", "--location", "--retry", "3", "--silent",
                    "--show-error", "--max-time", "1800", row["url"],
                    "--output", str(temporary),
                ], check=True)
                digest = sha256(temporary)
                expected = records.get(row["source_id"], {}).get("sha256")
                if expected and digest != expected:
                    raise ValueError(f"Changed source: {row['source_id']}")
                temporary.rename(path)
            digest = sha256(path)
            old = records.get(row["source_id"])
            if old and (digest != old["sha256"] or row["url"] != old["url"]):
                raise ValueError(f"Snapshot mismatch: {row['source_id']}")
            records[row["source_id"]] = old or {
                **row, "bytes": path.stat().st_size, "sha256": digest,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            }
            manifest.write_text(json.dumps(records, indent=2) + "\n")
            print(row["source_id"], path.stat().st_size, digest, flush=True)


if __name__ == "__main__":
    main()
