#!/usr/bin/env python3
"""Fetch and verify pinned NCBI source snapshots for diversity workflows."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path



ROOT = Path(__file__).resolve().parents[2]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "lariocidin-course-project/1.0 (NCBI-EUtils)"


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


class Cache:
    def __init__(self, raw_dir: Path, manifest_path: Path) -> None:
        self.raw_dir = raw_dir
        self.manifest_path = manifest_path
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.records = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        self.last_request = 0.0

    def _save_record(self, source_id: str, path: Path, url: str) -> None:
        record = {
            "url": url,
            "local_path": str(path.resolve().relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        }
        old = self.records.get(source_id)
        if old and any(old.get(k) != record[k] for k in ("url", "bytes", "sha256")):
            raise ValueError(f"Pinned snapshot changed: {source_id}")
        self.records[source_id] = old or record
        self.manifest_path.write_text(json.dumps(self.records, indent=2, sort_keys=True) + "\n")

    def get_json(self, source_id: str, endpoint: str, params: dict[str, object]) -> dict:
        path = self.raw_dir / f"{source_id}.json"
        url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
        if not path.exists():
            pause = 0.36 - (time.monotonic() - self.last_request)
            if pause > 0:
                time.sleep(pause)
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            self.last_request = time.monotonic()
            temporary = path.with_suffix(".json.download")
            temporary.write_bytes(payload)
            json.loads(payload)
            temporary.rename(path)
        self._save_record(source_id, path, url)
        return json.loads(path.read_text())

    def download(self, source_id: str, url: str, suffix: str) -> Path:
        path = self.raw_dir / f"{source_id}{suffix}"
        if not path.exists():
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            temporary = path.with_suffix(path.suffix + ".download")
            with urllib.request.urlopen(request, timeout=1800) as response, temporary.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            temporary.rename(path)
        self._save_record(source_id, path, url)
        return path


def result_uids(payload: dict, key: str = "esearchresult") -> list[str]:
    return [str(value) for value in payload[key].get("idlist", [])]


def assembly_summaries(cache: Cache, source_prefix: str, uids: list[str]) -> list[dict]:
    if not uids:
        return []
    payload = cache.get_json(
        f"{source_prefix}__assembly_summary",
        "esummary.fcgi",
        {"db": "assembly", "id": ",".join(uids), "retmode": "json"},
    )
    result = payload["result"]
    return [result[uid] for uid in result.get("uids", [])]


