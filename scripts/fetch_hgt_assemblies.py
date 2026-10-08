#!/usr/bin/env python3
"""Resolve exact BGC-carrier assemblies and cache bounded same-species controls."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
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


def linked_assemblies(cache: Cache, sample_id: str, nuccore_uids: list[str]) -> list[dict]:
    assembly_uids: list[str] = []
    for index, uid in enumerate(nuccore_uids):
        links = cache.get_json(
            f"{sample_id}__assembly_link_{index}",
            "elink.fcgi",
            {"dbfrom": "nuccore", "db": "assembly", "id": uid, "retmode": "json"},
        )
        for linkset in links.get("linksets", []):
            for database in linkset.get("linksetdbs", []):
                if database.get("dbto") == "assembly":
                    assembly_uids.extend(str(value) for value in database.get("links", []))
    return assembly_summaries(cache, sample_id, list(dict.fromkeys(assembly_uids)))


def resolve_contig(cache: Cache, sample_id: str, accession: str) -> tuple[dict, str]:
    search = cache.get_json(
        f"{sample_id}__nuccore_search",
        "esearch.fcgi",
        {"db": "nuccore", "term": f"{accession}[Accession]", "retmode": "json", "retmax": 5},
    )
    nuccore_uids = result_uids(search)
    summaries = linked_assemblies(cache, sample_id, nuccore_uids) if nuccore_uids else []
    method = "exact_nucleotide_accession"
    if not summaries:
        bare = accession.removeprefix("NZ_")
        match = re.match(r"^((?:[A-Z]{6}|[A-Z]{4})\d{2})\d+", bare)
        if not match:
            raise ValueError(f"No exact Nucleotide/Assembly link for {accession}")
        project = match.group(1)
        project_search = cache.get_json(
            f"{sample_id}__wgs_assembly_search",
            "esearch.fcgi",
            {"db": "assembly", "term": f"{project}[All Fields]", "retmode": "json", "retmax": 20},
        )
        summaries = assembly_summaries(
            cache, f"{sample_id}__wgs", result_uids(project_search)
        )
        method = f"wgs_project_fallback:{project}"
    if not summaries:
        raise ValueError(f"No Assembly link for {accession}")
    if len(summaries) > 1:
        raise ValueError(
            f"Ambiguous Assembly resolution for {accession}: "
            f"{[row.get('assemblyaccession') for row in summaries]}"
        )
    summaries.sort(key=lambda row: (not bool(row.get("ftppath_refseq")), row.get("assemblyaccession", "")))
    return summaries[0], method


def named_species(name: str) -> bool:
    lowered = name.lower()
    return bool(name) and " sp." not in lowered and "uncultured" not in lowered and " bacterium" not in lowered


def control_rank(row: dict) -> tuple:
    category = str(row.get("refseq_category", "")).lower()
    category_rank = 0 if "reference" in category else 1 if "representative" in category else 2
    status = str(row.get("assemblystatus", "")).lower()
    status_rank = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}.get(status, 4)
    return (category_rank, not bool(row.get("ftppath_refseq")), status_rank, row.get("assemblyaccession", ""))


def find_control(cache: Cache, sample_id: str, carrier: dict) -> dict | None:
    species = str(carrier.get("speciesname", ""))
    if not named_species(species):
        return None
    query = f'"{species}"[Organism] AND latest[filter]'
    search = cache.get_json(
        f"{sample_id}__species_assembly_search_v2",
        "esearch.fcgi",
        {"db": "assembly", "term": query, "retmode": "json", "retmax": 500},
    )
    summaries = assembly_summaries(cache, f"{sample_id}__species_v2", result_uids(search))
    carrier_accession = carrier.get("assemblyaccession")
    carrier_biosample = carrier.get("biosampleaccn")
    choices = [
        row for row in summaries
        if row.get("assemblyaccession") != carrier_accession
        and (not carrier_biosample or row.get("biosampleaccn") != carrier_biosample)
        and row.get("ftppath_refseq")
    ]
    return sorted(choices, key=control_rank)[0] if choices else None


def gbff_url(summary: dict) -> str:
    ftp = summary.get("ftppath_refseq") or summary.get("ftppath_genbank")
    if not ftp:
        raise ValueError(f"No FTP path for {summary.get('assemblyaccession')}")
    base = str(ftp).replace("ftp://", "https://")
    return f"{base}/{base.rsplit('/', 1)[-1]}_genomic.gbff.gz"


def validate_gbff(path: Path) -> tuple[int, int, dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records in {path}")
    cds = sum(1 for record in records for feature in record.features if feature.type == "CDS")
    source_values = {key: set() for key in ("strain", "isolate", "culture_collection", "type_material")}
    for record in records:
        for feature in record.features:
            if feature.type == "source":
                for key in source_values:
                    source_values[key].update(feature.qualifiers.get(key, []))
    metadata = {key: "|".join(sorted(values)) for key, values in source_values.items()}
    return len(records), cds, metadata


def output_row(
    sample_id: str, role: str, source_contig: str, summary: dict, path: Path,
    resolution_method: str,
) -> dict:
    record_count, cds_count, source_metadata = validate_gbff(path)
    return {
        "sample_id": sample_id,
        "genome_role": role,
        "source_contig_accession": source_contig,
        "resolution_method": resolution_method,
        "assembly_accession": summary.get("assemblyaccession", ""),
        "organism": summary.get("speciesname", ""),
        "biosample": summary.get("biosampleaccn", ""),
        "refseq_category": summary.get("refseq_category", ""),
        "assembly_level": summary.get("assemblystatus", ""),
        "source_strain": source_metadata["strain"],
        "source_isolate": source_metadata["isolate"],
        "culture_collection": source_metadata["culture_collection"],
        "type_material": source_metadata["type_material"],
        "local_gbff": str(path.relative_to(ROOT)),
        "record_count": record_count,
        "cds_count": cds_count,
        "sha256": digest(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    selected = selected.loc[selected["phylum"] == "Pseudomonadota"].copy()
    if len(selected) != 14 or selected["sample_id"].duplicated().any():
        raise ValueError("Expected 14 unique Pseudomonadota carriers")

    cache = Cache(args.raw_dir, args.manifest)
    rows: list[dict] = []
    for carrier in selected.itertuples(index=False):
        summary, resolution_method = resolve_contig(cache, carrier.sample_id, carrier.accession)
        carrier_path = cache.download(
            f"carrier__{carrier.sample_id}", gbff_url(summary), ".gbff.gz"
        )
        rows.append(output_row(
            carrier.sample_id, "carrier", carrier.accession, summary, carrier_path,
            resolution_method,
        ))
        control = find_control(cache, carrier.sample_id, summary)
        if control:
            control_id = f"control__{carrier.sample_id}"
            control_path = cache.download(control_id, gbff_url(control), ".gbff.gz")
            rows.append(output_row(
                control_id, "same_species_control", "", control, control_path,
                "same_species_rule",
            ))
        print(carrier.sample_id, summary.get("assemblyaccession"), "control" if control else "no_control", flush=True)

    table = pd.DataFrame(rows).sort_values(["genome_role", "sample_id"])
    if table.loc[table["genome_role"] == "carrier", "sample_id"].nunique() != 14:
        raise ValueError("Carrier assembly resolution is incomplete")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    print(f"resolved {sum(table.genome_role == 'carrier')} carriers and "
          f"{sum(table.genome_role == 'same_species_control')} controls")


if __name__ == "__main__":
    main()
