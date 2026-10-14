#!/usr/bin/env python3
"""Freeze primary NCBI metadata and genomes for the fixed H1 publication set."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from prepare_h1 import ROOT, H1, TABLES, read_tsv, write_tsv

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RAW = H1 / "ncbi"
MANIFEST = ROOT / "config/H1/H1_ncbi_download_manifest.json"
NUCLEOTIDE = re.compile(r"^(?:NZ_|NC_|NT_|NW_)?(?:[A-Z]{1,2}\d{5,}|[A-Z]{4}\d{8,}|[A-Z]{6}\d{9,})(?:\.\d+)?$")


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


class Snapshot:
    def __init__(self):
        RAW.mkdir(parents=True, exist_ok=True)
        self.entries = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
        self.last_request = 0.0

    def record(self, key, path, url):
        data = path.read_bytes()
        row = {"url": url, "local_path": str(path.relative_to(ROOT)), "bytes": len(data),
               "sha256": hashlib.sha256(data).hexdigest(), "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        previous = self.entries.get(key)
        if previous and any(previous.get(k) != row[k] for k in ["url", "sha256", "bytes"]):
            raise ValueError(f"Snapshot changed: {key}")
        self.entries[key] = previous or row
        atomic_json(MANIFEST, self.entries)

    def json(self, key, endpoint, parameters):
        path = RAW / f"{key}.json"
        url = EUTILS + endpoint + "?" + urlencode(parameters, doseq=True)
        if not path.exists():
            time.sleep(max(0, 0.38 - (time.monotonic() - self.last_request)))
            for attempt in range(3):
                try:
                    with urlopen(Request(url, headers={"User-Agent": "lariocidin-HGT/1.0"}), timeout=60) as response:
                        data = response.read()
                    json.loads(data)
                    path.write_bytes(data)
                    self.last_request = time.monotonic()
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(5 * (attempt + 1))
        self.record(key, path, url)
        return json.loads(path.read_text())


def chunks(values, size):
    for start in range(0, len(values), size):
        yield values[start:start + size]


def fetch_genome(item):
    accession, row = item
    folder = (row.get("ftp_path") or row.get("ftppath_refseq") or row.get("ftppath_genbank") or "").rstrip("/")
    folder = folder.replace("ftp://", "https://")
    if not folder or folder == "na":
        return accession, None, "no_genome_ftp_path"
    url = folder + "/" + folder.rsplit("/", 1)[-1] + "_genomic.gbff.gz"
    path = RAW / "genomes" / f"{accession}.gbff.gz"
    path.parent.mkdir(exist_ok=True)
    if path.exists():
        return accession, (path, url), ""
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"User-Agent": "lariocidin-HGT/1.0"}), timeout=90) as response:
                temporary = path.with_suffix(path.suffix + ".download")
                with temporary.open("wb") as handle:
                    while block := response.read(1024 * 1024):
                        handle.write(block)
            temporary.replace(path)
            return accession, (path, url), ""
        except Exception as error:
            if getattr(error, "code", None) == 404 or attempt == 2:
                return accession, None, str(error)
            time.sleep(5 * (attempt + 1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    snapshot = Snapshot()
    loci = read_tsv(TABLES / "H1_publication_loci.tsv")
    accessions = sorted({r["accession_evidence"] for r in loci if NUCLEOTIDE.fullmatch(r["accession_evidence"])})
    nucleotide_rows = {}
    for index, batch in enumerate(chunks(accessions, 40)):
        payload = snapshot.json(f"nuccore_batch_{index:03}", "esummary.fcgi", {"db": "nuccore", "id": ",".join(batch), "retmode": "json"})
        for uid in payload.get("result", {}).get("uids", []):
            row = payload["result"][uid]
            if "error" not in row:
                nucleotide_rows[uid] = row
    links = {}
    for index, batch in enumerate(chunks(sorted(nucleotide_rows), 25)):
        payload = snapshot.json(f"assembly_links_{index:03}", "elink.fcgi", {"dbfrom": "nuccore", "db": "assembly", "id": batch, "retmode": "json"})
        for linkset in payload.get("linksets", []):
            linked = [str(uid) for db in linkset.get("linksetdbs", []) if db.get("dbto") == "assembly" for uid in db.get("links", [])]
            for uid in linkset.get("ids", []):
                links[str(uid)] = linked
    assembly_uids = sorted({uid for values in links.values() for uid in values})
    assemblies_by_uid = {}
    for index, batch in enumerate(chunks(assembly_uids, 100)):
        payload = snapshot.json(f"assembly_summaries_{index:03}", "esummary.fcgi", {"db": "assembly", "id": ",".join(batch), "retmode": "json"})
        for uid in payload.get("result", {}).get("uids", []):
            if "error" not in payload["result"][uid]:
                assemblies_by_uid[uid] = payload["result"][uid]
    by_accession = {}
    for uid, row in nucleotide_rows.items():
        for name in [row.get("caption", ""), row.get("accessionversion", "")]:
            by_accession[name] = uid
            by_accession[name.removeprefix("NZ_")] = uid
    resolved = []
    chosen = {}
    for locus in loci:
        evidence = locus["accession_evidence"]
        uid = by_accession.get(evidence) or by_accession.get(evidence.split(".")[0])
        candidates = [assemblies_by_uid[x] for x in links.get(uid, []) if x in assemblies_by_uid]
        candidates.sort(key=lambda r: (not r["assemblyaccession"].startswith("GCF_"), r["assemblyaccession"]))
        # Keep the link set visible. Preference does not establish sequence identity.
        selected = candidates[0] if candidates else None
        if selected:
            chosen[selected["assemblyaccession"]] = selected
        resolved.append({"locus_id": locus["locus_id"], "accession_evidence": evidence,
                         "nuccore_uid": uid or "", "nucleotide_version": nucleotide_rows.get(uid, {}).get("accessionversion", ""),
                         "linked_assemblies": ";".join(r["assemblyaccession"] for r in candidates),
                         "preferred_linked_assembly": selected["assemblyaccession"] if selected else "",
                         "resolution_status": "nucleotide_assembly_link" if selected else "not_resolved_by_nucleotide_link"})
    write_tsv(TABLES / "H1_nucleotide_assembly_links.tsv", resolved)
    atomic_json(RAW / "assembly_reports_by_accession.json", chosen)
    candidates = {r["assembly_accession"]: r for r in read_tsv(TABLES / "H1_wgs_assembly_candidates.tsv")}
    npdc_path = TABLES / "H1_npdc_assembly_candidates.tsv"
    if npdc_path.exists():
        candidates.update({r["assembly_accession"]: r for r in read_tsv(npdc_path)})
    candidates.update(chosen)
    legacy = read_tsv(ROOT / "inputs/2026-10-13/legacy/config/hgt_assemblies.tsv")
    inputs = [{"genome_id": r["sample_id"], "assembly_accession": r["assembly_accession"],
               "path": str((H1 / "legacy_assemblies" / Path(r["local_gbff"]).name).relative_to(ROOT)),
               "input_scope": r["genome_role"], "biosample": r["biosample"], "organism": r["organism"],
               "assembly_level": r["assembly_level"], "source_strain": r["source_strain"], "taxid": ""} for r in legacy]
    old_accessions = {r["assembly_accession"] for r in inputs}
    candidates = {k: v for k, v in candidates.items() if k not in old_accessions}
    atomic_json(RAW / "assembly_candidate_metadata.json", candidates)
    print(json.dumps({"target_loci": len(loci), "nucleotide_queries": len(accessions), "nuccore_resolved": len(nucleotide_rows),
                      "linked_target_loci": sum(bool(r["preferred_linked_assembly"]) for r in resolved), "genome_candidates": len(candidates)}), flush=True)
    if args.metadata_only:
        return
    errors = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(fetch_genome, pair) for pair in sorted(candidates.items())]
        for index, future in enumerate(as_completed(futures), 1):
            accession, downloaded, error = future.result()
            if error:
                errors.append({"assembly_accession": accession, "error": error})
            else:
                path, url = downloaded
                snapshot.record("genome__" + accession, path, url)
                row = candidates[accession]
                inputs.append({"genome_id": accession, "assembly_accession": accession, "path": str(path.relative_to(ROOT)),
                               "input_scope": "publication_candidate", "biosample": row.get("biosample", row.get("biosampleaccn", "")),
                               "organism": row.get("organism_name", row.get("organism", "")),
                               "assembly_level": row.get("assembly_level", row.get("assemblystatus", "")),
                               "source_strain": row.get("infraspecific_name", ""), "taxid": row.get("taxid", "")})
            if index % 25 == 0:
                print(f"Genomes completed {index}/{len(futures)}; failed {len(errors)}", flush=True)
    inputs.append({"genome_id": "M2_CP169648_reference", "assembly_accession": "", "path": "inputs/2026-10-13/base/data/raw/CP169648.1.gb",
                   "input_scope": "single_reference_chromosome", "biosample": "", "organism": "Paenibacillus sp. M2",
                   "assembly_level": "single_complete_chromosome", "source_strain": "M2", "taxid": "3341793"})
    write_tsv(TABLES / "H1_assembly_inputs.tsv", sorted(inputs, key=lambda r: r["genome_id"]))
    atomic_json(ROOT / "results/tables/H1_genome_download_errors.json", errors)
    print(json.dumps({"genome_inputs": len(inputs), "download_errors": len(errors)}), flush=True)


if __name__ == "__main__":
    main()
