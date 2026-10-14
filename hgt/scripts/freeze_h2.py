#!/usr/bin/env python3
"""Freeze H2 provenance, standard HMMs, and NCBI taxonomy; no analysis."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import time
from urllib.request import Request, urlopen
from urllib.parse import urlencode

from prepare_h1 import ROOT, TABLES, read_tsv, write_tsv

H2 = ROOT / "inputs/H2_2026-10-13"
REFERENCES = H2 / "references"
MANIFEST = ROOT / "config/H2_input_manifest.json"
GTOTREE_COMMIT = "9ce125f69366b0fea043df5e91380fe8254b502e"
PFAMS = ["PF00733", "PF00310", "PF05402", "PF13471", "PF03412", "PF00583", "PF13302", "PF13508", "PF00005", "PF00664"]


def download(pair):
    name, url = pair
    path = REFERENCES / name
    if not path.exists():
        for attempt in range(3):
            try:
                with urlopen(Request(url, headers={"User-Agent": "lariocidin-HGT-H2/1.0"}), timeout=90) as response:
                    data = response.read()
                    headers = dict(response.headers)
                path.write_bytes(data)
                path.with_suffix(path.suffix + ".headers.json").write_text(json.dumps(headers, indent=2) + "\n")
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(3 * (attempt + 1))
    data = path.read_bytes()
    return {"path": str(path.relative_to(ROOT)), "url": url, "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def main():
    REFERENCES.mkdir(parents=True, exist_ok=True)
    copied = []
    sources = [TABLES / ("H1_" + name + ".tsv") for name in ["carrier_catalogue", "publication_loci", "publication_genes", "C_role_evidence", "component_candidates", "gene_genome_correspondence", "locus_genome_correspondence", "genome_quality", "relative_panel", "physical_loci"]]
    sources += [ROOT / "results/H1_sequences" / name for name in ["proteins.faa", "coding_sequences.fna", "loci.fna"]]
    sources += [ROOT / "config" / name for name in ["H1_lineage_review.tsv", "H1_named_anchor_review.tsv"]]
    for source in sources:
        destination = H2 / "H1" / source.name
        destination.parent.mkdir(exist_ok=True)
        if destination.exists() and hashlib.sha256(destination.read_bytes()).digest() != hashlib.sha256(source.read_bytes()).digest():
            raise ValueError(f"Existing H2 snapshot would change: {destination}")
        if not destination.exists():
            shutil.copyfile(source, destination)
        copied.append({"path": str(destination.relative_to(ROOT)), "source": str(source.relative_to(ROOT)),
                       "sha256": hashlib.sha256(destination.read_bytes()).hexdigest(), "bytes": destination.stat().st_size})
    base = f"https://raw.githubusercontent.com/AstrobioMike/GToTree/{GTOTREE_COMMIT}/"
    downloads = [("GToTree_HMM_sources.tsv", base + "hmm_sets/hmm-sources-and-info.tsv"),
                 ("GToTree.sh", base + "bin/GToTree"),
                 ("GToTree_User_Guide.html", "https://github.com/AstrobioMike/GToTree/wiki/User-Guide"),
                 ("GToTree_considerations.html", "https://github.com/AstrobioMike/GToTree/wiki/things-to-consider")]
    entries = [download(pair) for pair in downloads]
    source_info = read_tsv(REFERENCES / "GToTree_HMM_sources.tsv")
    bacterial = next(r for r in source_info if r["file"] == "Bacteria.hmm")
    entries.append(download(("Bacteria.hmm", bacterial["link"])))
    bacterial_bytes = (REFERENCES / "Bacteria.hmm").read_bytes()
    if hashlib.md5(bacterial_bytes).hexdigest() != bacterial["md5"]:
        raise ValueError("GToTree author's Bacteria HMM checksum differs")
    if sum(line.startswith(b"NAME ") for line in bacterial_bytes.splitlines()) != 74:
        raise ValueError("Expected published 74-marker panel")
    with ThreadPoolExecutor(max_workers=3) as pool:
        entries.extend(pool.map(download, [(x + ".hmm.gz", f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{x}/?annotation=hmm") for x in PFAMS]))
    profiles = []
    for accession in PFAMS:
        data = gzip.decompress((REFERENCES / (accession + ".hmm.gz")).read_bytes())
        if not data.startswith(b"HMMER3/") or b"GA " not in data:
            raise ValueError(f"Not a curated HMMER profile: {accession}")
        profiles.append(data)
    (REFERENCES / "BGC_domains.hmm").write_bytes(b"\n".join(profiles))
    taxids = sorted({x for r in read_tsv(TABLES / "H1_genome_quality.tsv") for x in r["taxid"].split(";") if x.isdigit()}, key=int)
    for start in range(0, len(taxids), 100):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({"db": "taxonomy", "id": ",".join(taxids[start:start + 100]), "retmode": "xml"})
        entries.append(download((f"NCBI_taxonomy_{start // 100:02}.xml", url)))
        time.sleep(0.4)
    assembly_tables = [TABLES / "H1_assembly_inputs.tsv", TABLES / "H1_supplement_assembly_inputs.tsv", TABLES / "H1_Nature_4BL1_assembly_inputs.tsv"]
    inputs = [r for path in assembly_tables for r in read_tsv(path)]
    write_tsv(H2 / "H1/H1_all_assembly_inputs.tsv", inputs)
    genome_entries = []
    for row in inputs:
        path = ROOT / row["path"]
        genome_entries.append({"genome_id": row["genome_id"], "path": row["path"],
                               "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "scope": "H2 frozen H1 inputs and primary method references", "GToTree_commit": GTOTREE_COMMIT,
                "genomes": genome_entries, "copied": copied, "downloaded": entries,
                "derived_profiles": {"path": str((REFERENCES / "BGC_domains.hmm").relative_to(ROOT)),
                                     "sha256": hashlib.sha256((REFERENCES / "BGC_domains.hmm").read_bytes()).hexdigest()}}
    if MANIFEST.exists():
        old = json.loads(MANIFEST.read_text())
        for key in ["genomes", "copied", "downloaded", "derived_profiles"]:
            if old[key] != manifest[key]:
                raise ValueError(f"Frozen H2 manifest differs: {key}")
    else:
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"genomes": len(inputs), "primary_references": len(entries), "SCG_models": 74, "component_domain_models": len(PFAMS)}))


if __name__ == "__main__":
    main()
