#!/usr/bin/env python3
"""Preserve catalogue selection while auditing primary BioSample provenance.

Cell culture is affirmative contrary evidence, whereas absence of a metagenome
flag is not proof of culture. Unknown/unretrieved samples remain unresolved.
"""
import json
import subprocess
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from urllib.parse import urlencode

import pandas as pd

from maize_genomes import digest
from profile_search import ROOT


def validate_accessions(samples, requested):
    received = [sample.get("accession") for sample in samples]
    if len(received) != len(set(received)) or set(received) != set(requested):
        raise ValueError("BioSample response does not match requested accession set")


def main():
    table = pd.read_csv(ROOT / "config/maize_s3_genomes.tsv", sep="\t")
    scope = table[table.in_scope]
    # Preserve source-specific XML and an unresolved category rather than infer
    # MAG/isolate status from the organism name or an empty attribute.
    accessions = sorted({s for s in scope.Sample_accession if s.startswith("SAMN")})
    manifest_path = ROOT / "config/maize_biosample_downloads.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    sample_data = {}
    for offset in range(0, len(accessions), 50):
        batch = accessions[offset:offset+50]
        path = ROOT / f"data/raw/s3/biosample_batch_{offset//50+1}.xml"
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
            "db": "biosample", "id": ",".join(batch), "retmode": "xml"})
        key = path.name
        old = manifest.get(key)
        if not path.exists():
            temporary = path.with_suffix(".xml.download")
            subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--silent", "--show-error",
                            "--max-time", "90", url, "--output", str(temporary)], check=True)
            if old and digest(temporary) != old["sha256"]:
                raise ValueError(f"Changed source {key}")
            temporary.rename(path)
            time.sleep(.4)
        checksum = digest(path)
        if old:
            assert old["sha256"] == checksum and old["url"] == url, key
        manifest[key] = old or {"url": url, "sha256": checksum, "bytes": path.stat().st_size,
                               "requested_accessions": batch, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        root = ET.parse(path).getroot()
        assert root.tag == "BioSampleSet", f"Unexpected API response: {root.tag}"
        validate_accessions(root.findall("BioSample"), batch)
        for sample in root.findall("BioSample"):
            attributes = {a.get("harmonized_name", a.get("attribute_name")): a.text or ""
                          for a in sample.findall("./Attributes/Attribute")}
            package = sample.findtext("Package", "")
            sample_type = attributes.get("sample_type", "")
            metagenomic = attributes.get("metagenomic", "")
            if sample_type.lower().strip() == "cell culture":
                status = "culture_conflicts_with_MAG_label"
            elif (metagenomic.lower() in ("true", "yes") or "MIMAG" in package.upper()
                  or "metagenome-assembled genome" in sample.findtext("./Description/Title", "").lower()
                  or attributes.get("environmental_package", "").lower() == "metagenome"):
                status = "primary_metadata_supports_metagenomic_origin"
            else:
                status = "retrieved_origin_not_resolved"
            sample_data[sample.get("accession")] = {
                "provenance_status": status, "biosample_package": package,
                "sample_type": sample_type, "metagenomic": metagenomic,
                "strain": attributes.get("strain", ""),
                "isolation_source": attributes.get("isolation_source", ""),
                "host": attributes.get("host", ""), "attributes_json": json.dumps(attributes, sort_keys=True),
            }
    ena_accessions = sorted({s for s in scope.Sample_accession if s.startswith("SAMEA")})
    for offset in range(0, len(ena_accessions), 50):
        batch = ena_accessions[offset:offset+50]
        path = ROOT / f"data/raw/s3/ena_sample_batch_{offset//50+1}.xml"
        url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(batch)
        old = manifest.get(path.name)
        if not path.exists():
            temporary = path.with_suffix(".xml.download")
            subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--silent", "--show-error",
                            "--max-time", "90", url, "--output", str(temporary)], check=True)
            if old and digest(temporary) != old["sha256"]:
                raise ValueError(f"Changed source {path.name}")
            temporary.rename(path)
        checksum = digest(path)
        if old:
            assert old["sha256"] == checksum and old["url"] == url
        manifest[path.name] = old or {"url": url, "sha256": checksum, "bytes": path.stat().st_size,
                                     "requested_accessions": batch, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        root = ET.parse(path).getroot()
        assert root.tag == "SAMPLE_SET", f"Unexpected ENA response: {root.tag}"
        validate_accessions(root.findall("SAMPLE"), batch)
        for sample in root.findall("SAMPLE"):
            attributes = {a.findtext("TAG"): a.findtext("VALUE", "")
                          for a in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
            title = sample.findtext("TITLE", "")
            sample_type = attributes.get("sample_type", attributes.get("sample type", ""))
            status = "retrieved_origin_not_resolved"
            if sample_type.lower().strip() == "cell culture":
                status = "culture_conflicts_with_MAG_label"
            elif "metagenome-assembled genome" in title.lower():
                status = "primary_metadata_supports_metagenomic_origin"
            sample_data[sample.get("accession")] = {
                "provenance_status": status, "biosample_package": "ENA_SAMPLE",
                "sample_type": sample_type, "metagenomic": "", "strain": attributes.get("strain", ""),
                "isolation_source": attributes.get("isolation_source", ""),
                "host": attributes.get("host", ""), "attributes_json": json.dumps({"title": title, **attributes}, sort_keys=True),
            }
    rows = [{"Genome": row.Genome, "Sample_accession": row.Sample_accession,
             "Study_accession": row.Study_accession, "catalogue_Genome_type": row.Genome_type,
             **sample_data.get(row.Sample_accession, {"provenance_status": "primary_origin_not_audited"})}
            for row in scope.itertuples()]
    result = pd.DataFrame(rows)
    result.to_csv(ROOT / "results/tables/maize_provenance_audit.tsv", sep="\t", index=False)
    print(result.provenance_status.value_counts().to_string())


if __name__ == "__main__":
    main()
