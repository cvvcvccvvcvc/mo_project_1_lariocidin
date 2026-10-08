#!/usr/bin/env python3
"""Cache and normalize primary BioSample metadata for the bounded S5 genomes."""

from __future__ import annotations

import argparse
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from fetch_hgt_assemblies import Cache, result_uids


FIELDS = (
    "sample_type", "isolation_source", "geo_loc_name", "country", "lat_lon",
    "collection_date", "host", "host_disease", "env_broad_scale", "env_local_scale",
    "env_medium", "environmental_sample",
)


def values_by_name(sample: ET.Element) -> dict[str, list[str]]:
    values: dict[str, list[str]] = {}
    for attribute in sample.findall("./Attributes/Attribute"):
        name = attribute.attrib.get("harmonized_name") or attribute.attrib.get("attribute_name", "")
        name = name.strip().lower().replace(" ", "_")
        value = (attribute.text or "").strip()
        if name and value:
            values.setdefault(name, []).append(value)
    return values


def environment_labels(text: str) -> str:
    lowered = text.lower()
    rules = {
        "human_associated": ("human", "patient", "clinical", "blood", "wound", "stool", "gut"),
        "plant_or_soil": ("plant", "soil", "root", "rhizosphere", "palm", "agricultur"),
        "invertebrate_associated": ("nematode", "insect", "steinernema", "heterorhabdit"),
        "aquatic": ("water", "marine", "aquatic", "seawater", "freshwater"),
    }
    labels = [label for label, words in rules.items() if any(word in lowered for word in words)]
    return "|".join(labels) if labels else "unknown"


def record_type(title: str, values: dict[str, list[str]]) -> str:
    text = " ".join([title, *[item for group in values.values() for item in group]]).lower()
    sample_type = "|".join(values.get("sample_type", [])).lower()
    if "cell culture" in sample_type:
        return "cell_culture"
    if any(word in text for word in ("metagenome", "metagenomic", " mag", "concoct")):
        return "metagenomic_assembly_source"
    return "isolate_or_unspecified"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--attributes-out", type=Path, required=True)
    args = parser.parse_args()

    assemblies = pd.read_csv(args.assemblies, sep="\t", keep_default_na=False)
    cache = Cache(args.raw_dir, args.manifest)
    summaries = []
    attributes = []
    for row in assemblies.itertuples(index=False):
        if not row.biosample:
            summaries.append({"sample_id": row.sample_id, "biosample": "", "metadata_status": "missing_accession"})
            continue
        search = cache.get_json(
            f"biosample_search__{row.sample_id}",
            "esearch.fcgi",
            {"db": "biosample", "term": f"{row.biosample}[Accession]", "retmode": "json", "retmax": 5},
        )
        uids = result_uids(search)
        if len(uids) != 1:
            raise ValueError(f"Expected one BioSample UID for {row.biosample}, found {uids}")
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "biosample", "id": uids[0], "retmode": "xml"
        })
        path = cache.download(f"biosample_v2__{row.sample_id}", url, ".xml")
        root = ET.parse(path).getroot()
        sample = root.find("./BioSample")
        if sample is None:
            raise ValueError(f"No BioSample in {path}")
        if sample.attrib.get("accession") != row.biosample:
            raise ValueError(
                f"BioSample accession mismatch for {row.sample_id}: "
                f"requested {row.biosample}, received {sample.attrib.get('accession')}"
            )
        title = sample.findtext("./Description/Title", default="").strip()
        organism = sample.findtext("./Description/Organism/OrganismName", default="").strip()
        values = values_by_name(sample)
        source_text = " ".join([
            title,
            *values.get("isolation_source", []), *values.get("host", []),
            *values.get("env_broad_scale", []), *values.get("env_local_scale", []),
            *values.get("env_medium", []),
        ])
        summary = {
            "sample_id": row.sample_id,
            "genome_role": row.genome_role,
            "assembly_accession": row.assembly_accession,
            "biosample": row.biosample,
            "metadata_status": "retrieved",
            "title": title,
            "organism": organism,
            "record_type": record_type(title, values),
            "reported_environment_category": environment_labels(source_text),
        }
        for field in FIELDS:
            summary[field] = "|".join(dict.fromkeys(values.get(field, [])))
        summaries.append(summary)
        for name, group in sorted(values.items()):
            for value in dict.fromkeys(group):
                attributes.append({
                    "sample_id": row.sample_id,
                    "biosample": row.biosample,
                    "attribute_name": name,
                    "attribute_value": value,
                })
        time.sleep(0.36)

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(args.summary_out, sep="\t", index=False)
    pd.DataFrame(attributes).to_csv(args.attributes_out, sep="\t", index=False)
    print(summary_table.groupby(["record_type", "reported_environment_category"]).size().to_string())


if __name__ == "__main__":
    main()
