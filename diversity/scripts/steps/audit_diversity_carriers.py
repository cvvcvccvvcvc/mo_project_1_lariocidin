#!/usr/bin/env python3
"""Check primary provenance for supported carriers, without changing BGC calls."""
import argparse
import json
import urllib.parse
import xml.etree.ElementTree as ET

import pandas as pd
from audit_maize_provenance import validate_accessions
from fetch_hgt_assemblies import Cache, ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="diversity_species")
    parser.add_argument("--refseq-metadata", default="data/processed/discovery/diversity_refseq_metadata.tsv")
    args = parser.parse_args()
    tables = ROOT / "results/tables"
    primary = pd.read_csv(tables / f"{args.prefix}_BGC_catalogue.tsv", sep="\t", keep_default_na=False)
    references = pd.read_csv(ROOT / args.refseq_metadata, sep="\t", keep_default_na=False).set_index("assembly")
    soil = primary[primary.source != "RefSeq"]
    cache = Cache(ROOT / "data/raw/diversity_2026-10-13", ROOT / "inputs/provenance/diversity_carrier_download_manifest.json")
    samples = {}
    for provider, prefix in (("NCBI", "SAMN"), ("ENA", "SAMEA")):
        accessions = sorted(set(x for x in soil.biosample if x.startswith(prefix)))
        if not accessions:
            continue
        if provider == "NCBI":
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db": "biosample", "id": ",".join(accessions), "retmode": "xml"})
        else:
            url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(accessions)
        path = cache.download(args.prefix + "_supported_" + provider.lower(), url, ".xml")
        root = ET.parse(path).getroot()
        elements = root.findall("BioSample" if provider == "NCBI" else "SAMPLE")
        validate_accessions(elements, accessions)
        for sample in elements:
            if provider == "NCBI":
                attrs = {a.get("harmonized_name") or a.get("attribute_name"): a.text or "" for a in sample.findall("./Attributes/Attribute")}
                title = sample.findtext("./Description/Title", "")
                package = sample.findtext("./Package", "")
                organism = sample.findtext("./Description/Organism/OrganismName", "")
            else:
                attrs = {a.findtext("TAG"): a.findtext("VALUE", "") for a in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")}
                title = sample.findtext("TITLE", "")
                package = "ENA_SAMPLE"
                organism = sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME", "")
            attrs = {str(k).strip().lower().replace(" ", "_"): v for k, v in attrs.items()}
            samples[sample.get("accession")] = {"metadata_provider": provider, "title": title, "package": package,
                                                "primary_organism": organism, "attributes_json": json.dumps(attrs, sort_keys=True), **attrs}
    output = []
    for row in primary.to_dict("records"):
        if row["source"] == "RefSeq":
            ref = references.loc[row["genome_id"]].to_dict()
            data = {"metadata_provider": "frozen_NCBI_Datasets_BioSample", "title": ref["biosample_title"],
                    "package": ref["biosample_package"], "sample_type": ref["sample_type"],
                    "metagenomic": ref["metagenomic"], "metagenome_source": ref["metagenome_source"],
                    "isolation_source": ref["isolation_source"], "host": ref["host"], "strain": ref["strain"]}
        else:
            data = samples.get(row["biosample"], {"metadata_provider": "unresolved"})
        sample_type = data.get("sample_type", "").lower().strip()
        package = data.get("package", "").upper()
        status = "retrieved_origin_not_resolved"
        basis = "retrieved_fields_do_not_resolve_origin"
        if sample_type in ("cell culture", "pure culture", "axenic culture") or package.startswith("MIGS.BA."):
            status = "culture_conflicts_with_MAG_label" if row["source_genome_type"] == "MAG" else "primary_metadata_supports_culture"
            basis = "explicit_culture_sample_type" if sample_type else "NCBI_MIGS_cultured_bacteria_archaea_package"
        elif ("MIMAG" in package or data.get("metagenomic", "").lower() in ("yes", "true", "1")
              or sample_type in ("metagenomic assembly", "metagenome-assembled genome")
              or data.get("title", "").lower().startswith("mag for ")
              or "metagenome-assembled genome" in data.get("title", "").lower()
              or data.get("environmental_package", "").lower() == "metagenome"):
            status = "primary_metadata_supports_metagenomic_origin"
            basis = "MIMAG_package_or_explicit_metagenomic_attribute_or_MAG_title"
        elif sample_type == "single cell":
            status = "primary_metadata_supports_single_cell_origin"
            basis = "explicit_single_cell_sample_type"
        output.append({"locus_id": row["locus_id"], "genome_id": row["genome_id"], "biosample": row["biosample"],
                       "catalogue_genome_type": row["source_genome_type"], "provenance_status": status,
                       "provenance_basis": basis, **data})
    result = pd.DataFrame(output).fillna("")
    result.to_csv(tables / f"{args.prefix}_carrier_provenance.tsv", sep="\t", index=False)
    print(result.provenance_status.value_counts().to_dict())


if __name__ == "__main__":
    main()
