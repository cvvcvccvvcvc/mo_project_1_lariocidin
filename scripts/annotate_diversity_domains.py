#!/usr/bin/env python3
"""Addressed, local Pfam annotation with published gathering thresholds.

This small domain panel tests particular role/AsnB hypotheses, not the whole
Pfam database. A negative panel result cannot exclude a role or an RRE.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from collections import defaultdict

from Bio import SeqIO
from fetch_hgt_assemblies import Cache, ROOT

PANEL = ("PF00005", "PF00310", "PF00583", "PF00664", "PF00733",
         "PF03415", "PF05402", "PF13471", "PF13537")


def main():
    processed, tables = ROOT / "data/processed", ROOT / "results/tables"
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "maize_direct_candidates.faa", "fasta")}
    original = {r.id: str(r.seq) for r in SeqIO.parse(processed / "all_cds.faa", "fasta")}
    with (tables / "pseudomonadota_core_roles.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequences[row["feature_id"]] = original[row["feature_id"]]
    with (tables / "maize_supported_core_proteins.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequences[row["protein_id"]] = row["sequence"]
    sequences.update({identifier: sequence for identifier, sequence in original.items()
                      if identifier.startswith("paenibacillus_m2__")})
    grouped = defaultdict(list)
    for identifier, sequence in sequences.items():
        grouped[sequence].append(identifier)
    query = "".join(f">{ids[0]}\n{sequence}\n" for sequence, ids in grouped.items())
    query_path = processed / "diversity_domain_queries.faa"
    query_path.write_text(query)
    cache = Cache(ROOT / "data/raw/diversity_2026-10-12", ROOT / "config/diversity_pfam_download_manifest.json")
    model_path = processed / "diversity_role_domains.hmm"
    versions = []
    with model_path.open("w") as dest:
        for accession in PANEL:
            url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/?annotation=hmm"
            path = cache.download(accession, url, ".hmm.gz")
            with gzip.open(path, "rt") as source:
                text = source.read()
            if not text.startswith("HMMER3/") or f"ACC   {accession}." not in text or "\nGA " not in text:
                raise ValueError(f"Expected a Pfam HMM with gathering thresholds: {accession}")
            versions.append(next(line.split()[1] for line in text.splitlines() if line.startswith("ACC ")))
            dest.write(text)
    config = {"method": "hmmsearch --cut_ga", "panel_versions": versions,
              "query_sha256": hashlib.sha256(query.encode()).hexdigest(),
              "models_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
              "query_count": len(grouped), "occurrences": len(sequences),
              "limitation": "Targeted panel only; negative result does not exclude domains absent from this panel or divergent RREs."}
    (ROOT / "config/diversity_domain_search.json").write_text(json.dumps(config, indent=2) + "\n")
    domtbl = tables / "diversity_pfam_hits.domtbl"
    command = [str(ROOT / ".env/bin/hmmsearch"), "--cpu", "4", "--seed", "42", "--noali", "--cut_ga",
               "--domtblout", str(domtbl), str(model_path), str(query_path)]
    with (tables / "diversity_pfam_run.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    domains = defaultdict(list)
    for line in domtbl.read_text().splitlines():
        if line and not line.startswith("#"):
            f = line.split()
            domains[f[0]].append(f)
    output = []
    for sequence, ids in grouped.items():
        hits = domains[ids[0]]
        for identifier in ids:
            output.append({"protein_id": identifier, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                           "domains": ";".join(sorted({f[3] for f in hits})),
                           "domain_accessions": ";".join(sorted({f[4] for f in hits})),
                           "hit_types": "Pfam_gathering_threshold" if hits else "",
                           "status": "domain_hits" if hits else "no_panel_hit_not_proof_of_no_domain"})
    with (tables / "diversity_domain_evidence.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(json.dumps({"unique_queries": len(grouped), "protein_occurrences": len(output),
                      "unique_queries_with_hits": sum(bool(domains[ids[0]]) for ids in grouped.values())}))


if __name__ == "__main__":
    main()
