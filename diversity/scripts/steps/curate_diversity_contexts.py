#!/usr/bin/env python3
"""Evidence tables and a conservative BGC/uncertain split for one context shard.

BLAST/domain hits are retained as evidence. They do not assign a metabolic
function to orphan C homologs or establish antibiotic production.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from profile_search import ROOT

BLAST_FIELDS = ["query_id", "subject_id", "pident", "alignment_length", "qlen", "slen", "evalue", "bitscore", "qstart", "qend", "sstart", "send"]


def read_blast(output):
    if not output.stat().st_size:
        return pd.DataFrame(columns=BLAST_FIELDS+["qcov", "scov"])
    hits = pd.read_csv(output, sep="\t", names=BLAST_FIELDS)
    hits["qcov"] = 100*(hits.qend-hits.qstart+1)/hits.qlen
    hits["scov"] = 100*(hits.send-hits.sstart+1)/hits.slen
    return hits.sort_values(["bitscore", "subject_id"], ascending=[False, True])


def blast(query, subject, output, cpu, short=False, database=False):
    command = [str(ROOT / ".env/bin/blastp"), "-query", str(query),
               "-db" if database else "-subject", str(subject), "-num_threads", str(cpu),
               "-evalue", "10" if short else "1e-3", "-max_target_seqs", "20",
               "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend sstart send", "-out", str(output)]
    if short:
        command += ["-task", "blastp-short", "-seg", "no"]
    with output.with_suffix(".log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    return read_blast(output)


def terminal_precursor_overlap(a, c, locus, window, sequence, allowed_lengths):
    """Recognize a complete upstream A with a short stop/start overlap with C."""
    overlap = min(a["end"], int(c.end))-max(a["start"], int(c.start))
    if overlap not in allowed_lengths or a["strand"] != int(c.strand):
        return False
    if (a["strand"] == 1 and not a["start"] < int(c.start) < a["end"] < int(c.end)) or (
            a["strand"] == -1 and not int(c.start) < a["start"] < int(c.end) < a["end"]):
        return False
    def coding(start, end):
        dna = window.seq[int(start)-locus["window_start"]:int(end)-locus["window_start"]]
        return dna if a["strand"] == 1 else dna.reverse_complement()
    dna = coding(a["start"], a["end"])
    c_dna = coding(c.start, c.end)
    code = CodonTable.unambiguous_dna_by_id[11]
    if len(dna)%3 or str(dna[:3]) not in code.start_codons or str(dna[-3:]) not in code.stop_codons:
        return False
    protein = str(dna[:-3].translate(table=11))
    protein = "M"+protein[1:]
    # Starts in different frames; internal translation fragments do not qualify.
    return (protein == sequence and str(c_dna[:3]) in code.start_codons
            and (abs(a["start"]-int(c.start)) if a["strand"] == 1 else abs(a["end"]-int(c.end)))%3 != 0)


def processing_pair(a, c, b1, b2, domains):
    low, high = sorted((a["start"], int(c.start)))
    for first in b1[(b1.start > low) & (b1.start < high)].to_dict("records"):
        for second in b2[(b2.start > low) & (b2.start < high)].to_dict("records"):
            if first["query_id"] != second["query_id"] or {"PqqD", "Transglut_core3"} <= domains[first["query_id"]]:
                return first, second
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen-name", required=True)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--reuse-evidence", action="store_true", help="Read saved BLAST/domain results without searching again")
    parser.add_argument("--output-name", help="Separate result name; preserves the original context evidence and catalogue")
    parser.add_argument("--protocol", type=Path, default=ROOT / "protocols/diversity_protocol.json")
    args = parser.parse_args()
    out = ROOT / "data/processed/discovery/contexts" / f"{args.screen_name}_contexts"
    tables = ROOT / "results/tables"
    protocol = json.loads(args.protocol.read_text())
    output_name = args.output_name or args.screen_name
    fragmented_c = (set(pd.read_csv(ROOT / protocol["C_integrity_review_table"], sep="\t").locus_id)
                    if protocol.get("C_integrity_review_table") else set())
    loci = pd.read_csv(out / "loci.tsv", sep="\t", keep_default_na=False)
    features = pd.read_csv(out / "features.tsv", sep="\t", keep_default_na=False)
    role_reference = pd.read_csv(ROOT / "data/processed/discovery/diversity_role_reference.tsv", sep="\t", keep_default_na=False).set_index("feature_id")
    precursor_reference = pd.read_csv(ROOT / "data/processed/literature/literature_precursor_occurrences.tsv", sep="\t", keep_default_na=False)
    precursor_groups = precursor_reference.groupby("precursor_id").source_group.agg(lambda x: set(x)).to_dict()
    cyclase_reference = pd.read_csv(ROOT / "data/processed/literature/literature_unique_cyclases.tsv", sep="\t", keep_default_na=False).set_index("cyclase_id")
    gene_sequences = {r.id: str(r.seq) for r in SeqIO.parse(out / "neighbors.faa", "fasta")}
    gene_sequences.update({r.id: str(r.seq) for r in SeqIO.parse(out / "precursor_queries.faa", "fasta")})
    c_ids = set(features.loc[features.is_C, "query_id"])
    if not args.reuse_evidence:
        SeqIO.write([r for r in SeqIO.parse(out / "neighbors.faa", "fasta") if r.id in c_ids], out / "C.faa", "fasta")

    domain_table = out / "domains.domtbl"
    if args.reuse_evidence:
        precursor_hits, role_hits, c_hits, reviewed_hits = [read_blast(out / name) for name in
            ["precursor_blast.tsv", "role_blast.tsv", "C_literature_blast.tsv", "C_reviewed_blast.tsv"]]
    else:
        precursor_hits = blast(out / "precursor_queries.faa", ROOT / "data/processed/literature/literature_unique_precursors.faa", out / "precursor_blast.tsv", args.cpu, short=True)
        role_hits = blast(out / "neighbors.faa", ROOT / "data/processed/discovery/diversity_role_reference.faa", out / "role_blast.tsv", args.cpu)
        c_hits = blast(out / "C.faa", ROOT / "data/processed/discovery/diversity_published_C.faa", out / "C_literature_blast.tsv", args.cpu)
        reviewed_hits = blast(out / "C.faa", ROOT / "data/interim/diversity_baseline/reviewed_bacterial.faa", out / "C_reviewed_blast.tsv", args.cpu, database=True)
        with (out / "domains.log").open("w") as log:
            subprocess.run([str(ROOT / ".env/bin/hmmsearch"), "--cpu", str(args.cpu), "--seed", "42", "--noali", "--cut_ga",
                            "--domtblout", str(domain_table), str(ROOT / "data/processed/discovery/diversity_role_domains.hmm"), str(out / "neighbors.faa")], stdout=log, stderr=log, check=True)
    domains = defaultdict(set)
    for line in domain_table.read_text().splitlines():
        if line and not line.startswith("#"):
            f = line.split()
            domains[f[0]].add(f[3])
    admissible_a = precursor_hits[(precursor_hits.pident >= protocol["precursor_min_identity_percent"])
                                 & (precursor_hits.alignment_length >= protocol["precursor_min_alignment_aa"])
                                 & (precursor_hits.qcov >= protocol["precursor_min_query_coverage_percent"])
                                 & (precursor_hits.scov >= protocol["precursor_min_subject_coverage_percent"])]
    if protocol.get("precursor_alignment_must_span_reference_core", False):
        sites = pd.read_csv(ROOT / protocol["precursor_core_site_table"], sep="\t")
        admissible_a = admissible_a.merge(sites, left_on="subject_id", right_on="precursor_id", how="inner")
        admissible_a = admissible_a[(admissible_a.sstart <= admissible_a.reference_core_start)
                                    & (admissible_a.send >= admissible_a.reference_core_end)]
    admissible_roles = role_hits[(role_hits.pident >= protocol["role_min_identity_percent"])
                                & (role_hits.qcov >= protocol["role_min_query_coverage_percent"])
                                & (role_hits.scov >= protocol["role_min_subject_coverage_percent"])]
    role_evidence = []
    supported_roles = defaultdict(set)
    role_hits_by_query = {q: data.to_dict("records") for q, data in admissible_roles.groupby("query_id")}
    for feature in features[features.origin == "annotated_CDS"].to_dict("records"):
        query = feature["query_id"]
        assignments = defaultdict(list)
        for hit in role_hits_by_query.get(query, []):
            ref = role_reference.loc[hit["subject_id"]]
            for role in ref["inferred_roles"].split(";"):
                assignments[role].append({**hit, "reference_basis": ref["role_evidence_basis"]})
        product = feature["product"].lower()
        for role, check in [("B1", "PqqD" in domains[query] or "pqqd" in product or "recognition element" in product),
                            ("B2", "Transglut_core3" in domains[query] or "biosynthesis b2" in product),
                            ("ABC", "ABC_tran" in domains[query]), ("GNAT", "Acetyltransf_1" in domains[query]),
                            ("C11", "Peptidase_C11" in domains[query])]:
            if check:
                supported_roles[query].add(role)
                assignments[role]  # Create a row even without a sequence hit.
        for role, matches in assignments.items():
            if matches:
                supported_roles[query].add(role)
            hit = matches[0] if matches else {}
            role_evidence.append({**feature, "proposed_role": role, "domain_panel": ";".join(sorted(domains[query])),
                                  "reference_id": hit.get("subject_id", ""), "reference_assignment_basis": hit.get("reference_basis", ""),
                                  "pident": hit.get("pident", ""), "qcov": hit.get("qcov", ""), "scov": hit.get("scov", ""),
                                  "evalue": hit.get("evalue", ""), "bitscore": hit.get("bitscore", ""),
                                  "function_status": "predicted_role_not_per_locus_function_experiment"})
    # Keep alternative starts until geometry/translation checks have run.
    # An invalid longer start must not hide a valid shorter precursor.
    a_rows = features.merge(admissible_a, on="query_id", how="inner")
    a_rows["stop_coordinate"] = a_rows["end"].where(a_rows.strand == 1, a_rows["start"])
    a_rows = a_rows.sort_values(["bitscore", "qcov", "query_id"], ascending=[False, False, True]).drop_duplicates(["locus_id", "query_id"])
    a_rows.to_csv(tables / f"{output_name}_precursor_evidence.tsv" if args.output_name else out / "precursor_evidence.tsv", sep="\t", index=False)
    c_best = c_hits.drop_duplicates("query_id").set_index("query_id")
    reviewed_best = reviewed_hits.drop_duplicates("query_id").set_index("query_id")
    annotated = features[features.origin == "annotated_CDS"]
    local_by_locus = {locus: data for locus, data in annotated.groupby("locus_id")}
    a_by_locus = {locus: data for locus, data in a_rows.groupby("locus_id")}
    windows = {r.id: r for r in SeqIO.parse(out / "windows.fna", "fasta")}
    output = []
    for locus in loci.to_dict("records"):
        local = local_by_locus[locus["locus_id"]]
        c = local[local.is_C].iloc[0]
        query = c.query_id
        candidates_a = a_by_locus.get(locus["locus_id"], a_rows.iloc[:0])
        if not protocol.get("canonical_opposite_A_with_processing_pair", False):
            candidates_a = candidates_a[candidates_a.strand == locus["C_strand"]]
        retained_a = []
        groups = set()
        processing = local[(local.strand == locus["C_strand"]) & ~local.is_C]
        b1 = processing[processing.query_id.map(lambda x: "B1" in supported_roles[x])]
        b2 = processing[processing.query_id.map(lambda x: "B2" in supported_roles[x])]
        for a in candidates_a.to_dict("records"):
            if any(x in gene_sequences[a["query_id"]] for x in "XBZJ"):
                continue
            a_groups = precursor_groups[a["subject_id"]] & {"canonical", "noncanonical"}
            if not a_groups or (a["strand"] != locus["C_strand"] and a_groups != {"canonical"}):
                continue
            if a["strand"] != locus["C_strand"] and not processing_pair(a, c, b1, b2, domains):
                continue
            terminal_overlap = terminal_precursor_overlap(a, c, locus, windows[locus["locus_id"]],
                gene_sequences[a["query_id"]], protocol.get("precursor_C_terminal_overlap_nt", []))
            if a["origin"] == "addressed_sORF_prediction" and any(max(a["start"], f["start"]) < min(a["end"], f["end"])
                       and not (f["is_C"] and terminal_overlap)
                       for f in local[local.aa_length > 200].to_dict("records")):
                continue
            if max(a["start"], locus["C_start"]) < min(a["end"], locus["C_end"]) and not terminal_overlap:
                continue
            retained_a.append(a)
        if retained_a:
            retained_a = (pd.DataFrame(retained_a).sort_values(["bitscore", "qcov", "query_id"], ascending=[False, False, True])
                          .drop_duplicates(["stop_coordinate", "strand"]).to_dict("records"))
            groups = set().union(*(precursor_groups[a["subject_id"]] & {"canonical", "noncanonical"} for a in retained_a))
        # A canonical processing pair must lie between a family precursor and
        # C. A single generic peptidase cannot count as both processing roles.
        canonical_pair = None
        for a in retained_a:
            canonical_pair = processing_pair(a, c, b1, b2, domains)
            if canonical_pair:
                break
        strong_c = locus["sequence_score"] >= protocol["strong_C_sequence_bits"] and locus["model_coverage"] >= protocol["strong_C_model_coverage"]
        status, family, reason = "homolog_only", "unresolved", "No admissible same-strand target-family precursor in addressed window"
        if protocol.get("canonical_opposite_A_with_processing_pair", False):
            reason = "No admissible target-family precursor in addressed window"
        if retained_a:
            status, reason = "partial_or_uncertain", "C/precursor evidence present; processing architecture or family relationship unresolved"
            if strong_c and locus["locus_id"] not in fragmented_c and locus["C_DNA_translation_match"] and not locus.get("C_partial_CDS", False) and not any(x in gene_sequences[query] for x in "XBZJ"):
                if groups == {"canonical"} and canonical_pair:
                    status, family, reason = "family_bgc_candidate", "canonical", "C and family precursor agree with local recognition/peptidase candidates; roles remain predictions"
                elif groups == {"noncanonical"}:
                    status, family, reason = "family_bgc_candidate", "noncanonical", "C and follower-precursor homology agree with compact noncanonical context; B1/B2 not required"
            if not strong_c:
                reason = "C/precursor evidence present, but C fails unchanged score/coverage requirement; completeness/function unresolved"
        if locus["locus_id"] in fragmented_c:
            status, reason = "partial_or_uncertain", "Adjacent CDSs align to consecutive portions of one C reference; intact C and function unresolved"
        if not locus["C_DNA_translation_match"]:
            status, reason = "partial_or_uncertain", "C annotation/assembly DNA translation conflict"
        elif locus.get("C_partial_CDS", False) or any(x in gene_sequences[query] for x in "XBZJ"):
            status, reason = "partial_or_uncertain", "Partial or ambiguous C annotation; complete sequence variant not established"
        chosen = []
        for a in retained_a:
            chosen.append((a["start"], a["end"], "A", a["query_id"], a["origin"]))
        chosen.append((int(c.start), int(c.end), "C", query, "annotated_CDS"))
        processing_roles = list(zip(("B1", "B2"), canonical_pair or ()))
        if canonical_pair and canonical_pair[0]["query_id"] == canonical_pair[1]["query_id"]:
            processing_roles = [("B1+B2", canonical_pair[0])]
        for role, row in processing_roles:
            if family == "canonical":
                chosen.append((row["start"], row["end"], role, row["query_id"], "annotated_CDS"))
        ordered = sorted(chosen, key=lambda x: x[0]*locus["C_strand"])
        core_start, core_end = min(x[0] for x in chosen), max(x[1] for x in chosen)
        window = windows[locus["locus_id"]]
        core = window.seq[core_start-locus["window_start"]:core_end-locus["window_start"]]
        if locus["C_strand"] == -1:
            core = core.reverse_complement()
        ambiguous_bases = sum(base not in "ACGT" for base in str(core).upper())
        if status == "family_bgc_candidate" and ambiguous_bases:
            status, reason = "partial_or_uncertain", "Coherent family architecture, but assembly ambiguity in defined DNA span"
        near = c_best.loc[query].to_dict() if query in c_best.index else {}
        review = reviewed_best.loc[query].to_dict() if query in reviewed_best.index else {}
        output.append({**locus, "status": status, "family": family, "interpretation_reason": reason,
                       "C_annotation_product": c["product"], "C_aa_length": c["aa_length"],
                       "C_query_id": query, "C_sha256": hashlib.sha256(gene_sequences[query].encode()).hexdigest(),
                       "C_nearest_reference": near.get("subject_id", ""), "C_reference_bitscore": near.get("bitscore", ""),
                       "C_reference_identity_percent": near.get("pident", ""), "C_reference_qcov": near.get("qcov", ""),
                       "C_reference_scov": near.get("scov", ""),
                       "C_reference_groups": (cyclase_reference.loc[near["subject_id"], "source_groups"]
                                              if near and near["subject_id"] in cyclase_reference.index else "canonical_LAR29" if near else ""),
                       "C_reviewed_reference": review.get("subject_id", ""), "C_reviewed_bitscore": review.get("bitscore", ""),
                       "C_strong": strong_c,
                       "C_integrity_evidence": "consecutive_C_segments" if locus["locus_id"] in fragmented_c else "no_split_flag_in_saved_alignments",
                       "C_domains": ";".join(sorted(domains[query])), "precursor_count": len(retained_a),
                       "precursor_C_strand_agreement": ";".join(str(a["strand"] == locus["C_strand"]) for a in retained_a),
                       "precursor_C_overlap_nt": ";".join(str(max(0,min(a["end"],int(c.end))-max(a["start"],int(c.start)))) for a in retained_a),
                       "precursor_query_ids": ";".join(a["query_id"] for a in retained_a),
                       "precursor_sha256": ";".join(hashlib.sha256(gene_sequences[a["query_id"]].encode()).hexdigest() for a in retained_a),
                       "precursor_origins": ";".join(a["origin"] for a in retained_a),
                       "architecture": "-".join(x[2] for x in ordered),
                       "chosen_core_gene_ids": ";".join(x[3] for x in ordered),
                       "core_start": core_start, "core_end": core_end,
                       "core_ambiguous_bases": ambiguous_bases,
                       "core_sha256": hashlib.sha256(str(core).upper().encode()).hexdigest(),
                       "edge_caveat": locus["left_edge_bp"] < 10000 or locus["right_edge_bp"] < 10000})
    pd.DataFrame(output).to_csv(tables / f"{output_name}_candidate_catalogue.tsv", sep="\t", index=False)
    pd.DataFrame(role_evidence).to_csv(tables / f"{output_name}_role_evidence.tsv", sep="\t", index=False)
    print(json.dumps({"contexts": len(output), "status_counts": pd.DataFrame(output).status.value_counts().to_dict()}))


if __name__ == "__main__":
    main()
