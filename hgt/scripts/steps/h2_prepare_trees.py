#!/usr/bin/env python3
"""Prepare separate component datasets and host concatenations with partitions."""
import argparse
from collections import defaultdict
import hashlib
import json
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_h2 import ROOT, OUT, read, write


def save_fasta(path, records):
    SeqIO.write((SeqRecord(Seq(s), id=key, description="") for key, s in sorted(records.items())), path, "fasta")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--components", action="store_true")
    parser.add_argument("--hosts", action="store_true")
    parser.add_argument("--actino", action="store_true")
    parser.add_argument("--fixed-host-models", action="store_true")
    args = parser.parse_args()
    destination = OUT / "actino" if args.actino else OUT
    directory = destination / "tree_inputs"
    directory.mkdir(exist_ok=True)
    jobs_path = destination / "tree_jobs.tsv"
    jobs = read(jobs_path) if jobs_path.exists() else []
    mag_hosts = set(json.loads((ROOT / "config/H2/H2_MAG_decision.json").read_text())["hosts"])
    qc = {r["host_id"]: r["host_tree_admission"] == "True" for r in read(destination / "marker_genome_QC.tsv")}
    equality = []
    if args.components:
        jobs = [r for r in jobs if r["kind"] != "component"]
        for path in sorted((OUT / "components").glob("*.faa")):
            family = path.stem
            records = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
            diagnostic = family.startswith("uncertain_")
            primary = {k: s for k, s in records.items() if k.split("_", 1)[0] in qc and (diagnostic or qc[k.split("_", 1)[0]])}
            if not primary:
                continue
            output = directory / (family + ".faa")
            save_fasta(output, primary)
            jobs.append({"kind": "component", "name": family, "input_path": str(output.relative_to(ROOT)), "partition_path": "", "filter": "gappyout", "dataset_role": "uncertain_Nature_diagnostic" if diagnostic else "primary_family_with_copy_number_caveats"})
            without = {k: s for k, s in primary.items() if k.split("_", 1)[0] not in mag_hosts}
            if without == primary:
                equality.append({"kind": "component", "dataset": family, "sensitivity": "without_two_MAGs", "status": "identical_inputs; both_MAGs_already_absent_or_QC_failed; no_duplicate_run"})
            else:
                output = directory / (family + "_noMAG.faa")
                save_fasta(output, without)
                jobs.append({"kind": "component", "name": family + "_noMAG", "input_path": str(output.relative_to(ROOT)), "partition_path": "", "filter": "gappyout", "dataset_role": "MAG_exclusion_sensitivity"})
        # Failed-marker MAGs are evaluated as an explicitly diagnostic addition.
        for family in ["canonical_C", "noncanonical_C"]:
            path = OUT / "components" / (family + ".faa")
            records = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
            augmented = {k: s for k, s in records.items() if qc.get(k.split("_", 1)[0], False) or (not args.actino and k.split("_", 1)[0] in mag_hosts)}
            primary = {k: s for k, s in records.items() if qc.get(k.split("_", 1)[0], False)}
            if augmented != primary:
                output = directory / (family + "_MAG_diagnostic.faa")
                save_fasta(output, augmented)
                jobs.append({"kind": "component", "name": family + "_MAG_diagnostic", "input_path": str(output.relative_to(ROOT)), "partition_path": "", "filter": "gappyout", "dataset_role": "below_threshold_MAGs_diagnostic_only"})
    if args.hosts:
        jobs = [r for r in jobs if r["kind"] != "host"]
        models = read(destination / "marker_models.tsv")
        model_choices = []
        for model in models:
            report = destination / "marker_trees" / model["marker_id"] / "ML_fastdiagnostic.iqtree"
            match = re.search(r"Best-fit model according to BIC:\s*(\S+)", report.read_text()) if report.exists() else None
            model_choices.append({"marker_id": model["marker_id"], "model": match.group(1) if match else "LG+G4",
                                  "basis": "BIC_from_independent_fast_marker_fit; parameters_reestimated_in_host_fit" if match else "fixed_LG_G4_fallback_for_uninformative_marker"})
        fixed_models = {r["marker_id"]: r["model"] for r in model_choices}
        if args.fixed_host_models:
            write(destination / "host_marker_model_choices.tsv", model_choices)
        alignments = {}
        for model in models:
            path = destination / "marker_trees" / model["marker_id"] / "trimmed_alignment.faa"
            if not path.exists():
                raise ValueError(f"Marker alignment incomplete: {model['marker_id']}")
            alignments[model["marker_id"]] = {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}
        panels = read(destination / "host_panels.tsv")
        for lineage in ["canonical", "noncanonical"]:
            primary = {r["host_id"] for r in panels if r["lineage"] == lineage and r["primary_admission"] == "True"}
            if not primary:
                continue
            requested_mags = {r["host_id"] for r in panels if r["lineage"] == lineage and r["MAG"] == "True"}
            sets = [(lineage, primary, "primary_independent_marker_tree")]
            if primary & mag_hosts:
                sets.append((lineage + "_noMAG", primary - mag_hosts, "MAG_exclusion_sensitivity"))
            else:
                equality.append({"kind": "host", "dataset": lineage, "sensitivity": "without_two_MAGs", "status": "identical_inputs; MAGs_not_admitted; no_duplicate_run"})
            if requested_mags - primary:
                sets.append((lineage + "_MAG_diagnostic", primary | requested_mags, "below_threshold_MAGs_diagnostic_only"))
            if args.actino and lineage == "noncanonical":
                uncertain = {r["host_id"] for r in read(OUT / "component_leaf_map.tsv") if r["family"] == "uncertain_Nature_system_C"}
                aliases = {r["genome_id"]: r["host_id"] for r in read(OUT / "host_aliases.tsv")}
                local = set(uncertain)
                from prepare_h2 import H1
                for relative in read(H1 / "H1_relative_panel.tsv"):
                    if any(aliases.get(g) in uncertain for g in relative["focal_genome_inputs"].split(";")):
                        h = aliases.get(relative["downloaded_genome_input"])
                        if h and qc.get(h, False):
                            local.add(h)
                sets.append(("uncertain_Nature", local, "four_uncertain_carriers_and_frozen_relatives; below_threshold_Kribbella_diagnostic_only"))
            for name, host_ids, role in sets:
                sequences = {host: "" for host in host_ids}
                partitions, position = [], 1
                for model in models:
                    seqs = alignments[model["marker_id"]]
                    length = len(next(iter(seqs.values())))
                    for host in host_ids:
                        sequences[host] += seqs.get(host, "-" * length)
                    partitions.append(f"  charset {model['marker_id']} = {position}-{position + length - 1};")
                    position += length
                output = directory / (name + "_hosts.faa")
                save_fasta(output, sequences)
                fixed = args.fixed_host_models and name != "uncertain_Nature"
                partition_path = directory / (name + ("_fixed_marker_partitions.nex" if fixed else "_partitions.nex"))
                assignments = ", ".join((fixed_models[model["marker_id"]] if fixed else "LG") + ":" + model["marker_id"] for model in models)
                partition_path.write_text("#nexus\nbegin sets;\n" + "\n".join(partitions) +
                                          "\n  charpartition H2 = " + assignments + ";\nend;\n")
                jobs.append({"kind": "host", "name": name, "input_path": str(output.relative_to(ROOT)), "partition_path": str(partition_path.relative_to(ROOT)), "filter": "pretrimmed_marker_partitions", "dataset_role": role})
    write(jobs_path, jobs)
    existing = read(destination / "MAG_sensitivity_identical_inputs.tsv") if (destination / "MAG_sensitivity_identical_inputs.tsv").exists() else []
    if equality:
        keys = {(r["kind"], r["dataset"]): r for r in existing + equality}
        write(destination / "MAG_sensitivity_identical_inputs.tsv", list(keys.values()))
    print(json.dumps({"component_jobs": sum(r["kind"] == "component" for r in jobs), "host_jobs": sum(r["kind"] == "host" for r in jobs), "identical_no_MAG_inputs": len(equality)}))


if __name__ == "__main__":
    main()
