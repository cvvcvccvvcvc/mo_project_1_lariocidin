#!/usr/bin/env python3
"""Search the published Actinobacteria panel in existing frozen native proteomes."""
import json
import os
import subprocess
from pathlib import Path

from prepare_h2 import ROOT, SNAP, OUT, read

hosts = [r for r in read(OUT / "marker_genome_QC.tsv") if r["phylum"] == "Actinomycetota"]
index = int(os.environ["SLURM_ARRAY_TASK_ID"])
for host in hosts[index::8]:
    directory = OUT / "native" / host["host_id"]
    command = [Path(os.environ["HGT_BIN"]) / "hmmsearch", "--cpu", "2", "--cut_ga",
               "--tblout", directory / "actino_markers.tbl", "--domtblout", directory / "actino_markers.domtbl",
               SNAP / "references/Actinobacteria.hmm", directory / "proteins.faa"]
    with (directory / "actino_markers.log").open("w") as handle:
        subprocess.run([str(x) for x in command], stdout=handle, stderr=subprocess.STDOUT, check=True)
print(json.dumps({"shard": index, "searched": len(hosts[index::8]), "panel": "published_Actinobacteria138"}))
