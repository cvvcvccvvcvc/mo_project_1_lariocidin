# Lariocidin-like BGCs: discovery and HGT analysis

This course research project consists of two independent analyses. [Diversity](diversity/README.md) searches for and curates lariocidin-like biosynthetic gene clusters (BGCs) in fixed public datasets. [HGT](hgt/README.md) tests phylogenetic hypotheses using its own input snapshot. The [research framework](docs/research.md) defines the research question and the limits of interpretation. A predicted BGC does not demonstrate product synthesis, and the current HGT evidence does not establish a transfer event.

## Reproducing the analysis

You need Linux, micromamba, access to public data sources and, for the full analysis, a Slurm cluster. Code and small frozen tables are tracked in Git. Genomes, environments, logs and large archives are excluded. Record the commit used for each run and keep its results separate from those produced with a different dataset version.

1. After cloning, enter `diversity/`. Create the environment with `micromamba create -y -p .env -f environment/environment.linux-64.lock` (on macOS, use `environment/environment.osx-arm64.lock`; the general environment specification is `environment/environment.yml`).
2. Run the reference analysis with `bash scripts/run_reference.sh`, then run `.env/bin/python -m unittest discover -s tests -v`. This workflow checks the installation, inputs and reference dataset. It does not replace the full search.
3. For the full search, obtain the two archives listed in the [external input registry](diversity/inputs/provenance/off_git_assets.tsv), place them in `diversity/inputs/archives/` and, from `diversity/`, verify them with `sha256sum -c inputs/provenance/off_git_assets.sha256`. The archives have not yet been published, so cloning the repository alone is insufficient to reproduce the full historical run.
4. On the cluster, from `diversity/`, submit `bash scripts/submit_baseline.sh`. After it completes successfully, submit `bash scripts/submit_species.sh`, `bash scripts/submit_additional.sh extension` and `bash scripts/submit_additional.sh controls`. Once all 17 shards are ready, submit `bash scripts/submit_diversity_review.sh`. These entry points configure the Slurm dependencies. The [diversity method](diversity/docs/method.md) documents the detailed sequence, inputs and outputs to verify.
5. Compare the final tables with the [recorded results](diversity/docs/results.md) and create a new manifest recording the commit, inputs, tool versions and SHA-256 checksums. [HGT](hgt/README.md) runs separately on its own snapshot. Do not silently replace that snapshot with a new diversity catalogue.

## Cluster access without personal details in Git

Obtain the login node address, username, key, any required jump host and working directory from the administrator. Keep these details outside the repository: a personal `mo-cluster` alias in `~/.ssh/config`, the private key in `~/.ssh/` and the clone path in `LARIOCIDIN_WORKDIR`. A typical session after setup is:

```sh
ssh mo-cluster
cd "$LARIOCIDIN_WORKDIR/diversity"
squeue -u "$USER"
bash scripts/submit_baseline.sh
```

Use the login node to prepare inputs and submit jobs. Run intensive computations through Slurm. Plan within a total budget of 64 CPUs, 512 GiB RAM and 250 GB of disk space, and check the current partitions, QOS and quota on the cluster itself. Do not commit hostnames, account names, keys, personal absolute paths or large datasets. Verify large inputs by SHA-256 checksum rather than filename.
