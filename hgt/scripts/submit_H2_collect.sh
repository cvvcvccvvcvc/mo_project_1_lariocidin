#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=${HGT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
marker_search_job=${1:?Pass the completed or running H2 marker search job ID}
cd "$HGT_ROOT"
export HGT_ROOT
HGT_MARKER_JOB=$(sbatch --parsable --dependency="afterok:$marker_search_job" \
    --output=results/slurm/H2_marker_QC_%j.out --error=results/slurm/H2_marker_QC_%j.err \
    scripts/jobs/H2_marker_collect.sbatch)
HGT_COMPONENT_JOB=$(sbatch --parsable --dependency="afterok:$HGT_MARKER_JOB" \
    --output=results/slurm/H2_components_%j.out --error=results/slurm/H2_components_%j.err \
    scripts/jobs/H2_components.sbatch)
printf '%s %s\n' "$HGT_MARKER_JOB" "$HGT_COMPONENT_JOB"
