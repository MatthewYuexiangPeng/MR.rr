#!/usr/bin/env bash
#SBATCH --job-name=mrrr_add_point_mrdag
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-160%20
#SBATCH --output=paper/output/additional_simulations/point_logs/mrdag_%A_%a.out
#SBATCH --error=paper/output/additional_simulations/point_logs/mrdag_%A_%a.err

set -euo pipefail

REPO_ROOT="${MRRR_REPO_ROOT:-${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}}"
R_MODULE="${MRRR_R_MODULE:-r/4.5.1-5zezbzn}"
MANIFEST="${MRRR_ADDITIONAL_POINT_MANIFEST:-paper/config/additional_point_tasks.csv}"
OVERWRITE="${MRRR_OVERWRITE:-false}"
DRY_RUN="${MRRR_DRY_RUN:-false}"

cd "${REPO_ROOT}"
mkdir -p paper/output/additional_simulations/point_chunks

if command -v module >/dev/null 2>&1; then
  module load "${R_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is unavailable. Set MRRR_R_MODULE to a valid R module." >&2
  exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Point-task manifest is missing: ${MANIFEST}" >&2
  exit 1
fi

echo "Repository root: ${REPO_ROOT}"
echo "Manifest: ${MANIFEST}"
echo "Slurm job: ${SLURM_ARRAY_JOB_ID:-unset}"
echo "MrDAG array index: ${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Git commit: $(git rev-parse HEAD)"
fi

Rscript --vanilla paper/scripts/24_dispatch_additional_point_task.R \
  --resource-class=mrdag \
  --array-index="${SLURM_ARRAY_TASK_ID}" \
  --manifest="${MANIFEST}" \
  --dry-run="${DRY_RUN}" \
  --overwrite="${OVERWRITE}"
