#!/usr/bin/env bash
#SBATCH --job-name=mrrr_sim_standard
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-160%20
#SBATCH --output=paper/output/full_run/logs/standard_%A_%a.out
#SBATCH --error=paper/output/full_run/logs/standard_%A_%a.err

set -euo pipefail

REPO_ROOT="${MRRR_REPO_ROOT:-${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}}"
R_MODULE="${MRRR_R_MODULE:-r/4.5.1-5zezbzn}"
MANIFEST="${MRRR_SIMULATION_MANIFEST:-paper/config/full_simulation_tasks.csv}"
OVERWRITE="${MRRR_OVERWRITE:-false}"

cd "${REPO_ROOT}"
mkdir -p paper/output/full_run/chunks

if command -v module >/dev/null 2>&1; then
  module load "${R_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is unavailable. Set MRRR_R_MODULE to a valid R module." >&2
  exit 1
fi

Rscript --vanilla paper/scripts/18_dispatch_simulation_task.R \
  --resource-class=standard \
  --array-index="${SLURM_ARRAY_TASK_ID}" \
  --manifest="${MANIFEST}" \
  --overwrite="${OVERWRITE}"
