#!/usr/bin/env bash
#SBATCH --job-name=mrrr_add_boot_std
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --array=1
#SBATCH --output=paper/output/additional_simulations/bootstrap_logs/standard_%A_%a.out
#SBATCH --error=paper/output/additional_simulations/bootstrap_logs/standard_%A_%a.err

# Submit from the repository root. Create bootstrap_logs BEFORE sbatch.
# Safe default: one array task. Explicitly set --array for the chosen manifest.
set -euo pipefail

MRRR_BOOT_ROOT="${MRRR_REPO_ROOT:-${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}}"
MRRR_BOOT_MODULE="${MRRR_R_MODULE:-r/4.5.1-5zezbzn}"
MRRR_BOOT_MANIFEST="${MRRR_ADDITIONAL_BOOTSTRAP_MANIFEST:-paper/config/additional_bootstrap_tasks.csv}"
MRRR_BOOT_OVERWRITE="${MRRR_OVERWRITE:-false}"
MRRR_BOOT_DRY_RUN="${MRRR_DRY_RUN:-false}"
MRRR_BOOT_INDEX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"

cd "${MRRR_BOOT_ROOT}"

if command -v module >/dev/null 2>&1; then
  module load "${MRRR_BOOT_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is unavailable. Set MRRR_R_MODULE to a valid R module." >&2
  exit 1
fi

if [[ ! -f "${MRRR_BOOT_MANIFEST}" ]]; then
  echo "Bootstrap-task manifest is missing: ${MRRR_BOOT_MANIFEST}" >&2
  exit 1
fi

# Parallelism is across bootstrap draws; prevent nested BLAS/OpenMP threads.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  Rscript --vanilla -e '
    a <- commandArgs(TRUE)
    m <- read.csv(a[1], stringsAsFactors = FALSE)
    stopifnot(all(c("resource_class", "array_index", "requested_cores") %in% names(m)))
    idx <- as.integer(a[3]); allocated <- as.integer(a[4])
    stopifnot(length(idx) == 1L, !is.na(idx), idx >= 1L,
              length(allocated) == 1L, !is.na(allocated), allocated >= 1L)
    r <- m[m$resource_class == a[2] & m$array_index == idx, , drop = FALSE]
    stopifnot(nrow(r) == 1L, !is.na(r$requested_cores),
              r$requested_cores >= 1L, r$requested_cores <= allocated)
    cat("Manifest cores / allocated CPUs:", r$requested_cores, "/", allocated, "\n")
  ' "${MRRR_BOOT_MANIFEST}" standard "${MRRR_BOOT_INDEX}" "${SLURM_CPUS_PER_TASK}"
fi

echo "Repository root: ${MRRR_BOOT_ROOT}"
echo "Manifest: ${MRRR_BOOT_MANIFEST}"
echo "Slurm job: ${SLURM_ARRAY_JOB_ID:-unset}"
echo "Standard bootstrap array index: ${MRRR_BOOT_INDEX}"

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Git commit: $(git rev-parse HEAD)"
fi

exec Rscript --vanilla paper/scripts/28_dispatch_additional_bootstrap_task.R \
  --resource-class=standard \
  --array-index="${MRRR_BOOT_INDEX}" \
  --manifest="${MRRR_BOOT_MANIFEST}" \
  --dry-run="${MRRR_BOOT_DRY_RUN}" \
  --overwrite="${MRRR_BOOT_OVERWRITE}"
