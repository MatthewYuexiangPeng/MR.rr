#!/usr/bin/env bash
#SBATCH --job-name=mrrr_spectral
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
# Resource requests and log destinations are supplied by the submitter.
set -euo pipefail

cd "${MRRR_SPECTRAL_ROOT:?MRRR_SPECTRAL_ROOT is unset}"
if command -v module >/dev/null 2>&1; then
  module load "${MRRR_R_MODULE:?MRRR_R_MODULE is unset}"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MRRR_ACTION="${1:-task}"
case "$MRRR_ACTION" in
  preflight)
    Rscript --vanilla paper/scripts/34_run_spectral_tasks.R \
      --action=verify --run-dir="${MRRR_SPECTRAL_RUN:?}"
    exec Rscript --vanilla paper/scripts/35_validate_spectral_workers.R \
      --cores=2 --output="${MRRR_SPECTRAL_RUN:?}/preflight"
    ;;
  merge)
    exec Rscript --vanilla paper/scripts/34_run_spectral_tasks.R \
      --action=merge --run-dir="${MRRR_SPECTRAL_RUN:?}"
    ;;
  task)
    MRRR_PHASE="${MRRR_SPECTRAL_PHASE:?MRRR_SPECTRAL_PHASE is unset}"
    MRRR_CLASS="${MRRR_SPECTRAL_CLASS:?MRRR_SPECTRAL_CLASS is unset}"
    MRRR_INDEX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset}"
    [[ "$MRRR_PHASE" == point || "$MRRR_PHASE" == bootstrap ]]
    [[ "$MRRR_CLASS" == standard || "$MRRR_CLASS" == mrdag ]]
    [[ "$MRRR_INDEX" =~ ^[1-9][0-9]*$ ]]
    mkdir -p "$MRRR_SPECTRAL_RUN/locks"
    # The OS releases this lock on exit, cancellation and node/process failure.
    exec flock --nonblock \
      "$MRRR_SPECTRAL_RUN/locks/${MRRR_PHASE}_${MRRR_CLASS}_${MRRR_INDEX}.lock" \
      Rscript --vanilla paper/scripts/34_run_spectral_tasks.R --action=task \
      --run-dir="$MRRR_SPECTRAL_RUN" --phase="$MRRR_PHASE" \
      --resource-class="$MRRR_CLASS" --array-index="$MRRR_INDEX" \
      --cores="${SLURM_CPUS_PER_TASK:-1}"
    ;;
  *) echo "Unknown worker action: $MRRR_ACTION" >&2; exit 2 ;;
esac
