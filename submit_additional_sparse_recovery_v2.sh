#!/usr/bin/env bash
# Run with bash; do not source this file into an interactive shell.
set -euo pipefail

export MRRR_REPO_ROOT="/home/peng.1276/MRrr-additional-cluster"
export MRRR_FULL_DIR="$MRRR_REPO_ROOT/paper/output/additional_simulations/cluster_full_b300_delta01"
cd "$MRRR_REPO_ROOT"
MRRR_RECEIPT="$MRRR_FULL_DIR/sparse_recovery_v2_job_id.txt"
MRRR_START="$MRRR_FULL_DIR/sparse_recovery_v2_start_date.txt"
MRRR_SUCCESS="$MRRR_FULL_DIR/sparse_cholesky_recovery/SUCCESS.txt"

case "${1:-}" in
  status)
    test -s "$MRRR_RECEIPT"
    MRRR_RECOVERY_JOB=$(cat "$MRRR_RECEIPT")
    [[ "$MRRR_RECOVERY_JOB" =~ ^[0-9]+$ ]]
    MRRR_START_DATE=$(cat "$MRRR_START")
    sacct -S "$MRRR_START_DATE" -j "$MRRR_RECOVERY_JOB" \
      --format=JobID%24,State,ExitCode,Elapsed,MaxRSS
    for MRRR_SUFFIX in out err; do
      MRRR_LOG="$MRRR_FULL_DIR/logs/sparse_recovery_v2_${MRRR_RECOVERY_JOB}.${MRRR_SUFFIX}"
      if [[ -f "$MRRR_LOG" ]]; then
        echo "Log: $MRRR_LOG"
        tail -n 40 "$MRRR_LOG"
      fi
    done
    if [[ -s "$MRRR_SUCCESS" ]]; then
      cat "$MRRR_SUCCESS"
    else
      echo "Final recovery success marker is not present yet. Check this job's state and logs."
    fi
    exit 0
    ;;
  submit) ;;
  *) echo "Usage: bash submit_additional_sparse_recovery_v2.sh {submit|status}" >&2; exit 1 ;;
esac

if [[ -e "$MRRR_RECEIPT" ]]; then
  echo "Version 2 was already submitted. Use the status command; no duplicate job submitted." >&2
  cat "$MRRR_RECEIPT"
  exit 1
fi
if [[ -d "$MRRR_FULL_DIR/point_chunks_recovered" || -e "$MRRR_SUCCESS" || \
      -e "$MRRR_FULL_DIR/point_merged_recovered/additional_point_results.rds" ]]; then
  echo "Recovery outputs already exist. No files overwritten and no job submitted." >&2
  exit 1
fi

test -s recover_additional_sparse_cholesky.R
test -s "$MRRR_FULL_DIR/point_tasks.csv"
module load r/4.5.1-5zezbzn
export R_LIBS="/home/peng.1276/R/mrrr-4.5:/home/peng.1276/R/library"
export MRRR_RSCRIPT="$(command -v Rscript)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# Run the specific RNG regression before allocating a recovery job. This does
# not fit models, load archived results, or write simulation outputs.
"$MRRR_RSCRIPT" --vanilla - <<'RS'
core <- new.env(parent = globalenv())
sys.source("paper/scripts/20_additional_simulation_core.R", envir = core)
source("recover_additional_sparse_cholesky.R")
original <- core$additional_make_generic_truth
recovery_test_truth_rng(original, recovery_truth_function(original))
RS

mkdir -p "$MRRR_FULL_DIR/logs"
date '+%Y-%m-%d' > "$MRRR_START"
MRRR_RECOVERY_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_REPO_ROOT" \
  --job-name=mrrr_sparse_recovery_v2 \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 \
  --output="$MRRR_FULL_DIR/logs/sparse_recovery_v2_%j.out" \
  --error="$MRRR_FULL_DIR/logs/sparse_recovery_v2_%j.err" \
  --wrap='exec "$MRRR_RSCRIPT" --vanilla "$MRRR_REPO_ROOT/recover_additional_sparse_cholesky.R"')
MRRR_RECOVERY_JOB="${MRRR_RECOVERY_JOB%%;*}"
[[ "$MRRR_RECOVERY_JOB" =~ ^[0-9]+$ ]]
printf '%s\n' "$MRRR_RECOVERY_JOB" > "$MRRR_RECEIPT"
sha256sum recover_additional_sparse_cholesky.R submit_additional_sparse_recovery_v2.sh \
  > "$MRRR_FULL_DIR/sparse_recovery_v2_source_sha256.txt"
echo "Submitted sparse recovery v2: $MRRR_RECOVERY_JOB"
echo "Run: bash submit_additional_sparse_recovery_v2.sh status"
