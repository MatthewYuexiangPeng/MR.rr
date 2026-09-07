#!/usr/bin/env bash
# Usage: bash paper/slurm/submit_spectral_simulations.sh submit|status|resume
# A new run prepares its configuration locally on the cluster, checks the native
# worker in a Slurm job, then releases four arrays and a dependent merger.
set -euo pipefail

MRRR_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export MRRR_SPECTRAL_ROOT=$(cd -- "$MRRR_SCRIPT_DIR/../.." && pwd)
cd "$MRRR_SPECTRAL_ROOT"
export MRRR_SPECTRAL_RUN="${MRRR_SPECTRAL_RUN:-$MRRR_SPECTRAL_ROOT/paper/output/spectral_rebuild/cluster_simulations_v1}"
[[ "$MRRR_SPECTRAL_RUN" == /* ]] || { echo 'MRRR_SPECTRAL_RUN must be absolute.' >&2; exit 2; }
MRRR_ACTION="${1:-}"
[[ $# -eq 1 && "$MRRR_ACTION" =~ ^(submit|status|resume)$ ]] || {
  echo 'Usage: bash paper/slurm/submit_spectral_simulations.sh submit|status|resume' >&2; exit 2;
}
export MRRR_R_MODULE="${MRRR_R_MODULE:-r/4.5.1-5zezbzn}"
if command -v module >/dev/null 2>&1; then module load "$MRRR_R_MODULE"; fi
export R_LIBS="${R_LIBS:-$HOME/R/mrrr-4.5:$HOME/R/library}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
command -v Rscript >/dev/null

if [[ "$MRRR_ACTION" == status ]]; then
  test -f "$MRRR_SPECTRAL_RUN/sealed_run.rds"
  if [[ -s "$MRRR_SPECTRAL_RUN/all_job_ids.txt" ]]; then
    MRRR_IDS=$(cat "$MRRR_SPECTRAL_RUN/all_job_ids.txt")
    echo 'Slurm history, including earlier attempts:'
    sacct -X --array -j "$MRRR_IDS" --format=JobID%28,State,ExitCode,Elapsed
  fi
  Rscript --vanilla paper/scripts/34_run_spectral_tasks.R \
    --action=inventory --fast=true --run-dir="$MRRR_SPECTRAL_RUN"
  if [[ -f "$MRRR_SPECTRAL_RUN/merged/STATUS.txt" ]]; then
    cat "$MRRR_SPECTRAL_RUN/merged/STATUS.txt"
  else
    echo "Merge status/logs: $MRRR_SPECTRAL_RUN/logs"
  fi
  exit 0
fi

command -v sbatch >/dev/null
command -v squeue >/dev/null
command -v flock >/dev/null
test -s paper/slurm/run_spectral_task.sbatch
mkdir -p "$(dirname -- "$MRRR_SPECTRAL_RUN")"
exec 9>"${MRRR_SPECTRAL_RUN}.submission.lock"
flock --nonblock 9 || { echo 'Another submission process is active.' >&2; exit 1; }
# Capture a fixed source checkout. Untracked local audit scripts are irrelevant.
git diff --quiet
git diff --cached --quiet
MRRR_COMMIT=$(git rev-parse HEAD)
if [[ "$MRRR_ACTION" == submit ]]; then
  [[ ! -e "$MRRR_SPECTRAL_RUN" ]] || {
    echo "Run exists: $MRRR_SPECTRAL_RUN. Use status or resume." >&2; exit 1;
  }
  Rscript --vanilla paper/scripts/32_prepare_spectral_simulations.R --output="$MRRR_SPECTRAL_RUN"
fi
test -d "$MRRR_SPECTRAL_RUN"
mkdir -p "$MRRR_SPECTRAL_RUN/logs"
if [[ -s "$MRRR_SPECTRAL_RUN/all_job_ids.txt" ]]; then
  # Query the current user's queue, then match recorded parent IDs. This also
  # works when completed jobs have aged out of the scheduler's active database.
  MRRR_ACTIVE=$(squeue --noheader --user="$(id -un)" --format='%i' |
    awk 'NR==FNR {wanted[$1]=1; next} {split($1,id,"_"); if (id[1] in wanted) print $1}' \
      "$MRRR_SPECTRAL_RUN/job_ids_lines.txt" -)
  [[ -z "$MRRR_ACTIVE" ]] || {
    echo "Existing jobs are still active: $MRRR_ACTIVE" >&2; exit 1;
  }
fi
Rscript --vanilla - "$MRRR_SPECTRAL_RUN" <<'RS'
b <- readRDS(file.path(commandArgs(TRUE)[1L], "simulation_bundle.rds"))
stopifnot(b$replicates == 1000L, b$bootstrap_size == 300L, nrow(b$catalog) == 100L)
RS
if [[ ! -f "$MRRR_SPECTRAL_RUN/sealed_run.rds" ]]; then
  Rscript --vanilla paper/scripts/34_run_spectral_tasks.R --action=seal --run-dir="$MRRR_SPECTRAL_RUN"
fi
Rscript --vanilla paper/scripts/34_run_spectral_tasks.R --action=verify --run-dir="$MRRR_SPECTRAL_RUN"
Rscript --vanilla paper/scripts/34_run_spectral_tasks.R --action=inventory --run-dir="$MRRR_SPECTRAL_RUN"
if [[ -f "$MRRR_SPECTRAL_RUN/merged/STATUS.txt" ]]; then
  cat "$MRRR_SPECTRAL_RUN/merged/STATUS.txt"
  exit 0
fi
if [[ -f "$MRRR_SPECTRAL_RUN/submission_commit.txt" ]]; then
  [[ "$(cat "$MRRR_SPECTRAL_RUN/submission_commit.txt")" == "$MRRR_COMMIT" ]] || {
    echo 'Checkout commit changed. Resume requires the original computation commit.' >&2; exit 1;
  }
else
  printf '%s\n' "$MRRR_COMMIT" > "$MRRR_SPECTRAL_RUN/submission_commit.txt"
fi
MRRR_BATCH=$(date -u +%Y%m%dT%H%M%S)_$$
MRRR_BATCH_IDS=()
record_job() {
  local label="$1" job="${2%%;*}"
  [[ "$job" =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $2" >&2; exit 1; }
  printf '%s\t%s\t%s\n' "$MRRR_BATCH" "$label" "$job" >> "$MRRR_SPECTRAL_RUN/job_history.tsv"
  printf '%s\n' "$job" >> "$MRRR_SPECTRAL_RUN/job_ids_lines.txt"
  paste -sd, "$MRRR_SPECTRAL_RUN/job_ids_lines.txt" > "$MRRR_SPECTRAL_RUN/all_job_ids.txt"
  MRRR_BATCH_IDS+=("$job")
  MRRR_LAST_JOB="$job"
  echo "Submitted $label: $job"
}
MRRR_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_SPECTRAL_ROOT" \
  --job-name=mrrr_spectral_check --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=00:30:00 \
  --output="$MRRR_SPECTRAL_RUN/logs/preflight_${MRRR_BATCH}_%j.out" \
  --error="$MRRR_SPECTRAL_RUN/logs/preflight_${MRRR_BATCH}_%j.err" \
  paper/slurm/run_spectral_task.sbatch preflight)
record_job preflight "$MRRR_JOB"
MRRR_PREFLIGHT="$MRRR_LAST_JOB"
submit_class() {
  local phase="$1" resource="$2" concurrent="$3" cpus="$4" memory="$5"
  local indices job
  indices=$(awk -F '\t' -v p="$phase" -v r="$resource" '$1==p && $2==r {print $3}' "$MRRR_SPECTRAL_RUN/pending_arrays.tsv")
  [[ -n "$indices" ]] || { echo 'Missing array specification.' >&2; exit 1; }
  [[ "$indices" != '-' ]] || return 0
  job=$(sbatch --parsable \
    --export="ALL,MRRR_SPECTRAL_PHASE=$phase,MRRR_SPECTRAL_CLASS=$resource" \
    --chdir="$MRRR_SPECTRAL_ROOT" --job-name="mrrr_${phase}_${resource}" \
    --dependency="afterok:$MRRR_PREFLIGHT" --kill-on-invalid-dep=yes \
    --nodes=1 --ntasks=1 --array="${indices}%${concurrent}" \
    --cpus-per-task="$cpus" --mem="$memory" --time=04:00:00 \
    --output="$MRRR_SPECTRAL_RUN/logs/${phase}_${resource}_%A_%a.out" \
    --error="$MRRR_SPECTRAL_RUN/logs/${phase}_${resource}_%A_%a.err" \
    paper/slurm/run_spectral_task.sbatch task)
  record_job "${phase}_${resource}" "$job"
}
# Same aggregate 408-CPU ceiling used by the preceding additional full run.
submit_class bootstrap mrdag "${MRRR_DAG_CONCURRENT:-8}" 32 32G
submit_class bootstrap standard "${MRRR_STD_CONCURRENT:-8}" 16 16G
submit_class point standard "${MRRR_POINT_STD_CONCURRENT:-16}" 1 8G
submit_class point mrdag "${MRRR_POINT_DAG_CONCURRENT:-8}" 1 8G
MRRR_DEPENDENCY=$(IFS=:; echo "${MRRR_BATCH_IDS[*]}")
MRRR_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_SPECTRAL_ROOT" \
  --job-name=mrrr_spectral_merge --dependency="afterok:$MRRR_DEPENDENCY" --kill-on-invalid-dep=yes \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:30:00 \
  --output="$MRRR_SPECTRAL_RUN/logs/merge_${MRRR_BATCH}_%j.out" \
  --error="$MRRR_SPECTRAL_RUN/logs/merge_${MRRR_BATCH}_%j.err" \
  paper/slurm/run_spectral_task.sbatch merge)
record_job merge "$MRRR_JOB"
printf '%s\n' "$MRRR_BATCH" > "$MRRR_SPECTRAL_RUN/SUBMISSION_COMPLETE"
echo "SPECTRAL SIMULATION JOBS SUBMITTED. Results: $MRRR_SPECTRAL_RUN"
