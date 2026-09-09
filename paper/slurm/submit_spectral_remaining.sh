#!/usr/bin/env bash
# Public actions: submit | status | resume. Other actions run only in Slurm.
set -euo pipefail
MRRR_ACTION="${1:-}"
[[ $# -eq 1 && "$MRRR_ACTION" =~ ^(submit|status|resume|controller|controller_resume|task|merge)$ ]] || {
  echo 'Usage: bash paper/slurm/submit_spectral_remaining.sh submit|status|resume' >&2; exit 2;
}
if [[ "$MRRR_ACTION" =~ ^(submit|status|resume)$ ]]; then
  MRRR_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  export MRRR_REMAIN_ROOT=$(cd -- "$MRRR_SCRIPT_DIR/../.." && pwd)
else
  # Slurm runs a spool copy: never infer the repository from the spool filename.
  : "${SLURM_JOB_ID:?Batch action requires Slurm}"
  : "${MRRR_REMAIN_ROOT:?Repository root was not exported by the submitter}"
fi
cd "$MRRR_REMAIN_ROOT"
export MRRR_REMAIN_RUN="${MRRR_REMAIN_RUN:-$MRRR_REMAIN_ROOT/paper/output/spectral_rebuild/cluster_remaining_v2}"
[[ "$MRRR_REMAIN_RUN" == /* ]] || { echo 'MRRR_REMAIN_RUN must be an absolute path.' >&2; exit 2; }
export MRRR_REMAIN_REFERENCE="${MRRR_REMAIN_REFERENCE:-$MRRR_REMAIN_ROOT/paper/output/spectral_rebuild/cluster_simulations_v1_blas_recovery/run/merged/spectral_simulation_results.rds}"
export MRRR_REMAIN_REFERENCE_BUNDLE="${MRRR_REMAIN_REFERENCE_BUNDLE:-$(dirname -- "$(dirname -- "$MRRR_REMAIN_REFERENCE")")/simulation_bundle.rds}"
MRRR_SCRIPT="$MRRR_REMAIN_ROOT/paper/slurm/submit_spectral_remaining.sh"
record_job() {
  local label="$1" job="${2%%;*}"
  [[ "$job" =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $2" >&2; exit 1; }
  printf '%s\t%s\t%s\n' "$MRRR_BATCH" "$label" "$job" >> "$MRRR_REMAIN_RUN/job_history.tsv"
  printf '%s\n' "$job" >> "$MRRR_REMAIN_RUN/job_ids_lines.txt"
  paste -sd, "$MRRR_REMAIN_RUN/job_ids_lines.txt" > "$MRRR_REMAIN_RUN/all_job_ids.txt"
  MRRR_LAST_JOB="$job"
  echo "Submitted $label: $job"
}
MRRR_SITE=()
[[ -z "${MRRR_REMAIN_PARTITION:-}" ]] || MRRR_SITE+=(--partition="$MRRR_REMAIN_PARTITION")
[[ -z "${MRRR_REMAIN_ACCOUNT:-}" ]] || MRRR_SITE+=(--account="$MRRR_REMAIN_ACCOUNT")
[[ -z "${MRRR_REMAIN_CONSTRAINT:-}" ]] || MRRR_SITE+=(--constraint="$MRRR_REMAIN_CONSTRAINT")
[[ -z "${MRRR_REMAIN_NODE:-}" ]] || MRRR_SITE+=(--nodelist="$MRRR_REMAIN_NODE")

if [[ "$MRRR_ACTION" == status ]]; then
  if [[ -s "$MRRR_REMAIN_RUN/all_job_ids.txt" ]]; then
    MRRR_IDS=$(cat "$MRRR_REMAIN_RUN/all_job_ids.txt")
    sacct -X --array -n -P -j "$MRRR_IDS" --format=State | sort | uniq -c
    sacct -X --array -n -P -j "$MRRR_IDS" --format=JobID%28,State%24,ExitCode,Elapsed |
      awk -F '|' '$2 !~ /^COMPLETED/ {print}'
  fi
  if [[ -f "$MRRR_REMAIN_RUN/merged/STATUS.txt" ]]; then cat "$MRRR_REMAIN_RUN/merged/STATUS.txt"; fi
  if [[ -s "$MRRR_REMAIN_RUN/controller_job_id.txt" ]]; then
    MRRR_CID=$(cat "$MRRR_REMAIN_RUN/controller_job_id.txt")
    for MRRR_LOG in "$MRRR_REMAIN_RUN/logs/controller_${MRRR_CID}.out" "$MRRR_REMAIN_RUN/logs/controller_${MRRR_CID}.err"; do
      if [[ -f "$MRRR_LOG" ]]; then echo "$MRRR_LOG"; tail -n 12 "$MRRR_LOG"; fi
    done
  fi
  echo "Run directory: $MRRR_REMAIN_RUN"
  exit 0
fi

if [[ "$MRRR_ACTION" == submit || "$MRRR_ACTION" == resume ]]; then
  command -v sbatch >/dev/null
  command -v flock >/dev/null
  mkdir -p "$(dirname -- "$MRRR_REMAIN_RUN")"
  exec 9>"${MRRR_REMAIN_RUN}.submission.lock"
  flock --nonblock 9 || { echo 'Another submission is active.' >&2; exit 1; }
  git diff --quiet; git diff --cached --quiet
  MRRR_COMMIT=$(git rev-parse HEAD)
  test -s "$MRRR_REMAIN_REFERENCE" || { echo "Recovered main result not found: $MRRR_REMAIN_REFERENCE" >&2; exit 1; }
  test -s "$MRRR_REMAIN_REFERENCE_BUNDLE" || { echo "Original main simulation bundle not found: $MRRR_REMAIN_REFERENCE_BUNDLE" >&2; exit 1; }
  if [[ "$MRRR_ACTION" == submit ]]; then
    [[ ! -e "$MRRR_REMAIN_RUN" ]] || { echo 'Run exists; use status or resume.' >&2; exit 1; }
    mkdir -p "$MRRR_REMAIN_RUN/logs"
    printf '%s\n' "$MRRR_COMMIT" > "$MRRR_REMAIN_RUN/submission_commit.txt"
  else
    test -f "$MRRR_REMAIN_RUN/submission_commit.txt"
    [[ "$(cat "$MRRR_REMAIN_RUN/submission_commit.txt")" == "$MRRR_COMMIT" ]] || {
      echo 'Resume requires the original computation commit.' >&2; exit 1;
    }
    if [[ -f "$MRRR_REMAIN_RUN/merged/STATUS.txt" ]]; then cat "$MRRR_REMAIN_RUN/merged/STATUS.txt"; exit 0; fi
    if [[ -s "$MRRR_REMAIN_RUN/job_ids_lines.txt" ]]; then
      MRRR_ACTIVE=$(squeue --noheader --user="$(id -un)" --format='%i' |
        awk 'NR==FNR {wanted[$1]=1; next} {split($1,id,"_"); if (id[1] in wanted) print $1}' "$MRRR_REMAIN_RUN/job_ids_lines.txt" -)
      [[ -z "$MRRR_ACTIVE" ]] || { echo "Existing jobs remain active: $MRRR_ACTIVE" >&2; exit 1; }
    fi
  fi
  export MRRR_BATCH=$(date -u +%Y%m%dT%H%M%S)_$$
  MRRR_NEXT=controller
  [[ "$MRRR_ACTION" != resume ]] || MRRR_NEXT=controller_resume
  MRRR_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_REMAIN_ROOT" "${MRRR_SITE[@]}" \
    --job-name=mrrr_remaining_prepare --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G --time=04:00:00 \
    --output="$MRRR_REMAIN_RUN/logs/controller_%j.out" --error="$MRRR_REMAIN_RUN/logs/controller_%j.err" \
    "$MRRR_SCRIPT" "$MRRR_NEXT")
  record_job "$MRRR_NEXT" "$MRRR_JOB"
  printf '%s\n' "$MRRR_LAST_JOB" > "$MRRR_REMAIN_RUN/controller_job_id.txt"
  echo "Remaining-paper workflow submitted: $MRRR_REMAIN_RUN"
  exit 0
fi

source "$MRRR_REMAIN_ROOT/paper/slurm/with_spectral_runtime.sh"
case "$MRRR_ACTION" in
  controller|controller_resume)
    # Serialize job-history updates with the login submitter and any later resume.
    exec 9>"${MRRR_REMAIN_RUN}.submission.lock"
    flock 9
    [[ "$(git rev-parse HEAD)" == "$(cat "$MRRR_REMAIN_RUN/submission_commit.txt")" ]]
    git diff --quiet; git diff --cached --quiet
    if [[ ! -f "$MRRR_REMAIN_RUN/seal.rds" ]]; then
      Rscript --vanilla paper/scripts/38_validate_spectral_remaining.R --cores=2 \
        --reference="$MRRR_REMAIN_REFERENCE" --reference-bundle="$MRRR_REMAIN_REFERENCE_BUNDLE" \
        --output="$MRRR_REMAIN_RUN/preflight_$MRRR_BATCH"
      Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=prepare --profile=full \
        --run="$MRRR_REMAIN_RUN" --reference="$MRRR_REMAIN_REFERENCE" \
        --reference-bundle="$MRRR_REMAIN_REFERENCE_BUNDLE"
    fi
    Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=verify --run="$MRRR_REMAIN_RUN"
    Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=inventory --run="$MRRR_REMAIN_RUN"
    MRRR_DEPS="$SLURM_JOB_ID"
    submit_group() {
      local group="$1" cpus="$2" memory="$3" time="$4" concurrent="$5" indices job
      indices=$(awk -F '\t' -v g="$group" '$1==g {print $2}' "$MRRR_REMAIN_RUN/pending_arrays.tsv")
      [[ -n "$indices" ]] || { echo "Missing task group: $group" >&2; exit 1; }
      [[ "$indices" != '-' ]] || return 0
      job=$(sbatch --parsable --export=ALL --chdir="$MRRR_REMAIN_ROOT" "${MRRR_SITE[@]}" \
        --job-name="mrrr_rem_$group" --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory" --time="$time" \
        --array="${indices}%${concurrent}" --dependency="afterok:$SLURM_JOB_ID" --kill-on-invalid-dep=yes \
        --output="$MRRR_REMAIN_RUN/logs/${group}_%A_%a.out" --error="$MRRR_REMAIN_RUN/logs/${group}_%A_%a.err" \
        "$MRRR_SCRIPT" task)
      record_job "$group" "$job"; MRRR_DEPS="$MRRR_DEPS:$MRRR_LAST_JOB"
    }
    submit_group rank 1 4G 00:30:00 "${MRRR_REMAIN_RANK_CONCURRENT:-8}"
    submit_group eta 1 8G 04:00:00 "${MRRR_REMAIN_ETA_CONCURRENT:-4}"
    submit_group bootstrap_standard 8 16G 04:00:00 "${MRRR_REMAIN_STD_CONCURRENT:-8}"
    submit_group bootstrap_mrdag 8 16G 04:00:00 "${MRRR_REMAIN_DAG_CONCURRENT:-8}"
    MRRR_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_REMAIN_ROOT" "${MRRR_SITE[@]}" \
      --job-name=mrrr_remaining_merge --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
      --dependency="afterok:$MRRR_DEPS" --kill-on-invalid-dep=yes \
      --output="$MRRR_REMAIN_RUN/logs/merge_%j.out" --error="$MRRR_REMAIN_RUN/logs/merge_%j.err" \
      "$MRRR_SCRIPT" merge)
    record_job merge "$MRRR_JOB"
    echo 'SPECTRAL REMAINING ARRAYS SUBMITTED'
    ;;
  task)
    MRRR_TASK="${SLURM_ARRAY_TASK_ID:?Missing global task ID}"
    [[ "$MRRR_TASK" =~ ^[1-9][0-9]*$ ]]
    mkdir -p "$MRRR_REMAIN_RUN/locks"
    exec flock --nonblock "$MRRR_REMAIN_RUN/locks/task_$MRRR_TASK.lock" \
      Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=task --run="$MRRR_REMAIN_RUN" \
      --task="$MRRR_TASK" --cores="${SLURM_CPUS_PER_TASK:-1}"
    ;;
  merge)
    exec 9>"${MRRR_REMAIN_RUN}.merge.lock"; flock --nonblock 9
    exec Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=merge --run="$MRRR_REMAIN_RUN"
    ;;
esac
