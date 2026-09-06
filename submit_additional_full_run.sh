#!/usr/bin/env bash
# Run with bash (do not source). Upload beside paper/ in the cluster checkout.
set -euo pipefail

export MRRR_REPO_ROOT="/home/peng.1276/MRrr-additional-cluster"
export MRRR_FULL_DIR="$MRRR_REPO_ROOT/paper/output/additional_simulations/cluster_full_b300_delta01"
cd "$MRRR_REPO_ROOT"

if [[ "${1:-}" == status ]]; then
  MRRR_IDS=$(cat "$MRRR_FULL_DIR/all_job_ids.txt")
  MRRR_START=$(cat "$MRRR_FULL_DIR/start_date.txt")
  sacct -X --array -n -P -S "$MRRR_START" -j "$MRRR_IDS" \
    --format=JobID,State,ExitCode | awk -F '|' '
      NF >= 3 && $1 != "" {
        count[$2]++; total++
        if ($2 !~ /^(COMPLETED|RUNNING|PENDING|CONFIGURING|COMPLETING)$/ ||
            ($2 == "COMPLETED" && $3 != "0:0"))
          print "CHECK:", $1, $2, $3
      }
      END {
        for (state in count) print state ":", count[state]
        print "Accounting records:", total, "(expected 842 when all tasks are recorded)"
      }'
  for MRRR_KIND in point bootstrap; do
    MRRR_LOG="$MRRR_FULL_DIR/logs/${MRRR_KIND}_merge.out"
    if [[ -f "$MRRR_LOG" ]]; then
      echo "Merge log: $MRRR_LOG"
      tail -n 16 "$MRRR_LOG"
    fi
    MRRR_ERR="$MRRR_FULL_DIR/logs/${MRRR_KIND}_merge.err"
    if [[ -s "$MRRR_ERR" ]]; then
      echo "Merge stderr: $MRRR_ERR"
      tail -n 12 "$MRRR_ERR"
    fi
  done
  exit 0
fi
if [[ $# -ne 0 ]]; then
  echo "Usage: bash submit_additional_full_run.sh [status]" >&2
  exit 1
fi

export MRRR_R_MODULE="r/4.5.1-5zezbzn"
module load "$MRRR_R_MODULE"
export R_LIBS="/home/peng.1276/R/mrrr-4.5:/home/peng.1276/R/library"
export MRRR_RSCRIPT="$(command -v Rscript)"
export MRRR_DRY_RUN=false MRRR_OVERWRITE=false
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

command -v sbatch >/dev/null
for MRRR_KIND in point bootstrap; do
  for MRRR_CLASS in standard mrdag; do
    test -s "paper/slurm/run_additional_${MRRR_KIND}_${MRRR_CLASS}.sbatch"
  done
done
test -s paper/scripts/25_merge_additional_point_chunks.R
test -s paper/scripts/29_merge_additional_bootstrap_chunks.R

# One fixed destination prevents accidental duplicate full submissions.
if [[ -e "$MRRR_FULL_DIR" ]]; then
  echo "Run directory already exists: $MRRR_FULL_DIR" >&2
  echo "No new jobs submitted. Use the status command; do not delete existing results." >&2
  exit 1
fi
mkdir -p "$(dirname "$MRRR_FULL_DIR")"
mkdir "$MRRR_FULL_DIR"
mkdir -p "$MRRR_FULL_DIR/logs" "$MRRR_FULL_DIR/provenance"
date '+%Y-%m-%d' > "$MRRR_FULL_DIR/start_date.txt"

# Retain the existing delta and all numeric configuration text. Lock only status.
"$MRRR_RSCRIPT" --vanilla - <<'RS'
stopifnot(getRversion() >= "4.5.0", getRversion() < "4.6.0")
packages <- c("MASS", "foreach", "doParallel", "mr.divw", "MrDAG", "CVXR", "osqp")
for (p in packages) {
  loadNamespace(p)
  cat(p, as.character(packageVersion(p)), find.package(p), "\n")
}
f <- "paper/config/approximate_low_rank_settings.csv"
d <- read.csv(f, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 4L, all(d$third_singular_value == 0.1),
          all(d$monte_carlo_replicates == 1000L), all(d$bootstrap_size == 300L),
          all(d$seed_base == 123L), all(d$status %in% c("provisional_delta", "locked")))
out <- file.path(Sys.getenv("MRRR_FULL_DIR"), "provenance")
stopifnot(file.copy(f, file.path(out, "approximate_config_before.csv")))
if (any(d$status != "locked")) {
  lines <- readLines(f, warn = FALSE)
  changed <- sub(",provisional_delta$", ",locked", lines)
  stopifnot(sum(changed != lines) == sum(d$status != "locked"))
  writeLines(changed, f, useBytes = TRUE)
}
stopifnot(all(read.csv(f, stringsAsFactors = FALSE)$status == "locked"))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
cat("Full-run configuration: delta=0.1, replicates=1000, B=300, status=locked\n")
RS

cp paper/config/approximate_low_rank_settings.csv "$MRRR_FULL_DIR/provenance/"
cp paper/config/rank_misspecification_settings.csv "$MRRR_FULL_DIR/provenance/"
git rev-parse HEAD > "$MRRR_FULL_DIR/provenance/git_commit.txt"
git diff HEAD -- paper/config/approximate_low_rank_settings.csv \
  > "$MRRR_FULL_DIR/provenance/configuration.patch"
sha256sum paper/scripts/2[0-9]_*.R paper/slurm/run_additional_*.sbatch \
  paper/config/approximate_low_rank_settings.csv paper/config/rank_misspecification_settings.csv \
  > "$MRRR_FULL_DIR/provenance/source_sha256.txt"

export MRRR_ADDITIONAL_POINT_MANIFEST="$MRRR_FULL_DIR/point_tasks.csv"
export MRRR_ADDITIONAL_BOOTSTRAP_MANIFEST="$MRRR_FULL_DIR/bootstrap_tasks.csv"

"$MRRR_RSCRIPT" --vanilla paper/scripts/23_make_additional_point_manifest.R \
  --total-replicates=1000 --standard-chunk-size=100 --mrdag-chunk-size=25 \
  --mrdag-niter=1000 --mrdag-burnin=200 --allow-provisional=false \
  --output="$MRRR_ADDITIONAL_POINT_MANIFEST" \
  --output-directory="$MRRR_FULL_DIR/point_chunks"

"$MRRR_RSCRIPT" --vanilla paper/scripts/27_make_additional_bootstrap_manifest.R \
  --total-replicates=1000 --standard-chunk-size=20 --mrdag-chunk-size=20 \
  --bootstrap-size=300 --allow-bootstrap-override=false \
  --standard-cores=16 --mrdag-cores=32 \
  --mrdag-niter=1000 --mrdag-burnin=200 --allow-provisional=false \
  --output="$MRRR_ADDITIONAL_BOOTSTRAP_MANIFEST" \
  --output-directory="$MRRR_FULL_DIR/bootstrap_chunks"

"$MRRR_RSCRIPT" --vanilla - <<'RS'
p <- read.csv(Sys.getenv("MRRR_ADDITIONAL_POINT_MANIFEST"))
b <- read.csv(Sys.getenv("MRRR_ADDITIONAL_BOOTSTRAP_MANIFEST"))
stopifnot(sum(p$resource_class == "standard") == 80L,
          sum(p$resource_class == "mrdag") == 160L,
          sum(b$resource_class == "standard") == 400L,
          sum(b$resource_class == "mrdag") == 200L,
          all(c(p$configuration_status, b$configuration_status) == "locked"),
          all(b$bootstrap_size == 300L), !any(b$bootstrap_size_override))
cat("Full manifests: 240 point chunks + 600 bootstrap chunks\n")
RS

record_job() {
  local label="$1" job="${2%%;*}"
  [[ "$job" =~ ^[0-9]+$ ]]
  printf '%s\n' "$job" > "$MRRR_FULL_DIR/${label}_job_id.txt"
  printf '%s\n' "$job" >> "$MRRR_FULL_DIR/job_ids_lines.txt"
  paste -sd, "$MRRR_FULL_DIR/job_ids_lines.txt" > "$MRRR_FULL_DIR/all_job_ids.txt"
  echo "Submitted $label: $job"
}

submit_array() {
  local kind="$1" class="$2" count="$3" concurrent="$4" cpus="$5" memory="$6"
  local job
  job=$(sbatch --parsable --export=ALL --chdir="$MRRR_REPO_ROOT" \
    --nodes=1 --ntasks=1 --array="1-${count}%${concurrent}" \
    --cpus-per-task="$cpus" --mem="$memory" --time=02:00:00 \
    --output="$MRRR_FULL_DIR/logs/${kind}_${class}_%A_%a.out" \
    --error="$MRRR_FULL_DIR/logs/${kind}_${class}_%A_%a.err" \
    "paper/slurm/run_additional_${kind}_${class}.sbatch")
  record_job "${kind}_${class}" "$job"
}

# Start the longest-running class first. Aggregate concurrency ceiling: 408 CPUs.
submit_array bootstrap mrdag 200 8 32 32G
submit_array bootstrap standard 400 8 16 16G
submit_array point standard 80 16 1 8G
submit_array point mrdag 160 8 1 8G

for MRRR_KIND in point bootstrap; do
  MRRR_STD=$(cat "$MRRR_FULL_DIR/${MRRR_KIND}_standard_job_id.txt")
  MRRR_DAG=$(cat "$MRRR_FULL_DIR/${MRRR_KIND}_mrdag_job_id.txt")
  if [[ "$MRRR_KIND" == point ]]; then
    MRRR_WRAP='exec "$MRRR_RSCRIPT" --vanilla paper/scripts/25_merge_additional_point_chunks.R --manifest="$MRRR_ADDITIONAL_POINT_MANIFEST" --output-dir="$MRRR_FULL_DIR/point_merged"'
  else
    MRRR_WRAP='exec "$MRRR_RSCRIPT" --vanilla paper/scripts/29_merge_additional_bootstrap_chunks.R --manifest="$MRRR_ADDITIONAL_BOOTSTRAP_MANIFEST" --output-dir="$MRRR_FULL_DIR/bootstrap_merged"'
  fi
  MRRR_JOB=$(sbatch --parsable --export=ALL --chdir="$MRRR_REPO_ROOT" \
    --job-name="mrrr_full_${MRRR_KIND}_merge" \
    --dependency="afterok:${MRRR_STD}:${MRRR_DAG}" \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --output="$MRRR_FULL_DIR/logs/${MRRR_KIND}_merge.out" \
    --error="$MRRR_FULL_DIR/logs/${MRRR_KIND}_merge.err" \
    --wrap="$MRRR_WRAP")
  record_job "${MRRR_KIND}_merge" "$MRRR_JOB"
done

touch "$MRRR_FULL_DIR/SUBMISSION_COMPLETE"
echo "FULL RUN SUBMITTED: 4 arrays + 2 dependent mergers."
echo "Results: $MRRR_FULL_DIR"
echo "Slurm now owns these jobs; you can close the terminal."
