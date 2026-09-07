# Spectral simulation cluster workflow

This workflow runs the twelve generic, sparse-loading and approximate-low-rank
simulation cells, including working-rank misspecification. There are 1,000
replicates, 300 bootstrap draws, 100 unique result keys, 108 displayed table
rows and 1,800 computation tasks. Rank-two Reg./Sparse MR-rr rows are computed
once and shared between the main generic and rank-sensitivity table views.
Standalone rank selection, support recovery, tuning paths and real-data
inference are later milestones.

## 1. Validate and commit the worker source on Windows

Run from the repository root after applying `paper_spectral_workers_step3.patch`:

```bash
Rscript --vanilla paper/scripts/35_validate_spectral_workers.R --base-only --cores=2
```

This check uses native R and two PSOCK workers. It excludes Sparse and external
comparators and therefore reports `PARTIAL_PASS`; that is the expected result
for this command. Its ten checks include serial/parallel equality, real
spectral point/bootstrap estimates, checkpoint recovery, invalid output
rejection and matching rank-two table summaries. Reports are under
`paper/output/spectral_rebuild/worker_validation/`. A failing check exits
nonzero. A new execution replaces the report and retains a new fixture run.

If all dependencies are installed, omitting `--base-only` runs the full eleven
checks, including Sparse, IVW, SRIVW and MrDAG. The cluster submitter always
runs that full mode on a compute node before releasing production arrays.

Commit these paths explicitly after the Windows check:

```bash
git add .gitattributes .gitignore \
  paper/SPECTRAL_REBUILD.md paper/SPECTRAL_CLUSTER.md \
  paper/lib/paper_external.R paper/lib/paper_worker.R \
  paper/scripts/34_run_spectral_tasks.R \
  paper/scripts/35_validate_spectral_workers.R \
  paper/slurm/run_spectral_task.sbatch \
  paper/slurm/submit_spectral_simulations.sh &&
git diff --cached --check &&
git commit -m "Add spectral simulation workers and resumable cluster workflow" &&
git push
```

The local `collect_mrrr_audit.R` is unrelated to the computation and need not
be added. Keep historical results and generated output out of this commit.

## 2. Prepare a fixed cluster checkout and submit

For the initial setup on the current cluster, create a separate worktree:

```bash
cd /home/peng.1276/MRrr-additional-cluster &&
git fetch origin &&
git worktree add --detach /home/peng.1276/MRrr-spectral-cluster origin/paper/spectral-rebuild
```

Do not modify or pull into that worktree during this computation. Subsequent
source development can continue in the Windows checkout. The worktree pins
the calculation to the commit fetched above. If this worktree already exists,
inspect its commit and run status; do not overwrite it with this setup command.

```bash
cd /home/peng.1276/MRrr-spectral-cluster &&
bash paper/slurm/submit_spectral_simulations.sh submit
```

The submitter loads the previously working module `r/4.5.1-5zezbzn`, uses
`R_LIBS` if already exported or defaults to
`$HOME/R/mrrr-4.5:$HOME/R/library`, and limits BLAS/OpenMP to one thread per R
process. `MRRR_R_MODULE` can override the module. It requires Rscript, Git,
Slurm and `flock`. Use the existing cluster package installation; no
dependencies are installed or upgraded by these commands.

The production bundle is newly generated **on the cluster** under
`paper/output/spectral_rebuild/cluster_simulations_v1/`. The Windows design
bundle is not needed. `MRRR_SPECTRAL_RUN` can select a fresh absolute run path.
The full submitter requires 1,000 replicates, B=300 and the complete catalog.

Submission performs the following steps automatically:

1. Validate the prepared design and seal its sources, configuration and runtime.
2. Submit a two-CPU, 8-GB, 30-minute native preflight. It first verifies the
   production seal on the compute node, then executes the full worker check
   with actual Sparse, comparator and PSOCK computations.
3. Submit four arrays with `afterok` on preflight. A failed preflight prevents
   them from running. Dependent jobs use `--kill-on-invalid-dep=yes`.
4. Submit a one-CPU, 16-GB, 90-minute merger with `afterok` on all four arrays.
   Individual job IDs are recorded as soon as Slurm accepts each submission.

| Array | Tasks | CPUs/task | Memory/task | Max concurrent | Time limit/task |
| --- | ---: | ---: | ---: | ---: | --- |
| Bootstrap MrDAG | 600 | 32 | 32 GB | 8 | 4 hours |
| Bootstrap standard | 600 | 16 | 16 GB | 8 | 4 hours |
| Point standard | 120 | 1 | 8 GB | 16 | 4 hours |
| Point MrDAG | 480 | 1 | 8 GB | 8 | 4 hours |

This keeps the previously used 408-CPU aggregate concurrency ceiling. Actual
concurrency depends on available resources and account/QoS limits. The four
limits can be changed using `MRRR_DAG_CONCURRENT`, `MRRR_STD_CONCURRENT`,
`MRRR_POINT_STD_CONCURRENT`, and `MRRR_POINT_DAG_CONCURRENT`. They change
scheduling, not RNG streams or the number of fits. Time limits are allocations,
not runtime forecasts.

An initial complete submission prints six job IDs: preflight, four arrays and
merge. There are 1,802 expected top-level task records including preflight and
merge. Once submission completes, closing the terminal does not cancel these
Slurm jobs.

## 3. Inspect completion and recover interruptions

```bash
cd /home/peng.1276/MRrr-spectral-cluster &&
bash paper/slurm/submit_spectral_simulations.sh status
```

`status` shows the entire recorded Slurm history plus a fast chunk inventory.
`present` in this fast inventory means a file exists, not that its numerical
content has passed validation. Only the strict merger writes
`merged/STATUS.txt` with `SPECTRAL SIMULATION MERGE: PASS`.

Logs and validation reports are under the run's `logs/` and `preflight/`.
After a failure, inspect that job's `.err` and `.out`. A failed replicate also
retains `checkpoints/task-ID/failed_fits.csv` and RDS diagnostics. A cancelled
merger after an upstream failure is expected; its dependency is unsatisfied.
Historical failed/cancelled job IDs remain in accounting even after recovery.

After all previously submitted jobs from this run have left the active queue,
recover interruptions or transient failures with:

```bash
bash paper/slurm/submit_spectral_simulations.sh resume
```

Resume validates every existing chunk, submits only missing tasks, reuses
completed replicate checkpoints, and schedules a fresh preflight and merger.
It also recovers a submission that stopped after only some job IDs were
accepted. Invalid existing chunks cause an explicit stop; they are not deleted.
Deterministic numerical errors are replayed with the same inputs and seeds;
repeated resume alone cannot repair them. Source/runtime changes require an
explicit new run rather than mixing implementations into existing chunks.

Once merge has succeeded, resume returns the retained merge status and submits
no new jobs. The final data are:

| Path inside the run | Contents |
| --- | --- |
| `sealed_run.rds` | Input/source/runtime identity and preparation provenance |
| `submission_commit.txt` | Git commit used by Slurm submission |
| `merged/spectral_simulation_results.rds` | Per-replicate point estimates, bootstrap summaries, diagnostics, truth and provenance |
| `merged/summary.csv` | One summary per unique result key, including diagnostic counts |
| `merged/entrywise.csv` | Metrics for each of the 27 effect entries |
| `merged/table_rows.csv` | The 108 displayed rows, including identical shared rank-two summaries |
| `merged/STATUS.txt` | Strict simulation merge result |

Preserve the whole run directory for audit and replay. LaTeX tables, final
figures and the submission release bundle will be generated in later steps
from these validated results and the remaining paper analyses.
