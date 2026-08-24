#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_bootstrap_one.R <setting_idx> <estimator_set>
#   setting_idx:
#     1 = (me=2.5, effect=0.25)
#     2 = (me=1.0, effect=0.25)
#     3 = (me=2.5, effect=1.00)
#     4 = (me=1.0, effect=1.00)
#   estimator_set: "no_mrdag" | "mrdag_only" | "all" | "mr_only"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_bootstrap_one.R <setting_idx> <estimator_set>")
setting_idx   <- as.integer(args[1])
estimator_set <- as.character(args[2])
stopifnot(setting_idx %in% 1:4)
stopifnot(estimator_set %in% c("no_mrdag", "mrdag_only", "all", "mr_only"))

settings <- list(
  list(me = "2.5", effect = "0.25"),
  list(me = "1",   effect = "0.25"),
  list(me = "2.5", effect = "1"),
  list(me = "1",   effect = "1")
)
me_w  <- settings[[setting_idx]]$me
eff_w <- settings[[setting_idx]]$effect

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
cat(sprintf(">>> Setting %d: me=%s effect=%s | %s | %d cores\n",
            setting_idx, me_w, eff_w, estimator_set, n_cores))

source("scripts/MR_rr_estimators.R")
source("scripts/bootstrap_core.R")


# Simulation 1
# sim_tag <- "regular_C"
# input_file <- "results/simulate_result_pred_260717_regularC.RData"

# Simulation 2
sim_tag <- "sparseB_C"
input_file <- "results/simulate_result_pred_260718_sparseC.RData"

load(input_file)


parameters_list <- simulate_result_prediction$parameters_list

set.seed(123 + setting_idx)

if (estimator_set == "mrdag_only") {
  result <- nonpara_bootstrap_parallel_mrdag(
    parameters_list = parameters_list,
    me_weight = me_w, effect_weight = eff_w,
    iteration = 1000, bootstrap_size = 300,
    n_cores = n_cores,
    mrdag_niter = 1000, mrdag_burnin = 200,
    save_estimates = TRUE)
} else {
  result <- nonpara_bootstrap_parallel(
    parameters_list = parameters_list,
    me_weight = me_w, effect_weight = eff_w,
    regularization_rate = NULL,
    bootstrap_size = 300, iteration = 1000,
    n_cores = n_cores,
    estimator_set = estimator_set,
    save_estimates = TRUE)
}

out_file <- sprintf(
  "results/CI_%s_idx%d_me%s_eff%s_%s.RData",
  estimator_set,
  setting_idx,
  me_w,
  eff_w,
  sim_tag
)

save(result, file = out_file)
cat(sprintf(">>> Saved to %s\n", out_file))