# Remaining paper computations. Values describe the archived numerical analyses.
spectral_remaining_config <- list(
  schema = "spectral-remaining-config-1", master_seed = 20260909L,
  rank_replicates = 1000L, rank_min = 1L, rank_alpha = 0.05,
  rank_designs = c("generic", "sparse_loading"), rank_chunk = 100L,
  # This diagnostic preserves signed eigenvalues as in rank_test_M0.
  # Undefined logarithms are recorded explicitly, never assigned a rank.
  rank_rule = "historical-signed-M0-with-explicit-undefined",
  pilot_replicates = 3L,
  simulation_eta_grid = sort(unique(c(10^seq(-5, -1, length.out = 30L), 1e-4, 1e-3))),
  simulation_eta_selected = c(1e-4, 1e-3, 1e-3, 1e-3),
  real_eta_grid = c(1e-4, 3e-4, 5e-4, 7e-4, 1e-3, 1.2e-3, 1.5e-3,
                    2e-3, 2.5e-3, 3e-3, 4e-3, 5e-3, 1e-2),
  real_ranks = c(1L, 2L), real_expected_selected_rank = 1L,
  real_eta_selected = c(`1` = 1.2e-3, `2` = 1e-3),
  real_bootstrap_size = 1000L, real_bootstrap_chunk = 20L,
  real_D_grid = 0:15, real_phi_c = 0.3,
  sparse_threshold = 1e-2, sparse_max_iter = 100L, sparse_tol = 1e-2,
  sparse_solver = "OSQP", tuning_zero_tol = 1e-2,
  mrdag_niter = 10000L, mrdag_burnin = 2000L,
  exposure_names = c("MMP12", "CNTN1", "FGL1", "MXRA8", "CNTFR", "SCG3", "HTRA1", "CLEC3B", "ANTXR2"),
  outcome_names = c("LAS", "CES", "SVS")
)
