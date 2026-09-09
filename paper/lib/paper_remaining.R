# Independent remaining-paper workflow; existing simulation seals stay untouched.
rem_assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
rem_md5 <- function(x) paper_object_md5(x)
rem_load <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  api <- new.env(parent = baseenv())
  for (f in c("paper_engine.R", "paper_simulation.R", "paper_external.R", "paper_worker.R", "paper_remaining.R"))
    sys.source(file.path(root, "paper/lib", f), envir = api)
  list(root = root, api = api, engine = api$paper_engine_load(root))
}
rem_config <- function(root, profile = "full") {
  rem_assert(profile %in% c("full", "development"), "Profile must be full or development.")
  e <- new.env(parent = baseenv()); sys.source(file.path(root, "paper/config/spectral_remaining.R"), e)
  c <- e$spectral_remaining_config
  rem_assert(identical(c$real_ranks, 1:2) && c$rank_min == 1L && c$rank_alpha == .05,
    "This manuscript stage requires real ranks 1/2 and archived rank search starting at one.")
  if (profile == "development") {
    c$rank_replicates <- 2L; c$rank_chunk <- 1L; c$pilot_replicates <- 1L
    c$simulation_eta_grid <- c(1e-4, 1e-3); c$real_eta_grid <- c(1e-3, 1.2e-3)
    c$real_bootstrap_size <- 3L; c$real_bootstrap_chunk <- 2L
    c$mrdag_niter <- 200L; c$mrdag_burnin <- 50L
  }
  c$profile <- profile; c
}
rem_catalog <- function() {
  x <- rbind(data.frame(method = c("ivw", "srivw", "mrdag"), rank = 0L),
    expand.grid(method = c("naive_mr_rr", "mr_rr", "regularized_mr_rr", "sparse_mr_rr"),
      rank = 1:2, stringsAsFactors = FALSE))
  x$key <- paste(x$method, x$rank, sep = "__")
  x$resource <- ifelse(x$method == "mrdag", "mrdag", "standard"); x
}
rem_sources <- function(ctx) {
  paths <- unique(c(ctx$engine$source_manifest$path, "paper/lib/paper_simulation.R",
    "paper/lib/paper_external.R", "paper/lib/paper_worker.R", "paper/lib/paper_remaining.R",
    "paper/config/spectral_simulations.csv", "paper/config/spectral_remaining.R",
    "paper/input/dat_1e-4.csv", "paper/input/rho_mat_1e-4.csv",
    "paper/scripts/37_spectral_remaining.R", "paper/scripts/38_validate_spectral_remaining.R",
    "paper/slurm/submit_spectral_remaining.sh", "paper/slurm/with_spectral_runtime.sh"))
  h <- unname(tools::md5sum(file.path(ctx$root, paths)))
  rem_assert(!anyNA(h), "A required remaining-paper source/input is missing.")
  data.frame(path = paths, md5 = h, stringsAsFactors = FALSE)
}
rem_runtime <- function(full = TRUE) {
  x <- paper_runtime(if (full) paper_sim_methods() else character())
  # This is a numerical guard in addition to the launcher's CPU/kernel policy.
  x$thread_environment <- Sys.getenv(c("OPENBLAS_CORETYPE", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"))
  x$numeric_probe <- paper_with_seed(271828L, {
    a <- matrix(stats::rnorm(177L * 9L), 177L, 9L)
    b <- matrix(stats::rnorm(9L * 3L), 9L, 3L)
    s <- crossprod(a) + diag(9L)
    rem_md5(list(a %*% b, s, chol(s), eigen(s, symmetric = TRUE)$values, solve(s)))
  })
  x
}
rem_rng <- function(cfg) {
  keys <- c("real/resample", paste0("real/method/", rem_catalog()$key),
    paste0("real/eta/", 1:2), unlist(lapply(1:4, function(s)
      paste0("pilot/", s, "/", seq_len(cfg$pilot_replicates)))))
  paper_with_seed(cfg$master_seed, {
    state <- get(".Random.seed", .GlobalEnv); out <- list()
    for (k in keys) {out[[k]] <- state; state <- parallel::nextRNGStream(state)}
    out
  })
}
rem_substate <- function(state, b) {
  if (b > 0L) for (j in seq_len(b)) state <- parallel::nextRNGSubStream(state)
  state
}
rem_covariances <- function(d) {
  n <- nrow(d$X)
  average <- function(se, cor) Reduce(`+`, lapply(seq_len(n), function(j) {
    v <- diag(se[j, ]); v %*% cor %*% v
  })) / n
  sx <- average(d$sx, d$cor_x); sy <- average(d$sy, d$cor_y)
  list(Sigma_X = sx, Sigma_Y = sy, W = solve(sy))
}
rem_real_data <- function(root, cfg) {
  x <- utils::read.csv(file.path(root, "paper/input/dat_1e-4.csv"), check.names = FALSE)
  rho <- as.matrix(utils::read.csv(file.path(root, "paper/input/rho_mat_1e-4.csv"), check.names = FALSE))
  rem_assert(nrow(x) == 177L && !anyDuplicated(x$SNP), "Expected 177 distinct harmonized SNPs.")
  scale <- sqrt(2 * x$ImpMAF * (1 - x$ImpMAF))
  d <- list(X = as.matrix(x[, paste0("gamma_exp", 1:9)]) * scale,
    Y = as.matrix(x[, paste0("gamma_out", 2:4)]) * scale,
    sx = as.matrix(x[, paste0("se_exp", 1:9)]) * scale,
    sy = as.matrix(x[, paste0("se_out", 2:4)]) * scale,
    cor_x = rho[1:9, 1:9], cor_y = rho[11:13, 11:13], snp = x$SNP,
    exposure_names = cfg$exposure_names, outcome_names = cfg$outcome_names)
  for (k in c("X", "Y", "sx", "sy", "cor_x", "cor_y"))
    rem_assert(all(is.finite(d[[k]])), paste("Non-finite real-data input:", k))
  rem_assert(all(d$sx > 0) && all(d$sy > 0), "SE values must be positive.")
  c(d, rem_covariances(d))
}
rem_resample_real <- function(d, idx) {
  out <- d
  for (k in c("X", "Y", "sx", "sy")) out[[k]] <- d[[k]][idx, , drop = FALSE]
  out$snp <- d$snp[idx]; cov <- rem_covariances(out)
  for (k in names(cov)) out[[k]] <- cov[[k]]
  out
}
rem_real_strength <- function(d, engine) {
  invroot <- solve(engine$estimators$.sqrt_matrix(d$cor_x))
  z <- t(vapply(seq_len(nrow(d$X)), function(i)
    as.vector(invroot %*% (d$X[i, ] / d$sx[i, ])), numeric(ncol(d$X))))
  M <- (crossprod(z) - nrow(z) * diag(ncol(z))) / sqrt(nrow(z))
  eigen((M + t(M)) / 2, symmetric = TRUE)$values
}
rem_rank <- function(d, engine, cfg) {
  n <- nrow(d$X); px <- ncol(d$X); py <- ncol(d$Y); full <- min(px, py)
  S <- crossprod(d$X) / n - d$Sigma_X; S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)$values
  result <- list(selected_rank = NA_integer_, archived_selected_rank = as.integer(full), valid = FALSE, error = "",
    min_exposure_eigenvalue = min(ev), indefinite = min(ev) <= 0,
    lambda = rep(NA_real_, full), statistic = rep(NA_real_, full - 1L),
    p_value = rep(NA_real_, full - 1L), df = (py - seq_len(full - 1L)) * (px - seq_len(full - 1L)))
  tryCatch({
    # Retain the archived signed inverse: no hidden PSD projection or regularization.
    Wroot <- engine$estimators$.sqrt_matrix(d$W)
    XY <- crossprod(d$X, d$Y) / n
    K <- Wroot %*% t(XY) %*% solve(S) %*% XY %*% Wroot
    lambda <- eigen((K + t(K)) / 2, symmetric = TRUE)$values
    result$lambda <- lambda
    rem_assert(all(is.finite(lambda)) && all(1 + lambda[-1L] > 0), "Rank statistic has an undefined log(1 + lambda).")
    result$statistic <- vapply(seq_len(full - 1L), function(r)
      (n - (px + py + 1) / 2) * sum(log1p(lambda[seq.int(r + 1L, full)])), numeric(1))
    result$p_value <- stats::pchisq(result$statistic, df = result$df, lower.tail = FALSE)
    rem_assert(all(is.finite(result$p_value)), "Non-finite rank p-value.")
    take <- which(result$p_value >= cfg$rank_alpha)
    result$selected_rank <- if (length(take)) as.integer(take[1L]) else as.integer(full)
    result$archived_selected_rank <- result$selected_rank
    result$valid <- TRUE; result
  }, error = function(e) {result$error <- conditionMessage(e); result})
}
rem_loss <- function(d, C, engine) {
  w <- engine$estimators$.sqrt_matrix(d$W)
  sum(((d$Y - d$X %*% t(C)) %*% w)^2) / nrow(d$X) - sum(diag(d$W %*% C %*% d$Sigma_X %*% t(C)))
}
rem_phi <- function(d, rank, engine, cfg, stage = "real_point") {
  S <- crossprod(d$X) / nrow(d$X) - d$Sigma_X
  invroot <- solve(engine$estimators$.sqrt_matrix(d$Sigma_X))
  A <- invroot %*% S %*% invroot
  mu <- min(eigen((A + t(A)) / 2, symmetric = TRUE)$values)
  phi <- mean(diag(d$Sigma_X))^2 / nrow(d$X) * exp(cfg$real_phi_c * (cfg$real_D_grid - sqrt(nrow(d$X)) * mu))
  rem_assert(all(is.finite(phi) & phi > 0), "Real-data phi grid is non-finite or underflowed.")
  fits <- lapply(phi, function(p) paper_engine_fit(engine, "regularized", stage, d$Y, d$X,
    rank, d$Sigma_X, d$W, regularization_rate = p))
  loss <- vapply(fits, function(f) rem_loss(d, f$AB, engine), numeric(1))
  rem_assert(all(is.finite(loss)), "Non-finite regularization objective.")
  i <- which.min(loss)
  list(fit = fits[[i]], selected_phi = phi[i], selected_D = cfg$real_D_grid[i],
    grid = data.frame(rank = rank, D = cfg$real_D_grid, phi = phi, objective = loss, selected = seq_along(phi) == i),
    homogeneous_scaled_siv = sqrt(nrow(d$X)) * mu)
}
rem_external_real <- function(method, d, cfg) {
  if (method == "mrdag") return(paper_external_fit(method, d$Y, d$X, d$Sigma_X, d$Sigma_Y, cfg))
  f <- getExportedValue("mr.divw", if (method == "ivw") "mvmr.ivw" else "mvmr.divw")
  C <- matrix(NA_real_, ncol(d$Y), ncol(d$X))
  for (j in seq_len(ncol(d$Y))) {
    # Preserve heteroskedastic SEs and their original trait/SNP order.
    args <- list(beta.exposure = d$X, se.exposure = d$sx,
      beta.outcome = as.vector(d$Y[, j]), se.outcome = d$sy[, j], gen_cor = d$cor_x)
    if (method == "srivw") args["phi_cand"] <- list(NULL)
    b <- do.call(f, args)$beta.hat
    rem_assert(length(b) == ncol(d$X), "Unexpected IVW/SRIVW coefficient dimensions.")
    C[j, ] <- b
  }
  list(AB = C)
}
rem_fit <- function(method, rank, d, ctx, cfg, state, point = NULL) {
  warnings <- character()
  stage <- if (is.null(point)) "real_point" else "real_bootstrap"
  ans <- paper_sim_with_state(state, withCallingHandlers({
    if (method %in% c("ivw", "srivw", "mrdag")) rem_external_real(method, d, cfg)
    else if (method == "regularized_mr_rr") {
      t <- rem_phi(d, rank, ctx$engine, cfg, stage); f <- t$fit; t$fit <- NULL; f$tuning <- t; f
    } else if (method == "sparse_mr_rr") {
      args <- list(engine = ctx$engine, method = "sparse", stage = if (is.null(point)) "real_point" else "real_bootstrap",
        Y = d$Y, X = d$X, rank = rank, Sigma_X = d$Sigma_X, W = d$W,
        sparse_max_iter = cfg$sparse_max_iter, sparse_tol = cfg$sparse_tol)
      if (is.null(point)) {
        args$sparse_lambda <- rep(cfg$real_eta_selected[as.character(rank)], 9L)
        args$sparse_threshold <- cfg$sparse_threshold; args$sparse_solver <- cfg$sparse_solver
      } else {
        args$sparse_support <- point$bootstrap_state$support
        # The frozen real-data bootstrap initializes from the original REFITTED A/B.
        args$sparse_A_init <- point$A; args$sparse_B_init <- point$B
      }
      do.call(paper_engine_fit, args)
    } else paper_engine_fit(ctx$engine, if (method == "mr_rr") "corrected" else "naive",
      stage, d$Y, d$X, rank, d$Sigma_X, d$W)
  }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}))
  rem_assert(is.matrix(ans$AB) && identical(dim(ans$AB), c(3L, 9L)) && all(is.finite(ans$AB)), "Invalid real-data C estimate.")
  if (method == "sparse_mr_rr" && !is.null(point))
    rem_assert(identical(ans$support, point$bootstrap_state$support) && all(ans$B[!ans$support] == 0), "Bootstrap support changed.")
  ans$recorded_warnings <- warnings; ans
}
rem_eta_path <- function(d, rank, grid, ctx, cfg, state) {
  rows <- fits <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    f <- paper_sim_with_state(state, paper_engine_fit(ctx$engine, "sparse", "tuning", d$Y, d$X,
      rank, d$Sigma_X, d$W, sparse_lambda = rep(grid[i], ncol(d$X)),
      sparse_threshold = cfg$sparse_threshold, sparse_max_iter = cfg$sparse_max_iter,
      sparse_tol = cfg$sparse_tol, sparse_solver = cfg$sparse_solver))
    rem_assert(!isTRUE(f$paper$sparse_refitted), "Tuning must use the penalized fit, without post-selection refit.")
    B <- f$B_raw; C <- f$AB_raw
    rem_assert(!is.null(B) && !is.null(C) && all(is.finite(B)) && all(is.finite(C)), "Missing unthresholded tuning coefficients.")
    rows[[i]] <- data.frame(rank = rank, eta = grid[i], objective = rem_loss(d, C, ctx$engine),
      zero_prop = mean(abs(B) < cfg$tuning_zero_tol), n_nonzero = sum(abs(B) >= cfg$tuning_zero_tol),
      converged = isTRUE(f$converged), projected = isTRUE(f$numerical_diagnostics$corrected_covariance_projected))
    fits[[i]] <- f
  }
  list(path = do.call(rbind, rows), fits = fits)
}
rem_task_plan <- function(cfg) {
  rows <- list(); add <- function(kind, design = "real", setting = 0L, first = 1L, last = first, resource = "standard") {
    i <- length(rows) + 1L
    rows[[i]] <<- data.frame(task_id = i, kind = kind, design = design, setting = setting,
      first = first, last = last, resource = resource, stringsAsFactors = FALSE)
  }
  for (d in cfg$rank_designs) for (s in 1:4) for (a in seq.int(1L, cfg$rank_replicates, cfg$rank_chunk))
    add("rank", d, s, a, min(a + cfg$rank_chunk - 1L, cfg$rank_replicates))
  for (s in 1:4) for (p in seq_len(cfg$pilot_replicates)) add("sim_eta", "generic", s, p)
  for (r in cfg$real_ranks) add("real_eta", first = r)
  for (res in c("standard", "mrdag")) for (a in seq.int(1L, cfg$real_bootstrap_size, cfg$real_bootstrap_chunk))
    add("real_bootstrap", first = a, last = min(a + cfg$real_bootstrap_chunk - 1L, cfg$real_bootstrap_size), resource = res)
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}
rem_input_path <- function(task) {
  if (task$kind == "rank") sprintf("inputs/rank_%s_s%d_%04d.rds", task$design, task$setting, task$first)
  else if (task$kind == "sim_eta") sprintf("inputs/pilot_s%d_p%d.rds", task$setting, task$first)
  else "inputs/real.rds"
}
rem_pilot_data <- function(bundle, setting, pilot, rng) {
  p <- bundle$parameters[[paste("generic", setting, sep = "/")]]
  paper_sim_with_state(rng[[paste("pilot", setting, pilot, sep = "/")]], {
    gauss <- function(S) matrix(stats::rnorm(177L * nrow(S)), 177L, nrow(S)) %*% chol(S)
    latent <- gauss(p$VX_tilde)
    list(X = latent + gauss(p$Sigma_X), Y = latent %*% t(p$C) + gauss(p$Sigma_Y),
      Sigma_X = p$Sigma_X, Sigma_Y = p$Sigma_Y, W = p$weight.matrix)
  })
}
rem_validate_reference_bundle <- function(ctx, ref, b, bundle_md5) {
  # The strict main merge authenticates the serialized design. Rebuilding its
  # SVD/calibration on a different BLAS is not an exact replay of that design.
  cells <- as.vector(outer(paper_sim_designs(), 1:4, paste, sep = "/"))
  rem_assert(identical(ref$schema, "spectral-merged-1") &&
    identical(ref$seal$schema, "spectral-run-1") && identical(ref$seal$validation_only, FALSE) &&
    identical(ref$run_id, ref$seal$run_id) && setequal(names(ref$data_md5), cells) &&
    all(lengths(ref$data_md5) == 1000L) && all(grepl("^[0-9a-f]{32}$", unlist(ref$data_md5))),
    "Reference is not a complete production strict merge.")
  rem_assert(identical(bundle_md5, ref$seal$bundle_md5),
    "Reference simulation_bundle.rds does not match the original seal; do not regenerate or replace it.")
  rem_assert(identical(b$schema, "spectral-simulation-design-1") &&
    identical(b$replicates, 1000L) && identical(b$bootstrap_size, 300L),
    "Reference design is not the required N=1000, B=300 main simulation.")
  for (key in c("config", "truths", "parameters", "prediction_exposure", "catalog", "table_map"))
    rem_assert(identical(b[[key]], ref[[key]]), paste("Reference design/merge disagree:", key))
  rem_assert(identical(b$config, paper_sim_read_config(ctx$root)),
    "Main simulation configuration differs from the reference design.")
  paths <- c(ctx$engine$source_manifest$path, "paper/config/spectral_simulations.csv",
    "paper/input/dat_1e-4.csv", "paper/input/rho_mat_1e-4.csv", "paper/lib/paper_simulation.R")
  # Script 32 appends its own fingerprint when it serializes a production bundle.
  recorded_paths <- b$source_manifest$path
  rem_assert(!anyDuplicated(recorded_paths) && all(paths %in% recorded_paths) &&
    identical(b$source_manifest$md5, unname(tools::md5sum(file.path(ctx$root, recorded_paths)))) &&
    identical(b$source_manifest$md5, ref$seal$sources$md5[match(recorded_paths, ref$seal$sources$path)]),
    "Main design source/input fingerprints differ from the original sealed design.")
  rem_assert(identical(b$rng$version, "cmrg-24-streams-per-replicate-v1") &&
    identical(b$rng$master_seed, b$config$master_seed[1L]) &&
    identical(b$rng$streams, paper_sim_make_registry(b$rng$master_seed)),
    "Reference design does not use the canonical main-simulation RNG registry.")
  invisible(TRUE)
}
rem_reference_inputs <- function(ctx, reference = "", reference_bundle = "") {
  if (!nzchar(reference)) {
    rem_assert(!nzchar(reference_bundle), "--reference-bundle requires --reference.")
    return(list(bundle = paper_sim_build_bundle(ctx$root, ctx$engine), ref = NULL,
      provenance = list(used = FALSE), files = character(), file_hashes = character()))
  }
  reference <- normalizePath(reference, winslash = "/", mustWork = TRUE)
  if (!nzchar(reference_bundle))
    reference_bundle <- file.path(dirname(dirname(reference)), "simulation_bundle.rds")
  rem_assert(file.exists(reference_bundle), paste("Original simulation_bundle.rds is required:",
    reference_bundle, "Supply --reference-bundle=PATH if stored elsewhere; no fresh-design fallback is allowed."))
  reference_bundle <- normalizePath(reference_bundle, winslash = "/", mustWork = TRUE)
  files <- c(reference = reference, bundle = reference_bundle)
  before <- tools::md5sum(files)
  ref <- readRDS(reference)
  # Object hashes serialize an R writer-version header as well as the data.
  # Keep the original exact policy; report a version mismatch explicitly.
  rem_assert(identical(ref$seal$runtime$R, R.version.string), paste(
    "Exact main-data replay requires the original R version:", ref$seal$runtime$R,
    "; current:", R.version.string))
  body <- ref$seal; body$run_id <- NULL
  rem_assert(identical(ref$seal$run_id, rem_md5(body)), "Original main-run seal checksum mismatch.")
  b <- readRDS(reference_bundle)
  rem_validate_reference_bundle(ctx, ref, b, unname(before[2L]))
  rem_assert(identical(before, tools::md5sum(files)), "Reference files changed while being read.")
  list(bundle = b, ref = ref, files = files, file_hashes = before,
    provenance = list(used = TRUE, run_id = ref$run_id, file_md5 = unname(before[1L]),
      bundle_file_md5 = unname(before[2L]),
      design_source = "original sealed simulation_bundle.rds; no recalibration",
      original_preparation = b$provenance))
}
rem_check_reference_data <- function(x, ref, design, setting, replicate) {
  actual <- rem_md5(list(X = x$X, Y = x$Y))
  if (!is.null(ref)) {
    expected <- unname(ref$data_md5[[paste(design, setting, sep = "/")]][replicate])
    rem_assert(identical(actual, expected), paste("Exact main-data replay failed:", design, setting, replicate,
      "expected", expected, "actual", actual,
      "The original design bundle has already been verified. Inspect R/BLAS/runtime diagnostics; no hash was bypassed."))
  }
  actual
}
rem_prepare <- function(ctx, run, profile = "full", reference = "", reference_bundle = "") {
  dir.create(run, recursive = TRUE, showWarnings = FALSE)
  run <- normalizePath(run, winslash = "/", mustWork = TRUE)
  rem_assert(!file.exists(file.path(run, "seal.rds")), "Run already prepared. Use resume, or a new run directory.")
  rem_assert(profile != "full" || nzchar(reference), "Full remaining runs require the strict main-simulation reference.")
  cfg <- rem_config(ctx$root, profile); runtime <- rem_runtime(TRUE); src <- rem_sources(ctx)
  main <- rem_reference_inputs(ctx, reference, reference_bundle)
  b <- main$bundle; ref <- main$ref
  if (!is.null(ref)) cat("Original sealed main design and RNG registry: PASS\nReference bundle MD5:",
    main$provenance$bundle_file_md5, "\n")
  task <- rem_task_plan(cfg); rng <- rem_rng(cfg)
  reference_hashes <- main$file_hashes
  checked <- 0L; write <- function(x, p) paper_write_rds(x, file.path(run, p))
  for (j in which(task$kind == "rank")) {
    t <- task[j, ]; p <- b$parameters[[paste(t$design, t$setting, sep = "/")]]
    data <- lapply(seq.int(t$first, t$last), function(i) {
      x <- paper_sim_data(b, t$design, t$setting, i)
      h <- rem_check_reference_data(x, ref, t$design, t$setting, i)
      if (!is.null(ref)) {
        checked <<- checked + 1L
      }
      list(X = x$X, Y = x$Y, data_md5 = h, replicate = i)
    })
    write(list(data = data, parameters = p), rem_input_path(t))
    if (j %% 10L == 0L) cat("Materialized rank tasks:", j, "/", sum(task$kind == "rank"), "\n")
  }
  for (j in which(task$kind == "sim_eta")) {
    t <- task[j, ]; d <- rem_pilot_data(b, t$setting, t$first, rng)
    write(d, rem_input_path(t))
  }
  real <- rem_real_data(ctx$root, cfg); rank <- rem_rank(real, ctx$engine, cfg)
  rem_assert(rank$valid && rank$selected_rank == cfg$real_expected_selected_rank,
    "Real-data selected rank differs from the expected primary rank one. Inspect before labeling manuscript results.")
  cat("Real-data working rank:", rank$selected_rank, "\n")
  catalog <- rem_catalog(); point <- list()
  for (i in seq_len(nrow(catalog))) {
    row <- catalog[i, ]; cat("Real point:", row$key, "\n")
    point[[row$key]] <- rem_fit(row$method, row$rank, real, ctx, cfg,
      rng[[paste0("real/method/", row$key)]])
  }
  indices <- vapply(seq_len(cfg$real_bootstrap_size), function(i)
    paper_sim_with_state(rem_substate(rng[["real/resample"]], i), sample.int(177L, 177L, replace = TRUE)), integer(177L))
  write(list(data = real, rank = rank, point = point, indices = indices,
    iv_strength = rem_real_strength(real, ctx$engine)), "inputs/real.rds")
  task$input_path <- vapply(seq_len(nrow(task)), function(i) rem_input_path(task[i, ]), character(1))
  task$input_md5 <- unname(tools::md5sum(file.path(run, task$input_path)))
  task$output_path <- sprintf("chunks/task_%04d.rds", task$task_id)
  files <- unique(task[, c("input_path", "input_md5")]); rownames(files) <- NULL
  write(task, "tasks.rds"); utils::write.csv(task, file.path(run, "tasks.csv"), row.names = FALSE)
  seal <- list(schema = "spectral-remaining-run-1", config = cfg, runtime = runtime, sources = src,
    tasks_md5 = unname(tools::md5sum(file.path(run, "tasks.rds"))), inputs = files, rng = rng,
    catalog = catalog, simulation_config = b$config, truths = b$truths,
    reference = c(main$provenance, list(matching_datasets = checked)),
    preparation_version = "remaining-cluster-v2-saved-design",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  rem_assert(identical(src, rem_sources(ctx)), "Computation source changed during preparation.")
  if (length(main$files)) rem_assert(identical(reference_hashes, tools::md5sum(main$files)),
    "Original reference files changed during preparation.")
  seal$run_id <- rem_md5(seal); write(seal, "seal.rds")
  writeLines(c("SPECTRAL REMAINING PREPARATION: PASS", paste("Run ID:", seal$run_id),
    paste("Profile:", profile), paste("Tasks:", nrow(task)),
    paste("Main-simulation datasets checked:", checked),
    if (length(main$files)) "Original design and merged reference files unchanged: PASS",
    "All rank/pilot datasets and real-data bootstrap indices are materialized once."), file.path(run, "PREPARATION.txt"))
  cat("SPECTRAL REMAINING PREPARATION: PASS\nRun ID:", seal$run_id, "\nTasks:", nrow(task), "\n")
  invisible(seal)
}
rem_open <- function(ctx, run, runtime = TRUE) {
  run <- normalizePath(run, winslash = "/", mustWork = TRUE)
  s <- readRDS(file.path(run, "seal.rds")); cp <- s; cp$run_id <- NULL
  rem_assert(identical(s$schema, "spectral-remaining-run-1") && identical(s$run_id, rem_md5(cp)), "Sealed metadata changed.")
  rem_assert(identical(s$sources, rem_sources(ctx)), "Source/input revision differs from the sealed run.")
  rem_assert(identical(s$tasks_md5, unname(tools::md5sum(file.path(run, "tasks.rds")))), "Task manifest changed.")
  if (runtime) rem_assert(identical(s$runtime, rem_runtime(TRUE)), "R/packages/BLAS numerical probe differ from the preparation environment.")
  t <- readRDS(file.path(run, "tasks.rds")); expected <- rem_task_plan(s$config)
  rem_assert(identical(t[, names(expected)], expected), "Missing, overlapping or altered task coverage.")
  list(run = run, seal = s, tasks = t)
}
rem_input <- function(state, task) {
  f <- file.path(state$run, task$input_path)
  rem_assert(identical(unname(tools::md5sum(f)), task$input_md5), paste("Materialized input changed:", task$input_path))
  readRDS(f)
}
rem_unit_key <- function(t, i) paste(t$kind, t$design, t$setting, t$resource, i, sep = "_")
rem_compute_unit <- function(ctx, state, task, input, id) {
  cfg <- state$seal$config; rng <- state$seal$rng
  if (task$kind == "rank") {
    d <- input$data[[id - task$first + 1L]]; p <- input$parameters
    d$Sigma_X <- p$Sigma_X; d$W <- p$weight.matrix
    result <- rem_rank(d, ctx$engine, cfg)
    return(list(diagnostic = result, data_md5 = d$data_md5))
  }
  if (task$kind == "sim_eta") return(rem_eta_path(input, 2L, cfg$simulation_eta_grid,
    ctx, cfg, rng[[paste("pilot", task$setting, id, sep = "/")]]))
  if (task$kind == "real_eta") return(rem_eta_path(input$data, id, cfg$real_eta_grid,
    ctx, cfg, rng[[paste0("real/eta/", id)]]))
  idx <- input$indices[, id]; d <- rem_resample_real(input$data, idx)
  catalog <- state$seal$catalog; catalog <- catalog[catalog$resource == task$resource, ]
  fits <- list()
  for (i in seq_len(nrow(catalog))) {
    c <- catalog[i, ]; fits[[c$key]] <- rem_fit(c$method, c$rank, d, ctx, cfg,
      rem_substate(rng[[paste0("real/method/", c$key)]], id), point = input$point[[c$key]])
  }
  list(fits = fits, indices_md5 = rem_md5(idx), data_md5 = rem_md5(list(X = d$X, Y = d$Y)),
    covariance_md5 = rem_md5(list(Sigma_X = d$Sigma_X, Sigma_Y = d$Sigma_Y, W = d$W)))
}
rem_unit <- function(ctx, state, task, input, id) {
  payload <- tryCatch(list(ok = TRUE, result = rem_compute_unit(ctx, state, task, input, id)),
    error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  r <- list(schema = "spectral-remaining-unit-1", run_id = state$seal$run_id,
    unit_key = rem_unit_key(task, id), task_id = task$task_id, unit_id = id,
    input_md5 = task$input_md5, payload = payload)
  r$body_md5 <- rem_md5(r); r
}
rem_validate_unit <- function(r, state, t, id) {
  cp <- r; cp$body_md5 <- NULL
  rem_assert(identical(r$schema, "spectral-remaining-unit-1") &&
    identical(r$run_id, state$seal$run_id) && identical(r$unit_key, rem_unit_key(t, id)) &&
    identical(r$task_id, t$task_id) && identical(r$unit_id, id) && identical(r$input_md5, t$input_md5) &&
    identical(r$body_md5, rem_md5(cp)), "Checkpoint identity or body hash mismatch.")
  rem_assert(isTRUE(r$payload$ok), paste("Work unit failed:", r$unit_key, r$payload$error))
  p <- r$payload$result
  if (t$kind == "rank") {
    rem_assert(is.list(p$diagnostic) && length(p$data_md5) == 1L, "Missing rank diagnostic.")
    if (p$diagnostic$valid) rem_assert(p$diagnostic$selected_rank %in% 1:3 && all(is.finite(p$diagnostic$p_value)), "Invalid rank output.")
    else rem_assert(nzchar(p$diagnostic$error), "Undefined rank needs an explicit reason.")
  } else if (t$kind %in% c("sim_eta", "real_eta")) {
    grid <- if (t$kind == "sim_eta") state$seal$config$simulation_eta_grid else state$seal$config$real_eta_grid
    rem_assert(identical(p$path$eta, grid) && all(is.finite(p$path$objective)) && all(is.finite(p$path$zero_prop)), "Incomplete tuning path.")
  } else {
    keys <- state$seal$catalog$key[state$seal$catalog$resource == t$resource]
    rem_assert(identical(names(p$fits), keys), "Missing or extra real-data fit.")
    for (f in p$fits) rem_assert(is.matrix(f$AB) && identical(dim(f$AB), c(3L, 9L)) && all(is.finite(f$AB)), "Non-finite bootstrap C.")
  }
  invisible(TRUE)
}
rem_validate_chunk <- function(x, state, task) {
  ids <- seq.int(task$first, task$last)
  rem_assert(identical(x$run_id, state$seal$run_id) && identical(x$task_id, task$task_id) &&
    identical(x$unit_ids, ids) && length(x$records) == length(ids), "Chunk coverage/identity mismatch.")
  for (j in seq_along(ids)) rem_validate_unit(x$records[[j]], state, task, ids[j])
  invisible(TRUE)
}
rem_run_task <- function(ctx, state, task_id, cores = 1L) {
  rem_assert(task_id %in% state$tasks$task_id, "Unknown task ID.")
  t <- state$tasks[state$tasks$task_id == task_id, ]; path <- file.path(state$run, t$output_path)
  input <- rem_input(state, t)
  if (file.exists(path)) {x <- readRDS(path); rem_validate_chunk(x, state, t); cat("SKIP valid task:", task_id, "\n"); return(invisible(x))}
  ids <- seq.int(t$first, t$last)
  checkpoints <- file.path(state$run, "checkpoints", paste0(vapply(ids, function(i) rem_unit_key(t, i), character(1)), ".rds"))
  records <- vector("list", length(ids)); missing <- integer()
  for (j in seq_along(ids)) {
    if (file.exists(checkpoints[j])) {
      r <- readRDS(checkpoints[j])
      if (isTRUE(r$payload$ok)) {rem_validate_unit(r, state, t, ids[j]); records[[j]] <- r} else missing <- c(missing, j)
    } else missing <- c(missing, j)
  }
  pool <- NULL
  if (cores > 1L && t$kind == "real_bootstrap" && length(missing) > 1L) {
    pool <- parallel::makePSOCKcluster(min(cores, length(missing)))
    on.exit(parallel::stopCluster(pool), add = TRUE)
    parallel::clusterCall(pool, function(root, run, task, data) {
      a <- new.env(parent = baseenv()); sys.source(file.path(root, "paper/lib/paper_remaining.R"), a)
      ctx <- a$rem_load(root); st <- ctx$api$rem_open(ctx, run)
      assign(".remaining_job", list(ctx = ctx, state = st, task = task, input = data), .GlobalEnv)
      TRUE
    }, ctx$root, state$run, t, input)
  }
  # Persist after each bounded wave; retries reuse completed individual draws.
  waves <- split(missing, ceiling(seq_along(missing) / max(1L, cores)))
  for (wave in waves) {
    values <- if (is.null(pool)) lapply(ids[wave], function(i) rem_unit(ctx, state, t, input, i)) else
      parallel::parLapply(pool, ids[wave], function(i) {
        z <- get(".remaining_job", .GlobalEnv); z$ctx$api$rem_unit(z$ctx, z$state, z$task, z$input, i)
      })
    for (k in seq_along(wave)) {
      j <- wave[k]; records[[j]] <- values[[k]]; paper_write_rds(values[[k]], checkpoints[j])
    }
    for (k in seq_along(wave)) {
      j <- wave[k]
      rem_validate_unit(values[[k]], state, t, ids[j])
      cat("Task", task_id, "unit", ids[j], "PASS\n")
    }
  }
  x <- list(run_id = state$seal$run_id, task_id = t$task_id, unit_ids = ids, records = records)
  rem_validate_chunk(x, state, t); paper_write_rds(x, path)
  cat("SPECTRAL REMAINING TASK: PASS", task_id, "\n"); invisible(x)
}
rem_inventory <- function(state, write = TRUE) {
  result <- state$tasks; result$status <- "missing"; result$detail <- ""
  for (i in seq_len(nrow(result))) {
    f <- file.path(state$run, result$output_path[i]); if (!file.exists(f)) next
    err <- tryCatch({rem_validate_chunk(readRDS(f), state, state$tasks[i, ]); NULL}, error = function(e) conditionMessage(e))
    result$status[i] <- if (is.null(err)) "valid" else "invalid"
    if (!is.null(err)) result$detail[i] <- err
  }
  if (write) utils::write.csv(result, file.path(state$run, "inventory.csv"), row.names = FALSE)
  result
}
rem_merge <- function(ctx, state) {
  inv <- rem_inventory(state); print(table(inv$kind, inv$resource, inv$status))
  rem_assert(all(inv$status == "valid"), "Merge refused: missing/invalid remaining-paper tasks.")
  for (i in seq_len(nrow(state$seal$inputs)))
    rem_assert(identical(unname(tools::md5sum(file.path(state$run, state$seal$inputs$input_path[i]))),
      state$seal$inputs$input_md5[i]), "Materialized input changed before merge.")
  target <- file.path(state$run, "merged")
  rem_assert(!dir.exists(target), "Merged output already exists; retain it and inspect STATUS.txt.")
  cfg <- state$seal$config; rank_rows <- eta_sim <- eta_real <- list()
  real <- readRDS(file.path(state$run, "inputs/real.rds")); boot <- list(standard = vector("list", cfg$real_bootstrap_size), mrdag = vector("list", cfg$real_bootstrap_size))
  for (i in seq_len(nrow(state$tasks))) {
    t <- state$tasks[i, ]; chunk <- readRDS(file.path(state$run, t$output_path))
    rem_validate_chunk(chunk, state, t)
    for (r in chunk$records) {
      p <- r$payload$result
      if (t$kind == "rank") {
        z <- p$diagnostic
        rank_rows[[length(rank_rows) + 1L]] <- data.frame(design = t$design, setting = t$setting,
          replicate = r$unit_id, rank = z$selected_rank, archived_rank = z$archived_selected_rank,
          valid = z$valid, indefinite = z$indefinite,
          min_exposure_eigenvalue = z$min_exposure_eigenvalue, p_rank1 = z$p_value[1L], p_rank2 = z$p_value[2L],
          error = z$error, data_md5 = p$data_md5, stringsAsFactors = FALSE)
      } else if (t$kind == "sim_eta") eta_sim[[length(eta_sim) + 1L]] <- cbind(setting = t$setting, pilot = r$unit_id, p$path)
      else if (t$kind == "real_eta") eta_real[[length(eta_real) + 1L]] <- p$path
      else boot[[t$resource]][[r$unit_id]] <- p
    }
  }
  ranks <- do.call(rbind, rank_rows); sim <- do.call(rbind, eta_sim); eta <- do.call(rbind, eta_real)
  rem_assert(nrow(ranks) == 8L * cfg$rank_replicates && !anyDuplicated(ranks[, c("design", "setting", "replicate")]), "Rank experiment incomplete.")
  for (b in seq_len(cfg$real_bootstrap_size)) {
    a <- boot$standard[[b]]; d <- boot$mrdag[[b]]
    rem_assert(!is.null(a) && !is.null(d) && identical(a$indices_md5, rem_md5(real$indices[, b])) &&
      identical(a$indices_md5, d$indices_md5) && identical(a$data_md5, d$data_md5) &&
      identical(a$covariance_md5, d$covariance_md5), "Real-data bootstrap inputs differ across resource classes.")
    original <- rem_resample_real(real$data, real$indices[, b])
    rem_assert(identical(a$data_md5, rem_md5(list(X = original$X, Y = original$Y))) &&
      identical(a$covariance_md5, rem_md5(list(Sigma_X = original$Sigma_X, Sigma_Y = original$Sigma_Y, W = original$W))),
      "Real-data bootstrap replay differs from the materialized original data.")
  }
  catalog <- state$seal$catalog; summaries <- arrays <- list()
  for (i in seq_len(nrow(catalog))) {
    row <- catalog[i, ]; fs <- lapply(boot[[row$resource]], function(b) b$fits[[row$key]])
    C <- vapply(fs, function(f) as.vector(f$AB), numeric(27L))
    rem_assert(identical(dim(C), c(27L, cfg$real_bootstrap_size)) && all(is.finite(C)), "Missing/non-finite bootstrap draws.")
    arrays[[row$key]] <- C
    q <- t(apply(C, 1L, stats::quantile, probs = c(.025, .975), names = FALSE, type = 7L))
    summaries[[row$key]] <- data.frame(key = row$key, method = row$method, rank = row$rank,
      outcome = rep(cfg$outcome_names, 9L), exposure = rep(cfg$exposure_names, each = 3L),
      estimate = as.vector(real$point[[row$key]]$AB), se = apply(C, 1L, stats::sd),
      lower = q[, 1L], upper = q[, 2L], excludes_zero = q[, 1L] > 0 | q[, 2L] < 0)
    if (row$method == "sparse_mr_rr") for (f in fs)
      rem_assert(identical(f$support, real$point[[row$key]]$bootstrap_state$support) && all(f$B[!f$support] == 0), "A sparse bootstrap fit changed support.")
  }
  rank_summary <- do.call(rbind, lapply(split(ranks, interaction(ranks$design, ranks$setting, drop = TRUE)), function(x) {
    data.frame(design = x$design[1], setting = x$setting[1], N = nrow(x),
      rank1_percent = 100 * sum(x$rank == 1L, na.rm = TRUE) / nrow(x),
      rank2_percent = 100 * sum(x$rank == 2L, na.rm = TRUE) / nrow(x),
      rank3_percent = 100 * sum(x$rank == 3L, na.rm = TRUE) / nrow(x),
      archived_rank1_percent = 100 * mean(x$archived_rank == 1L),
      archived_rank2_percent = 100 * mean(x$archived_rank == 2L),
      archived_rank3_percent = 100 * mean(x$archived_rank == 3L),
      undefined_percent = 100 * sum(!x$valid) / nrow(x), indefinite_percent = 100 * mean(x$indefinite))
  }))
  sim_summary <- stats::aggregate(sim[, c("objective", "zero_prop", "n_nonzero")], sim[, c("setting", "eta")], mean)
  real_summary <- do.call(rbind, summaries)
  phi_point <- do.call(rbind, lapply(real$point[c("regularized_mr_rr__1", "regularized_mr_rr__2")], function(f) f$tuning$grid))
  phi_boot <- do.call(rbind, lapply(seq_len(cfg$real_bootstrap_size), function(b) do.call(rbind, lapply(1:2, function(r) {
    f <- boot$standard[[b]]$fits[[paste0("regularized_mr_rr__", r)]]$tuning
    data.frame(bootstrap = b, rank = r, phi = f$selected_phi, D = f$selected_D, homogeneous_scaled_siv = f$homogeneous_scaled_siv)
  }))))
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(catalog)), function(i) {
    row <- catalog[i, ]; fs <- lapply(boot[[row$resource]], function(b) b$fits[[row$key]])
    data.frame(key = row$key, point_nonconverged = identical(real$point[[row$key]]$converged, FALSE),
      bootstrap_nonconverged = sum(vapply(fs, function(f) identical(f$converged, FALSE), logical(1))),
      bootstrap_projected = sum(vapply(fs, function(f) isTRUE(f$numerical_diagnostics$corrected_covariance_projected), logical(1))),
      bootstrap_warnings = sum(vapply(fs, function(f) length(f$recorded_warnings), integer(1))))
  }))
  output <- list(schema = "spectral-remaining-merged-1", run_id = state$seal$run_id,
    production = cfg$profile == "full", seal = state$seal, rank_replicates = ranks, rank_summary = rank_summary,
    simulation_eta = sim, simulation_eta_summary = sim_summary, real_eta = eta,
    real_data = real$data, real_rank = real$rank, real_point = real$point,
    real_iv_strength = real$iv_strength,
    real_bootstrap_indices = real$indices, real_bootstrap = arrays,
    real_summary = real_summary, real_phi_point = phi_point, real_phi_bootstrap = phi_boot,
    diagnostics = diagnostics)
  stage <- tempfile("remaining_merge_", tmpdir = state$run); dir.create(stage)
  paper_write_rds(output, file.path(stage, "spectral_remaining_results.rds"))
  for (name in c("rank_replicates", "rank_summary", "simulation_eta", "simulation_eta_summary", "real_eta",
    "real_summary", "real_phi_point", "real_phi_bootstrap", "diagnostics"))
    utils::write.csv(output[[name]], file.path(stage, paste0(name, ".csv")), row.names = FALSE)
  status <- c("SPECTRAL REMAINING STRICT MERGE: PASS", paste("Production:", output$production),
    paste("Run ID:", output$run_id), paste("Rank datasets:", nrow(ranks)),
    paste("Undefined rank diagnostics (explicitly retained):", sum(!ranks$valid)),
    paste("Real-data bootstrap draws per unique method/rank:", cfg$real_bootstrap_size),
    "Identical original data, SNP resamples and recalculated covariances across resource classes: PASS",
    "Real-data IVW/SRIVW/MrDAG fits are shared between the rank-one and rank-two analyses.",
    "All real-data bootstrap estimates finite; percentile intervals use all draws.",
    "Sparse support fixed; original-data refitted A/B used for bootstrap initialization.",
    "Regularized MR-rr uses the spectral solver and reselects phi in every real-data bootstrap.",
    "Archived rank fallback is also saved explicitly for comparison (undefined statistics formerly defaulted to full rank).",
    "Review rank undefined counts and convergence diagnostics before updating manuscript claims.")
  writeLines(status, file.path(stage, "STATUS.txt"))
  rem_assert(file.rename(stage, target), "Could not install merged output.")
  cat(paste(status, collapse = "\n"), "\nMerged output:", target, "\n")
  invisible(output)
}
