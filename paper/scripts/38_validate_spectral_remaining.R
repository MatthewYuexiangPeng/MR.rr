#!/usr/bin/env Rscript
# Base tests are genuinely executable without CVXR/mr.divw/MrDAG.
# The full native test also executes the actual packages and PSOCK workers.
remaining_check_comparator_coefficients <- function(saved, direct, expected_length,
                                                   label, tol = 1e-10) {
  # The audited mr.divw functions return K-by-1 beta.hat matrices, whereas
  # indexing a saved C row returns a vector. all.equal(check.attributes=FALSE)
  # still rejects numeric-vs-matrix classes even when every coefficient agrees.
  # Check the documented shape before flattening; do not relax numeric checks.
  valid <- is.numeric(saved) && is.null(dim(saved)) &&
    is.numeric(direct) && (is.null(dim(direct)) ||
      identical(dim(direct), c(as.integer(expected_length), 1L))) &&
    length(saved) == expected_length && length(direct) == expected_length &&
    all(is.finite(saved)) && all(is.finite(direct))
  if (!valid) stop("Invalid comparator coefficient shape or non-finite values: ", label, call. = FALSE)
  same <- all.equal(unname(saved), as.numeric(direct), tolerance = tol, check.attributes = FALSE)
  if (!isTRUE(same)) stop("Comparator coefficient mismatch: ", label,
    "; maximum absolute difference = ", format(max(abs(saved - as.numeric(direct))), digits = 16L),
    "; ", paste(same, collapse = "; "), call. = FALSE)
  invisible(TRUE)
}
validate_remaining <- function(root = getwd(), output = tempfile("remaining_validation_"),
                               base_only = FALSE, cores = 2L, reference = "", reference_bundle = "") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  a <- new.env(parent = baseenv()); sys.source(file.path(root, "paper/lib/paper_remaining.R"), a)
  ctx <- a$rem_load(root); api <- ctx$api; cfg <- api$rem_config(root)
  pass <- function(x) cat("PASS:", x, "\n")
  fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
  close <- function(x, y, tol = 1e-10) stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = tol, check.attributes = FALSE)))
  coefficients <- c(-.4, .2, 0, .7, -.1, .3, .8, -.6, .05)
  remaining_check_comparator_coefficients(coefficients, matrix(coefficients, 9L, 1L), 9L, "matrix fixture")
  remaining_check_comparator_coefficients(coefficients, coefficients, 9L, "vector fixture")
  fails(remaining_check_comparator_coefficients(coefficients, matrix(coefficients, 3L, 3L), 9L, "wrong shape"))
  fails(remaining_check_comparator_coefficients(coefficients, coefficients[-1L], 9L, "wrong length"))
  fails(remaining_check_comparator_coefficients(coefficients, c(Inf, coefficients[-1L]), 9L, "non-finite"))
  fails(remaining_check_comparator_coefficients(coefficients, matrix(rev(coefficients), 9L, 1L), 9L, "wrong order"))
  fails(remaining_check_comparator_coefficients(coefficients, matrix(coefficients + 1e-6, 9L, 1L), 9L, "changed values"))
  pass("Comparator validation accepts equal vector/column-matrix values and rejects wrong shape, order or coefficients")
  d <- api$rem_real_data(root, cfg)
  raw <- utils::read.csv(file.path(root, "paper/input/dat_1e-4.csv"))
  v <- sqrt(2 * raw$ImpMAF * (1 - raw$ImpMAF))
  for (j in 1:9) {close(d$X[, j], raw[[paste0("gamma_exp", j)]] * v); close(d$sx[, j], raw[[paste0("se_exp", j)]] * v)}
  for (j in 1:3) {close(d$Y[, j], raw[[paste0("gamma_out", j + 1L)]] * v); close(d$sy[, j], raw[[paste0("se_out", j + 1L)]] * v)}
  pass("Real beta and SE columns retain every trait and SNP in order")
  idx <- rep(1:89, length.out = 177L); db <- api$rem_resample_real(d, idx)
  hand <- matrix(0, 9L, 9L)
  for (i in idx) hand <- hand + outer(d$sx[i, ], d$sx[i, ]) * d$cor_x
  close(db$Sigma_X, hand / 177L); stopifnot(!identical(d$Sigma_X, db$Sigma_X))
  close(db$W %*% db$Sigma_Y, diag(3L))
  pass("Real bootstrap resamples beta/SE together and recalculates covariance/weights")
  rank <- api$rem_rank(d, ctx$engine, cfg); stopifnot(rank$valid, rank$selected_rank == 1L)
  W <- ctx$engine$estimators$.sqrt_matrix(d$W); XY <- crossprod(d$X, d$Y) / 177L
  lambda <- eigen(W %*% t(XY) %*% solve(crossprod(d$X) / 177L - d$Sigma_X) %*% XY %*% W)$values
  p <- vapply(1:2, function(r) 1 - stats::pchisq((177 - 6.5) * sum(log(1 + lambda[(r + 1L):3L])), (3L-r)*(9L-r)), numeric(1))
  close(rank$p_value, p)
  stopifnot(round(min(api$rem_real_strength(d, ctx$engine)), 2) == 36.86)
  pass("Real rank-one diagnostic matches archived M0 and heteroskedastic SIV is 36.86")
  X <- rbind(sqrt(60) * diag(9L), matrix(0, 51L, 9L))
  fixture <- list(X = X, Y = X[, 1:3] / 10, Sigma_X = diag(c(rep(2, 3L), rep(.5, 6L))), W = diag(3L))
  z <- api$rem_rank(fixture, ctx$engine, cfg)
  stopifnot(z$valid, z$indefinite, all(z$lambda < 0), z$selected_rank == 1L)
  fixture$Y <- X[, 1:3] * 10; z <- api$rem_rank(fixture, ctx$engine, cfg)
  stopifnot(!z$valid, is.na(z$selected_rank), nzchar(z$error), z$archived_selected_rank == 3L)
  pass("Signed rank rule is retained; undefined logarithms are recorded without assigning a rank")
  for (r in 1:2) {
    point <- api$rem_phi(d, r, ctx$engine, cfg); boot <- api$rem_phi(db, r, ctx$engine, cfg, "real_bootstrap")
    legacy <- ctx$engine$estimators$mr_rr_regularized(d$Y, d$X, r, d$Sigma_X,
      regularization_rate = point$selected_phi, W = d$W, implementation = "legacy")
    close(point$fit$AB, legacy$AB, 1e-8)
    stopifnot(point$selected_D == point$grid$D[which.min(point$grid$objective)],
      !identical(point$grid$phi, boot$grid$phi), boot$fit$paper$stage == "real_bootstrap")
  }
  pass("Spectral phi fits match stable legacy examples; bootstrap rebuilds the tuning grid")
  RNGkind("Mersenne-Twister"); set.seed(171L); old <- .Random.seed; kind <- RNGkind()
  rng <- api$rem_rng(cfg)
  states <- c(rng, lapply(1:3, function(b) api$rem_substate(rng[["real/resample"]], b)))
  stopifnot(!anyDuplicated(vapply(states, api$rem_md5, character(1))))
  draw <- function(b) api$paper_sim_with_state(api$rem_substate(rng[["real/resample"]], b), sample.int(177L, 177L, TRUE))
  first <- lapply(1:3, draw); stats::runif(3L); reverse <- lapply(3:1, draw)
  stopifnot(identical(first, rev(reverse)))
  set.seed(171L); api$rem_rng(cfg)
  fails(api$paper_sim_with_state(rng[[1L]], stop("intentional")))
  stopifnot(identical(old, .Random.seed), identical(kind, RNGkind()))
  pass("Data/method streams are distinct; bootstrap prefixes replay and caller RNG is restored")
  t <- api$rem_task_plan(cfg); cat <- api$rem_catalog()
  stopifnot(nrow(t) == 194L, sum(t$kind == "rank") == 80L, sum(t$kind == "sim_eta") == 12L,
    nrow(cat) == 11L, sum(cat$rank == 0L) == 3L)
  for (design in cfg$rank_designs) for (s in 1:4) {
    sub <- t[t$kind == "rank" & t$design == design & t$setting == s, ]
    ids <- unlist(Map(seq.int, sub$first, sub$last)); stopifnot(identical(ids, 1:1000))
  }
  for (res in c("standard", "mrdag")) {
    sub <- t[t$kind == "real_bootstrap" & t$resource == res, ]
    stopifnot(identical(unlist(Map(seq.int, sub$first, sub$last)), 1:1000))
  }
  pass("Production task coverage is complete, with rank-independent real methods shared once")
  bundle <- api$paper_sim_build_bundle(root, ctx$engine)
  for (s in 1:4) {
    pilot <- api$rem_pilot_data(bundle, s, 1L, rng)
    stopifnot(identical(dim(pilot$X), c(177L, 9L)), identical(dim(pilot$Y), c(177L, 3L)),
      identical(pilot, api$rem_pilot_data(bundle, s, 1L, rng)),
      identical(pilot$Sigma_X, bundle$parameters[[paste0("generic/", s)]]$Sigma_X))
    for (design in cfg$rank_designs) {
      x <- api$paper_sim_data(bundle, design, s, 1L); p <- bundle$parameters[[paste(design, s, sep = "/")]]
      rr <- api$rem_rank(c(x[c("X", "Y")], list(Sigma_X = p$Sigma_X, W = p$weight.matrix)), ctx$engine, cfg)
      stopifnot(rr$archived_selected_rank %in% 1:3, rr$valid || nzchar(rr$error))
    }
  }
  pass("Materialized pilot data replay and all eight rank-design cells use the canonical main DGP")
  # A temporary metadata-only fixture exercises the production reference reader.
  # It contains no estimator results and is never passed to prepare or merge.
  # Match the extra provenance/fingerprint added by the actual preparation script.
  extra_source <- "paper/scripts/32_prepare_spectral_simulations.R"
  bundle$source_manifest <- rbind(bundle$source_manifest, data.frame(path = extra_source,
    md5 = unname(tools::md5sum(file.path(root, extra_source))), stringsAsFactors = FALSE))
  bundle$provenance <- list(mode = "REFERENCE_READER_FIXTURE_ONLY")
  fr <- file.path(output, "reference_reader_fixture_only"); dir.create(file.path(fr, "merged"), recursive = TRUE)
  bf <- file.path(fr, "simulation_bundle.rds"); rf <- file.path(fr, "merged/reference.rds")
  saveRDS(bundle, bf, version = 2L)
  ss <- list(schema = "spectral-run-1", validation_only = FALSE,
    bundle_md5 = unname(tools::md5sum(bf)), sources = bundle$source_manifest,
    runtime = list(R = R.version.string))
  ss$run_id <- api$rem_md5(ss)
  hashes <- lapply(names(bundle$parameters), function(cell) {
    v <- strsplit(cell, "/", fixed = TRUE)[[1L]]
    x <- api$paper_sim_data(bundle, v[1L], as.integer(v[2L]), 1L)
    rep(api$rem_md5(x[c("X", "Y")]), 1000L)
  }); names(hashes) <- names(bundle$parameters)
  rr <- c(list(schema = "spectral-merged-1", seal = ss, run_id = ss$run_id, data_md5 = hashes),
    bundle[c("config", "truths", "parameters", "prediction_exposure", "catalog", "table_map")])
  saveRDS(rr, rf, version = 2L); before <- tools::md5sum(c(bf, rf))
  builder <- api$paper_sim_build_bundle
  loaded <- tryCatch({
    api$paper_sim_build_bundle <- function(...) stop("Reference replay attempted to recalibrate the main design.")
    api$rem_reference_inputs(ctx, rf)
  }, finally = {api$paper_sim_build_bundle <- builder})
  stopifnot(identical(loaded$bundle, bundle), identical(loaded$ref, rr),
    identical(loaded$provenance$bundle_file_md5, ss$bundle_md5))
  for (design in cfg$rank_designs) for (s in 1:4) {
    x <- api$paper_sim_data(loaded$bundle, design, s, 1L)
    api$rem_check_reference_data(x, loaded$ref, design, s, 1L)
  }
  fails(api$rem_reference_inputs(ctx, rf, file.path(fr, "missing.rds")))
  bad <- bundle; bad$parameters[[1L]]$C[1L] <- bad$parameters[[1L]]$C[1L] + 1e-8
  badfile <- file.path(fr, "wrong_bundle.rds"); saveRDS(bad, badfile, version = 2L)
  fails(api$rem_reference_inputs(ctx, rf, badfile))
  bad <- rr; bad$parameters[[1L]]$C[1L] <- bad$parameters[[1L]]$C[1L] + 1e-8
  fails(api$rem_validate_reference_bundle(ctx, bad, bundle, ss$bundle_md5))
  bad <- rr; bad$seal$sources$md5[1L] <- strrep("0", 32L)
  fails(api$rem_validate_reference_bundle(ctx, bad, bundle, ss$bundle_md5))
  bad <- rr; bad$seal$runtime$R <- "Wrong R version"; badfile <- file.path(fr, "wrong_R.rds")
  saveRDS(bad, badfile, version = 2L); fails(api$rem_reference_inputs(ctx, badfile, bf))
  x <- api$paper_sim_data(bundle, "generic", 1L, 1L); x$Y[1L] <- x$Y[1L] + 1e-12
  fails(api$rem_check_reference_data(x, rr, "generic", 1L, 1L))
  stopifnot(identical(before, tools::md5sum(c(bf, rf))))
  pass("Saved-design replay never recalibrates; wrong bundles, sources, R versions and changed data are rejected")
  # Exercise real rank workers/checkpoints in a clearly marked, nonproduction fixture.
  dir.create(file.path(output, "fixture"), showWarnings = FALSE)
  frun <- normalizePath(file.path(output, "fixture"), winslash = "/")
  task <- t[1L, ]; task$first <- task$last <- 1L
  task$input_path <- "inputs/test.rds"; task$output_path <- "chunks/task_0001.rds"
  input <- list(data = list(list(X = d$X, Y = d$Y, data_md5 = api$rem_md5(list(X = d$X, Y = d$Y)))),
    parameters = list(Sigma_X = d$Sigma_X, weight.matrix = d$W))
  api$paper_write_rds(input, file.path(frun, task$input_path)); task$input_md5 <- unname(tools::md5sum(file.path(frun, task$input_path)))
  st <- list(run = frun, tasks = task, seal = list(run_id = "VALIDATION_FIXTURE_ONLY", config = cfg, rng = rng))
  x <- api$rem_run_task(ctx, st, 1L); hash <- unname(tools::md5sum(file.path(frun, task$output_path)))
  api$rem_run_task(ctx, st, 1L); stopifnot(identical(hash, unname(tools::md5sum(file.path(frun, task$output_path)))))
  bad <- x; bad$records[[1]]$payload$result$data_md5 <- "corrupt"; fails(api$rem_validate_chunk(bad, st, task))
  bad <- x; bad$unit_ids <- 1:2; fails(api$rem_validate_chunk(bad, st, task))
  removed <- file.path(frun, "removed_chunk.rds"); stopifnot(file.rename(file.path(frun, task$output_path), removed))
  api$rem_run_task(ctx, st, 1L); stopifnot(identical(hash, unname(tools::md5sum(file.path(frun, task$output_path)))))
  input$data[[1]]$X[1, 1] <- 123; api$paper_write_rds(input, file.path(frun, task$input_path))
  fails(api$rem_input(st, task))
  pass("Checkpoints resume exactly; changed inputs/results and overlapping chunk coverage are rejected")
  if (nzchar(reference)) {
    main <- api$rem_reference_inputs(ctx, reference, reference_bundle)
    for (design in cfg$rank_designs) for (s in 1:4) for (i in 1:2) {
      x <- api$paper_sim_data(main$bundle, design, s, i)
      api$rem_check_reference_data(x, main$ref, design, s, i)
    }
    stopifnot(identical(main$file_hashes, tools::md5sum(main$files)))
    pass("Original sealed main bundle replays all 16 development datasets exactly; reference files unchanged")
  }
  if (!base_only) {
    api$paper_external_require(api$paper_sim_methods())
    run <- file.path(output, "native_development")
    api$rem_prepare(ctx, run, "development", reference, reference_bundle)
    state <- api$rem_open(ctx, run); real <- readRDS(file.path(run, "inputs/real.rds"))
    for (m in c("ivw", "srivw")) {
      f <- getExportedValue("mr.divw", if (m == "ivw") "mvmr.ivw" else "mvmr.divw")
      for (j in 1:3) {
        args <- list(beta.exposure = d$X, se.exposure = d$sx,
          beta.outcome = d$Y[, j], se.outcome = d$sy[, j], gen_cor = d$cor_x)
        if (m == "srivw") args["phi_cand"] <- list(NULL)
        direct <- do.call(f, args)$beta.hat
        remaining_check_comparator_coefficients(real$point[[paste0(m, "__0")]]$AB[j, ],
          direct, ncol(d$X), paste(m, "outcome", j))
      }
    }
    pass("Real IVW/SRIVW adapter agrees with direct native calls including heteroskedastic SE and gen_cor")
    for (id in state$tasks$task_id) api$rem_run_task(ctx, state, id, cores)
    # Compare each resource class with direct serial calls for identical canonical inputs.
    for (res in c("standard", "mrdag")) {
      tt <- state$tasks[state$tasks$kind == "real_bootstrap" & state$tasks$resource == res, ][1L, ]
      input <- api$rem_input(state, tt); saved <- readRDS(file.path(run, tt$output_path))
      for (j in seq_along(saved$unit_ids)) {
        serial <- api$rem_unit(ctx, state, tt, input, saved$unit_ids[j])
        stopifnot(identical(serial, saved$records[[j]]))
      }
    }
    pass("Serial and PSOCK bootstrap results agree for actual sparse and external comparator fits")
    bad <- state; bad$seal$run_id <- "wrong-run"
    stopifnot(all(api$rem_inventory(bad, write = FALSE)$status == "invalid"))
    stopifnot(all(api$rem_inventory(state)$status == "valid"))
    merged <- api$rem_merge(ctx, state)
    stopifnot(!merged$production, length(merged$real_bootstrap) == 11L,
      all(vapply(merged$real_bootstrap, function(x) identical(dim(x), c(27L, 3L)), logical(1))))
    pass("Native development run strictly merges complete rank/tuning/real bootstrap outputs")
  }
  status <- if (base_only) "SPECTRAL REMAINING VALIDATION: PARTIAL_PASS (native packages and PSOCK require cluster preflight)"
    else "SPECTRAL REMAINING VALIDATION: PASS (native packages and PSOCK exercised)"
  writeLines(c(status, paste("R:", R.version.string)), file.path(output, "STATUS.txt"))
  utils::write.csv(api$rem_sources(ctx), file.path(output, "source_fingerprints.csv"), row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
  cat(status, "\nReport:", output, "\n"); invisible(TRUE)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); opt <- function(k, default) {
    x <- args[startsWith(args, paste0("--", k, "="))]
    if (length(x)) sub("^[^=]+=", "", x[1L]) else default
  }
  validate_remaining(opt("root", getwd()), opt("output", file.path(getwd(), "paper/output/spectral_rebuild/remaining_validation")),
    "--base-only" %in% args, as.integer(opt("cores", "2")), opt("reference", ""), opt("reference-bundle", ""))
}
