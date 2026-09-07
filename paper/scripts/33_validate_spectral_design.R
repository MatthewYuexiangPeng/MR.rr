#!/usr/bin/env Rscript
# Base-R numerical/data/RNG/manifest checks; no CVXR or production fits required.
paper_validate_spectral_design <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  api <- new.env(parent = baseenv())
  for (name in c("paper_engine.R", "paper_simulation.R"))
    sys.source(file.path(root, "paper/lib", name), envir = api)
  engine <- api$paper_engine_load(root)
  bundle <- api$paper_sim_build_bundle(root, engine, replicates = 3L, bootstrap_size = 5L)
  reference <- readRDS(file.path(root, "paper/tests/fixtures/simulation_design_reference.rds"))
  output <- file.path(root, "paper/output/spectral_rebuild/design_validation")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  check <- function(label, expression) {
    error <- tryCatch({force(expression); NULL}, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <<- data.frame(check = label,
      status = if (is.null(error)) "PASS" else "FAIL", detail = if (is.null(error)) "" else error)
    cat(if (is.null(error)) "PASS: " else "FAIL: ", label,
        if (!is.null(error)) paste0(" - ", error) else "", "\n", sep = "")
  }
  close_to <- function(x, y, tolerance = 1e-10) {
    stopifnot(identical(dim(x), dim(y)), all(is.finite(x)),
              max(abs(x - y)) <= tolerance * max(1, max(abs(y))))
  }
  throws <- function(expression) stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
  api$paper_with_seed(9062026L, {
    check("Raw-input calibration matches archived population matrices", {
      close_to(bundle$calibration$sigma_x_unweighted, reference$Sigma_X)
      close_to(bundle$calibration$sigma_y, reference$Sigma_Y)
      close_to(bundle$calibration$sigma_gg_unweighted, reference$VX_tilde)
      stopifnot(bundle$calibration$pz == 177L, bundle$calibration$px == 9L, bundle$calibration$py == 3L)
    })
    check("Generic and sparse-loading truth match archived designs", {
      close_to(bundle$truths$generic$C, reference$generic_C)
      close_to(bundle$truths$sparse_loading$C, reference$sparse_C)
      S <- bundle$truths$sparse_loading
      close_to(crossprod(S$A, bundle$calibration$weight_matrix %*% S$A), diag(2))
      stopifnot(all(colSums(S$B != 0) == 1L))
    })
    check("Approximate truth has singular values 1, 1, 0.1", {
      close_to(svd(bundle$truths$approximate_low_rank$C)$d, c(1, 1, .1))
      close_to(bundle$truths$generic$U, bundle$truths$approximate_low_rank$U)
      close_to(bundle$truths$generic$V, bundle$truths$approximate_low_rank$V)
    })
    check("Four population SIV values and covariance factors", {
      actual <- vapply(bundle$parameters[paste("generic", 1:4, sep = "/")],
                        function(x) x$iv_strength, numeric(1))
      close_to(unname(actual), unname(reference$generic_siv))
      for (p in bundle$parameters) for (name in c("Sigma_X", "Sigma_Y", "VX_tilde"))
        close_to(crossprod(chol(p[[name]])), p[[name]])
    })
    check("One rank-two result key is used by both tables", {
      stopifnot(nrow(bundle$catalog) == 100L, !anyDuplicated(bundle$catalog$result_key),
                sum(bundle$catalog$bootstrap) == 80L, nrow(bundle$table_map) == 108L)
      main <- bundle$table_map[bundle$table_map$table_id == "main_generic" &
        bundle$table_map$method %in% c("regularized_mr_rr", "sparse_mr_rr"), ]
      sensitivity <- bundle$table_map[bundle$table_map$table_id == "rank_misspecification" &
        bundle$table_map$working_rank == 2L, ]
      both <- merge(main, sensitivity, by = c("setting", "method"))
      stopifnot(nrow(both) == 8L,
                identical(both$point_result_key.x, both$point_result_key.y),
                identical(both$bootstrap_result_key.x, both$bootstrap_result_key.y))
    })
    check("Data replay survives task reordering and intervening RNG use", {
      first <- api$paper_sim_data(bundle, "generic", 4L, 3L)
      api$paper_sim_data(bundle, "sparse_loading", 1L, 1L)
      stats::rnorm(100)
      second <- api$paper_sim_data(bundle, "generic", 4L, 3L)
      stopifnot(identical(first, second))
      close_to(first$Y, first$latent_X %*% t(bundle$truths$generic$C) + first$error_Y)
      close_to(first$X, first$latent_X + first$error_X)
    })
    check("CMRG streams are distinct across data, resampling, methods and ranks", {
      states <- list()
      for (d in api$paper_sim_designs()) for (s in c(1L, 4L)) for (r in c(1L, 3L)) {
        states[[length(states) + 1L]] <- api$paper_sim_stream(bundle, d, s, r, "data")
        states[[length(states) + 1L]] <- api$paper_sim_stream(bundle, d, s, r, "resample")
        for (rank in 1:3) for (method in c("regularized_mr_rr", "sparse_mr_rr"))
          states[[length(states) + 1L]] <- api$paper_sim_stream(bundle, d, s, r, "method", method, rank)
      }
      stopifnot(!anyDuplicated(vapply(states, paste, collapse = ",", character(1))))
    })
    check("Bootstrap resamples are shared and preserve prefixes when B changes", {
      a <- api$paper_sim_resamples(bundle, "generic", 2L, 1L, B = 5L)
      api$paper_sim_method_states(bundle, "generic", 2L, 1L, "regularized_mr_rr", 1L)
      b <- api$paper_sim_resamples(bundle, "generic", 2L, 1L, B = 3L)
      stopifnot(identical(a[, 1:3, drop = FALSE], b), all(a >= 1L & a <= 177L))
      t1 <- api$paper_sim_method_states(bundle, "generic", 2L, 1L, "regularized_mr_rr", 2L, 5L)
      t2 <- api$paper_sim_method_states(bundle, "generic", 2L, 1L, "regularized_mr_rr", 2L, 3L)
      stopifnot(identical(t1[, 1:4, drop = FALSE], t2))
    })
    check("Data and bootstrap generation restore the caller RNG", {
      before <- get(".Random.seed", envir = .GlobalEnv)
      kinds <- RNGkind()
      api$paper_sim_data(bundle, "generic", 1L, 1L)
      api$paper_sim_resamples(bundle, "generic", 1L, 1L)
      stopifnot(identical(before, get(".Random.seed", envir = .GlobalEnv)), identical(kinds, RNGkind()))
    })
    check("Spectral point and bootstrap estimates use one shared rank-two object", {
      dat <- api$paper_sim_data(bundle, "generic", 4L, 1L)
      cfg <- bundle$config[bundle$config$design == "generic" & bundle$config$setting == 4L, ]
      p <- bundle$parameters[["generic/4"]]
      key <- bundle$catalog$result_key[bundle$catalog$design == "generic" &
        bundle$catalog$setting == 4L & bundle$catalog$method == "regularized_mr_rr" &
        bundle$catalog$working_rank == 2L]
      fit <- function(X, Y) api$paper_engine_fit(engine, "regularized", "simulation", Y, X, 2L,
        p$Sigma_X, p$weight.matrix, regularization_rate = cfg$regularization_rate)$AB
      indices <- api$paper_sim_resamples(bundle, "generic", 4L, 1L)
      results <- setNames(list(list(point = fit(dat$X, dat$Y), bootstrap = lapply(1:5, function(b)
        fit(dat$X[indices[, b], , drop = FALSE], dat$Y[indices[, b], , drop = FALSE])))), key)
      keys <- bundle$table_map$point_result_key[bundle$table_map$setting == 4L &
        bundle$table_map$design == "generic" & bundle$table_map$method == "regularized_mr_rr" &
        bundle$table_map$working_rank == 2L]
      stopifnot(length(keys) == 2L, identical(results[[keys[1L]]], results[[keys[2L]]]))
    })
    check("Development task coverage is unchanged by chunk sizes", {
      api$paper_sim_validate_tasks(bundle, api$paper_sim_tasks(bundle))
      api$paper_sim_validate_tasks(bundle, api$paper_sim_tasks(bundle, 1L, 2L, 2L, 3L))
    })
    full <- bundle
    full$replicates <- 1000L
    full$bootstrap_size <- 300L
    full$config$replicates <- 1000L
    full$config$bootstrap_size <- 300L
    tasks <- api$paper_sim_tasks(full)
    check("Full manifests cover every required fit exactly once", {
      api$paper_sim_validate_tasks(full, tasks)
      stopifnot(nrow(tasks) == 1800L)
    })
    check("Missing and overlapping tasks are rejected", {
      small <- api$paper_sim_tasks(bundle)
      throws(api$paper_sim_validate_tasks(bundle, small[-1L, ]))
      bad <- small
      bad$result_keys[1L] <- paste(bad$result_keys[1L], strsplit(bad$result_keys[1L], ";", fixed = TRUE)[[1]][1], sep = ";")
      throws(api$paper_sim_validate_tasks(bundle, bad))
    })
    check("Prepared bundle survives save/read and does not require frozen results", {
      f <- tempfile(fileext = ".rds")
      on.exit(unlink(f), add = TRUE)
      saveRDS(bundle, f, version = 2)
      restored <- readRDS(f)
      stopifnot(identical(restored, bundle),
                !any(grepl("freeze|\\.RData$", bundle$source_manifest$path)))
    })
  })
  out <- do.call(rbind, rows)
  utils::write.csv(out, file.path(output, "checks.csv"), row.names = FALSE)
  utils::write.csv(bundle$source_manifest, file.path(output, "sources.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
  ok <- all(out$status == "PASS")
  writeLines(if (ok) "SPECTRAL DESIGN VALIDATION: PASS" else "SPECTRAL DESIGN VALIDATION: FAIL",
             file.path(output, "STATUS.txt"))
  if (!ok) stop("Design validation failed; inspect checks.csv.")
  cat("SPECTRAL DESIGN VALIDATION: PASS\nReport:", normalizePath(output, winslash = "/"), "\n")
  invisible(out)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) stop("Run without arguments from the repository root.")
  paper_validate_spectral_design()
}
