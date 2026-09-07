# Calibrated Gaussian simulation design for the spectral paper rebuild.
# Calibration/population formulas are retained from the audited baseline.
# Runtime inputs are paper/input CSV files; no frozen code/results are read.

paper_sim_assert_finite <- function(x, label) {
  if (!is.numeric(x) || any(!is.finite(x))) stop(label, " must be finite numeric data.", call. = FALSE)
}

paper_sim_second_moment <- function(value) {
  crossprod(value) / nrow(value)
}


paper_sim_population_siv <- function(parameters, pz = 177L) {
  sigma_x <- (parameters$Sigma_X + t(parameters$Sigma_X)) / 2
  decomposition <- eigen(sigma_x, symmetric = TRUE)

  if (any(!is.finite(decomposition$values)) ||
      any(decomposition$values <= 0)) {
    stop(
      "Sigma_X is not positive definite.",
      call. = FALSE
    )
  }

  sigma_x_inverse_sqrt <-
    decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values)) %*%
    t(decomposition$vectors)

  standardized_strength <-
    sigma_x_inverse_sqrt %*%
    parameters$VX_tilde %*%
    sigma_x_inverse_sqrt
  standardized_strength <-
    (standardized_strength + t(standardized_strength)) / 2

  eigenvalues <- eigen(
    standardized_strength,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  sqrt(pz) * min(eigenvalues)
}


paper_sim_calibration <- function(lip_data, lip_correlation, pz = 177L) {
  px <- 9L
  py <- 3L

  if (!identical(nrow(lip_data), pz)) {
    stop(
      "The paper simulation calibration requires exactly 177 instruments; ",
      "the raw input contains ",
      nrow(lip_data),
      ".",
      call. = FALSE
    )
  }

  required_columns <- c(
    "ImpMAF",
    paste0("gamma_exp", seq_len(px)),
    paste0("se_exp", seq_len(px)),
    paste0("se_out", 2:4)
  )
  missing_columns <- setdiff(required_columns, names(lip_data))

  if (length(missing_columns) > 0L) {
    stop(
      "The raw simulation input is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(lip_correlation) < 13L ||
      ncol(lip_correlation) < 13L) {
    stop(
      "The raw correlation matrix must be at least 13 x 13.",
      call. = FALSE
    )
  }

  variant_variance <-
    2 * lip_data$ImpMAF * (1 - lip_data$ImpMAF)
  sqrt_variant_variance <- sqrt(variant_variance)

  if (any(!is.finite(variant_variance)) || any(variant_variance <= 0)) {
    stop("ImpMAF must be finite and strictly between zero and one.", call. = FALSE)
  }
  if (any(!is.finite(as.matrix(lip_correlation))) ||
      max(abs(as.matrix(lip_correlation) - t(as.matrix(lip_correlation)))) > 1e-10) {
    stop("Correlation input must be finite and symmetric.", call. = FALSE)
  }
  gamma_observed <- as.matrix(
    lip_data[, paste0("gamma_exp", seq_len(px))]
  )
  gamma_star_calibration <-
    gamma_observed * sqrt_variant_variance

  se_exposure <- as.matrix(
    lip_data[, paste0("se_exp", seq_len(px))]
  )
  se_exposure <- se_exposure * sqrt_variant_variance
  correlation_x <- as.matrix(
    lip_correlation[seq_len(px), seq_len(px)]
  )

  sigma_x_by_instrument <- lapply(
    seq_len(pz),
    function(instrument) {
      diagonal <- diag(se_exposure[instrument, ])
      diagonal %*% correlation_x %*% diagonal
    }
  )
  sigma_x_unweighted <-
    Reduce(`+`, sigma_x_by_instrument) / pz
  sigma_gg_unweighted <-
    paper_sim_second_moment(gamma_star_calibration) - sigma_x_unweighted

  se_outcome <- as.matrix(
    lip_data[, paste0("se_out", 2:4)]
  )
  se_outcome <- se_outcome * sqrt_variant_variance
  correlation_y <- as.matrix(
    lip_correlation[11:13, 11:13]
  )

  sigma_y_by_instrument <- lapply(
    seq_len(pz),
    function(instrument) {
      diagonal <- diag(se_outcome[instrument, ])
      diagonal %*% correlation_y %*% diagonal
    }
  )
  sigma_y <- Reduce(`+`, sigma_y_by_instrument) / pz

  values_to_check <- list(
    variant_variance = variant_variance,
    sigma_x_unweighted = sigma_x_unweighted,
    sigma_gg_unweighted = sigma_gg_unweighted,
    sigma_y = sigma_y
  )

  for (value_name in names(values_to_check)) {
    paper_sim_assert_finite(values_to_check[[value_name]], value_name)
  }
  for (name in c("sigma_x_unweighted", "sigma_gg_unweighted", "sigma_y")) {
    tryCatch(chol(values_to_check[[name]]), error = function(e) {
      stop("Calibration covariance is not positive definite: ", name, call. = FALSE)
    })
  }

  list(
    pz = pz,
    px = px,
    py = py,
    variant_variance = variant_variance,
    sigma_x_unweighted = sigma_x_unweighted,
    sigma_gg_unweighted = sigma_gg_unweighted,
    sigma_y = sigma_y,
    weight_matrix = solve(sigma_y)
  )
}



paper_sim_parameters <- function(
    causal_effect,
    exposure_error_weight,
    genetic_effect_weight,
    calibration,
    estimators,
    rank = 2L) {
  sigma_x <-
    exposure_error_weight * calibration$sigma_x_unweighted
  sigma_gg <-
    genetic_effect_weight * calibration$sigma_gg_unweighted
  sigma_y <- calibration$sigma_y
  weight_matrix <- solve(sigma_y)

  sigma_xx <- sigma_x + sigma_gg
  sigma_yx <- causal_effect %*% sigma_gg
  sigma_xy <- t(sigma_yx)
  weight_sqrt <- estimators$.sqrt_matrix(weight_matrix)
  weight_sqrt_inverse <- solve(weight_sqrt)

  naive_target <-
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_xx) %*%
    sigma_xy %*%
    weight_sqrt
  naive_vectors <- eigen(naive_target)$vectors[
    , seq_len(rank), drop = FALSE
  ]
  A <- weight_sqrt_inverse %*% naive_vectors
  B <-
    t(naive_vectors) %*%
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_xx)

  corrected_target <-
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_gg) %*%
    sigma_xy %*%
    weight_sqrt
  corrected_vectors <- eigen(corrected_target)$vectors[
    , seq_len(rank), drop = FALSE
  ]
  A_corrected <- weight_sqrt_inverse %*% corrected_vectors
  B_corrected <-
    t(corrected_vectors) %*%
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_gg)

  parameters <- list(
    py = calibration$py,
    px = calibration$px,
    var_Z = calibration$variant_variance,
    VX_tilde = sigma_gg,
    Sigma_X = sigma_x,
    Sigma_Y = sigma_y,
    weight.matrix = weight_matrix,
    SigmaXX = sigma_xx,
    SigmaYX = sigma_yx,
    SigmaXY = sigma_xy,
    C = causal_effect,
    A = A,
    B = B,
    A_d = A_corrected,
    B_d = B_corrected,
    C_r = A_corrected %*% B_corrected,
    r_RR = rank,
    VY_tilde = causal_effect %*% sigma_gg %*% t(causal_effect)
  )

  parameters$iv_strength <- paper_sim_population_siv(
    parameters,
    pz = calibration$pz
  )
  parameters
}


paper_sim_designs <- function() c("generic", "sparse_loading", "approximate_low_rank")
paper_sim_methods <- function() c("ivw", "srivw", "naive_mr_rr", "mr_rr",
                                  "regularized_mr_rr", "sparse_mr_rr", "mrdag")
paper_sim_integer <- function(x, label, minimum = 1L, maximum = 1000L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < minimum || x > maximum) stop("Invalid ", label, call. = FALSE)
  as.integer(x)
}

paper_sim_read_config <- function(root) {
  cfg <- utils::read.csv(file.path(root, "paper/config/spectral_simulations.csv"),
                         stringsAsFactors = FALSE)
  required <- c("design", "setting", "measurement_error_weight", "genetic_effect_weight",
    "truth_rank", "third_singular_value", "regularization_rate", "sparse_lambda", "replicates",
    "bootstrap_size", "master_seed", "truth_seed", "mrdag_niter", "mrdag_burnin",
    "sparse_threshold", "sparse_max_iter", "sparse_tol", "sparse_solver")
  if (!setequal(names(cfg), required) || nrow(cfg) != 12L ||
      anyDuplicated(paste(cfg$design, cfg$setting))) stop("Invalid simulation configuration schema.")
  for (name in setdiff(required, c("design", "sparse_solver")))
    paper_sim_assert_finite(cfg[[name]], name)
  expected <- expand.grid(setting = 1:4, design = paper_sim_designs(), stringsAsFactors = FALSE)
  if (!setequal(paste(cfg$design, cfg$setting), paste(expected$design, expected$setting)))
    stop("Configuration must contain all three designs and four settings.")
  cfg <- cfg[order(match(cfg$design, paper_sim_designs()), cfg$setting), ]
  rownames(cfg) <- NULL
  stopifnot(all(cfg$measurement_error_weight == rep(c(2.5, 1, 2.5, 1), 3)),
            all(cfg$genetic_effect_weight == rep(c(.25, .25, 1, 1), 3)),
            all(cfg$regularization_rate >= 0), all(cfg$sparse_lambda >= 0),
            all(cfg$sparse_threshold >= 0), all(cfg$sparse_tol > 0),
            all(cfg$sparse_solver == "OSQP"))
  for (name in c("replicates", "bootstrap_size", "master_seed", "truth_seed")) {
    if (length(unique(cfg[[name]])) != 1L) stop("All cells must share ", name)
    limit <- if (name %in% c("master_seed", "truth_seed")) .Machine$integer.max else 1000L
    paper_sim_integer(cfg[[name]][1L], name, maximum = limit)
  }
  for (i in seq_len(nrow(cfg))) {
    paper_sim_integer(cfg$sparse_max_iter[i], "sparse_max_iter", maximum = 100000L)
    paper_sim_integer(cfg$mrdag_niter[i], "mrdag_niter", maximum = 10000000L)
    paper_sim_integer(cfg$mrdag_burnin[i], "mrdag_burnin", minimum = 0L,
                       maximum = cfg$mrdag_niter[i] - 1L)
  }
  stopifnot(all(cfg$bootstrap_size >= 2L))
  approx <- cfg$design == "approximate_low_rank"
  stopifnot(all(cfg$third_singular_value[!approx] == 0), all(cfg$truth_rank[!approx] == 2),
            all(cfg$truth_rank[approx] == 3), all(cfg$third_singular_value[approx] > 0),
            all(cfg$third_singular_value[approx] < 1),
            length(unique(cfg$third_singular_value[approx])) == 1L)
  cfg
}

paper_sim_truth <- function(design, calibration, engine, seed = 123L, delta = 0) {
  design <- match.arg(design, paper_sim_designs())
  paper_with_seed(seed, {
    if (design == "sparse_loading") {
      Q <- qr.Q(qr(matrix(stats::rnorm(6), 3, 2)))
      A <- solve(engine$estimators$.sqrt_matrix(calibration$weight_matrix)) %*% Q
      B <- matrix(0, 2, 9)
      for (j in 1:9) {
        selected_row <- sample.int(2L, 1L)
        B[selected_row, j] <- stats::rnorm(1, sd = 25)
      }
      list(C = A %*% B, A = A, B = B)
    } else {
      d <- svd(matrix(stats::rnorm(27), 3, 9))
      values <- c(1, 1, if (design == "generic") 0 else delta)
      list(C = d$u %*% diag(values) %*% t(d$v), U = d$u, V = d$v,
           singular_values = values)
    }
  }, kind = "Mersenne-Twister")
}

paper_sim_result_catalog <- function(cfg) {
  rows <- list()
  for (i in seq_len(nrow(cfg))) for (method in paper_sim_methods()) {
    r <- if (method %in% c("ivw", "srivw", "mrdag")) 0L else 2L
    if (cfg$design[i] == "generic" && method %in% c("regularized_mr_rr", "sparse_mr_rr")) r <- 1:3
    for (rank in r) {
      key <- paste(cfg$design[i], paste0("setting-", cfg$setting[i]), method,
                    paste0("rank-", rank), sep = "__")
      rows[[length(rows) + 1L]] <- data.frame(
        result_key = key, design = cfg$design[i], setting = cfg$setting[i],
        method = method, working_rank = rank, truth_rank = cfg$truth_rank[i],
        resource_class = if (method == "mrdag") "mrdag" else "standard",
        bootstrap = method != "sparse_mr_rr", stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

paper_sim_table_map <- function(catalog) {
  main <- catalog[catalog$working_rank %in% c(0L, 2L), ]
  main$table_id <- c(generic = "main_generic", sparse_loading = "main_sparse_loading",
                     approximate_low_rank = "approximate_low_rank")[main$design]
  sensitivity <- catalog[catalog$design == "generic" &
    catalog$method %in% c("regularized_mr_rr", "sparse_mr_rr"), ]
  sensitivity$table_id <- "rank_misspecification"
  out <- rbind(main, sensitivity)
  out$point_result_key <- out$result_key
  out$bootstrap_result_key <- ifelse(out$bootstrap, out$result_key, NA_character_)
  rownames(out) <- NULL
  out
}

paper_sim_make_registry <- function(master_seed) {
  # Fixed layout independent of requested replicate count, B, task order or cores.
  # 3 designs x 4 settings x 1000 replicate slots x 24 independent CMRG streams.
  paper_with_seed(master_seed, {
    state <- get(".Random.seed", envir = .GlobalEnv)
    registry <- matrix(NA_integer_, nrow = 7L, ncol = 3L * 4L * 1000L * 24L)
    for (j in seq_len(ncol(registry))) {
      registry[, j] <- state
      state <- parallel::nextRNGStream(state)
    }
    registry
  })
}

paper_sim_stream <- function(bundle, design, setting, replicate, purpose, method = NULL, rank = 0L) {
  d <- match(design, paper_sim_designs())
  if (is.na(d)) stop("Unknown design.")
  s <- paper_sim_integer(setting, "setting", maximum = 4L)
  rep <- paper_sim_integer(replicate, "replicate", maximum = 1000L)
  purpose <- match.arg(purpose, c("data", "resample", "method"))
  slot <- switch(purpose, data = 0L, resample = 1L, method = {
    m <- match(method, paper_sim_methods())
    if (length(m) != 1L || is.na(m)) stop("Unknown stochastic method key.")
    if (method %in% c("ivw", "srivw", "mrdag")) {
      paper_sim_integer(rank, "non-rank method rank", minimum = 0L, maximum = 0L)
      r <- 1L
    } else r <- paper_sim_integer(rank, "working rank", maximum = 3L)
    2L + (m - 1L) * 3L + (r - 1L)
  })
  index <- (((d - 1L) * 4L + s - 1L) * 1000L + rep - 1L) * 24L + slot + 1L
  bundle$rng$streams[, index]
}

paper_sim_with_state <- function(state, expression) {
  if (!is.integer(state) || length(state) != 7L || anyNA(state) || state[1L] != 10407L)
    stop("Expected an Inversion/Rejection L'Ecuyer-CMRG state.")
  paper_with_seed(123L, {
    assign(".Random.seed", state, envir = .GlobalEnv)
    force(expression)
  })
}

paper_sim_substreams <- function(state, B) {
  B <- paper_sim_integer(B, "bootstrap size", maximum = 1000L)
  out <- matrix(NA_integer_, 7, B + 1L)
  out[, 1L] <- state
  for (b in seq_len(B)) out[, b + 1L] <- parallel::nextRNGSubStream(out[, b])
  out
}

paper_sim_data <- function(bundle, design, setting, replicate) {
  cell <- paste(design, setting, sep = "/")
  parameters <- bundle$parameters[[cell]]
  if (is.null(parameters)) stop("Unknown design/setting cell.")
  state <- paper_sim_stream(bundle, design, setting, replicate, "data")
  paper_sim_with_state(state, {
    n <- length(parameters$var_Z)
    gaussian <- function(S) matrix(stats::rnorm(n * nrow(S)), n, nrow(S)) %*% chol(S)
    latent <- gaussian(parameters$VX_tilde)
    error_x <- gaussian(parameters$Sigma_X)
    error_y <- gaussian(parameters$Sigma_Y)
    list(X = latent + error_x, Y = latent %*% t(parameters$C) + error_y,
         latent_X = latent, latent_Y = latent %*% t(parameters$C),
         error_X = error_x, error_Y = error_y, rng_state = state,
         design = design, setting = as.integer(setting), replicate = as.integer(replicate))
  })
}

paper_sim_resamples <- function(bundle, design, setting, replicate, B = bundle$bootstrap_size) {
  n <- bundle$calibration$pz
  states <- paper_sim_substreams(paper_sim_stream(bundle, design, setting, replicate, "resample"), B)
  vapply(seq_len(B), function(b) paper_sim_with_state(states[, b + 1L],
    sample.int(n, n, replace = TRUE)), integer(n))
}

paper_sim_method_states <- function(bundle, design, setting, replicate, method, rank,
                                    B = bundle$bootstrap_size) {
  paper_sim_substreams(paper_sim_stream(bundle, design, setting, replicate,
                                        "method", method, rank), B)
}

paper_sim_build_bundle <- function(root, engine, replicates = NULL, bootstrap_size = NULL) {
  cfg <- paper_sim_read_config(root)
  if (!is.null(replicates)) cfg$replicates <- paper_sim_integer(replicates, "replicates")
  if (!is.null(bootstrap_size)) cfg$bootstrap_size <- paper_sim_integer(bootstrap_size, "bootstrap size", 2L)
  dat <- utils::read.csv(file.path(root, "paper/input/dat_1e-4.csv"), check.names = FALSE)
  rho <- utils::read.csv(file.path(root, "paper/input/rho_mat_1e-4.csv"), check.names = FALSE)
  calibration <- paper_sim_calibration(dat, rho)
  truths <- stats::setNames(lapply(paper_sim_designs(), function(d) {
    row <- cfg[cfg$design == d, ][1L, ]
    paper_sim_truth(d, calibration, engine, row$truth_seed, row$third_singular_value)
  }), paper_sim_designs())
  parameters <- stats::setNames(lapply(seq_len(nrow(cfg)), function(i) {
    paper_sim_parameters(truths[[cfg$design[i]]]$C, cfg$measurement_error_weight[i],
      cfg$genetic_effect_weight[i], calibration, engine$estimators, rank = 2L)
  }), paste(cfg$design, cfg$setting, sep = "/"))
  catalog <- paper_sim_result_catalog(cfg)
  input_paths <- c("paper/config/spectral_simulations.csv", "paper/input/dat_1e-4.csv",
                   "paper/input/rho_mat_1e-4.csv", "paper/lib/paper_simulation.R")
  structure(list(schema = "spectral-simulation-design-1", config = cfg,
    replicates = as.integer(cfg$replicates[1L]), bootstrap_size = as.integer(cfg$bootstrap_size[1L]),
    calibration = calibration, truths = truths, parameters = parameters,
    prediction_exposure = paper_with_seed(123L, stats::rnorm(9), kind = "Mersenne-Twister"),
    catalog = catalog, table_map = paper_sim_table_map(catalog),
    rng = list(version = "cmrg-24-streams-per-replicate-v1", master_seed = cfg$master_seed[1L],
      streams = paper_sim_make_registry(cfg$master_seed[1L])),
    source_manifest = rbind(engine$source_manifest, data.frame(path = input_paths,
      md5 = unname(tools::md5sum(file.path(root, input_paths))), stringsAsFactors = FALSE))
  ), class = "mrrr_simulation_bundle")
}

paper_sim_tasks <- function(bundle, point_standard_chunk = 100L, point_mrdag_chunk = 25L,
                            bootstrap_standard_chunk = 20L, bootstrap_mrdag_chunk = 20L) {
  sizes <- c(point_standard_chunk, point_mrdag_chunk, bootstrap_standard_chunk, bootstrap_mrdag_chunk)
  for (x in sizes) paper_sim_integer(x, "chunk size")
  rows <- list()
  for (phase in c("point", "bootstrap")) for (resource in c("standard", "mrdag")) {
    chunk <- sizes[match(paste(phase, resource), c("point standard", "point mrdag",
                                                "bootstrap standard", "bootstrap mrdag"))]
    array_index <- 0L
    for (i in seq_len(nrow(bundle$config))) {
      cfg <- bundle$config[i, ]
      catalog <- bundle$catalog[bundle$catalog$design == cfg$design &
        bundle$catalog$setting == cfg$setting & bundle$catalog$resource_class == resource, ]
      if (phase == "bootstrap") catalog <- catalog[catalog$bootstrap, ]
      if (!nrow(catalog)) next
      for (start in seq.int(1L, bundle$replicates, by = chunk)) {
        array_index <- array_index + 1L
        end <- min(start + chunk - 1L, bundle$replicates)
        key <- sprintf("%s__%s__%s__setting-%d__rep-%04d-%04d", phase, resource,
                       cfg$design, cfg$setting, start, end)
        rows[[length(rows) + 1L]] <- data.frame(
          task_id = length(rows) + 1L, phase = phase, resource_class = resource,
          array_index = array_index, design = cfg$design, setting = cfg$setting,
          replicate_start = start, replicate_end = end,
          bootstrap_size = if (phase == "bootstrap") bundle$bootstrap_size else 0L,
          result_keys = paste(catalog$result_key, collapse = ";"),
          output_relative_path = paste0("chunks/", key, ".rds"),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

paper_sim_validate_tasks <- function(bundle, tasks) {
  stopifnot(nrow(tasks) > 0L, !anyDuplicated(tasks$task_id),
            identical(tasks$task_id, seq_len(nrow(tasks))),
            !anyDuplicated(tasks$output_relative_path),
            all(tasks$phase %in% c("point", "bootstrap")),
            all(tasks$resource_class %in% c("standard", "mrdag")))
  for (phase in c("point", "bootstrap")) for (resource in c("standard", "mrdag")) {
    selected <- tasks[tasks$phase == phase & tasks$resource_class == resource, ]
    stopifnot(identical(selected$array_index, seq_len(nrow(selected))))
  }
  counts <- matrix(0L, nrow(bundle$catalog), bundle$replicates)
  for (phase in c("point", "bootstrap")) {
    counts[,] <- 0L
    for (i in which(tasks$phase == phase)) {
      row <- tasks[i, ]
      ids <- seq.int(row$replicate_start, row$replicate_end)
      stopifnot(row$replicate_start >= 1L, row$replicate_end <= bundle$replicates,
                row$replicate_start <= row$replicate_end,
                row$bootstrap_size == if (phase == "bootstrap") bundle$bootstrap_size else 0L)
      keys <- strsplit(row$result_keys, ";", fixed = TRUE)[[1L]]
      expected <- bundle$catalog[bundle$catalog$design == row$design &
        bundle$catalog$setting == row$setting & bundle$catalog$resource_class == row$resource_class, ]
      if (phase == "bootstrap") expected <- expected[expected$bootstrap, ]
      stopifnot(length(keys) > 0L, !anyDuplicated(keys), setequal(keys, expected$result_key))
      index <- match(keys, bundle$catalog$result_key)
      counts[index, ids] <- counts[index, ids, drop = FALSE] + 1L
    }
    target <- if (phase == "point") rep(1L, nrow(counts)) else as.integer(bundle$catalog$bootstrap)
    stopifnot(all(counts == target))
  }
  invisible(TRUE)
}
