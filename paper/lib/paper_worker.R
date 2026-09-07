# Execution, checkpoints and validated merging for the spectral simulation run.

paper_worker_load <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  api <- new.env(parent = baseenv())
  for (f in c("paper_engine.R", "paper_simulation.R", "paper_external.R", "paper_worker.R"))
    sys.source(file.path(root, "paper/lib", f), envir = api)
  list(api = api, engine = api$paper_engine_load(root), root = root)
}

paper_object_md5 <- function(x) {
  path <- tempfile("mrrr_hash_")
  on.exit(unlink(path), add = TRUE)
  con <- file(path, open = "wb")
  tryCatch(serialize(x, con, version = 2L), finally = close(con))
  unname(tools::md5sum(path))
}

paper_write_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".mrrr_write_", tmpdir = dirname(path))
  on.exit(unlink(temp), add = TRUE)
  saveRDS(x, temp, version = 2L)
  # Retain every replaced incomplete/checkpoint version. A completed checkpoint
  # is validated and skipped by the caller, not overwritten.
  if (file.exists(path)) {
    previous <- tempfile(paste0(basename(path), ".previous_"), tmpdir = dirname(path))
    if (!file.rename(path, previous)) stop("Cannot preserve previous output: ", path)
  }
  if (!file.rename(temp, path)) stop("Cannot install output: ", path)
  invisible(path)
}

paper_runtime <- function(methods) {
  required <- paper_external_require(methods)
  pkgs <- unique(c("stats", "parallel", required))
  if (length(required)) {
    installed <- utils::installed.packages()
    dependency <- tools::package_dependencies(required, db = installed,
      which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
    pkgs <- unique(c(pkgs, unlist(dependency, use.names = FALSE)))
  }
  pkgs <- sort(setdiff(pkgs, "R"))
  versions <- do.call(rbind, lapply(pkgs, function(p) {
    d <- utils::packageDescription(p)
    data.frame(package = p, version = as.character(d$Version),
      remote_sha = if (is.null(d$RemoteSha)) "" else d$RemoteSha,
      built = if (is.null(d$Built)) "" else d$Built, stringsAsFactors = FALSE)
  }))
  code <- do.call(rbind, lapply(required, function(p) {
    folder <- find.package(p)
    files <- list.files(folder, recursive = TRUE, full.names = FALSE)
    files <- files[grepl("^(DESCRIPTION|NAMESPACE)$|^R/.*\\.(rdb|rdx)$|^libs/", files)]
    data.frame(package = p, path = files,
      md5 = unname(tools::md5sum(file.path(folder, files))), stringsAsFactors = FALSE)
  }))
  list(R = R.version.string, platform = R.version$platform,
    BLAS = unname(extSoftVersion()["BLAS"]), packages = versions, package_code = code)
}

paper_run_seal <- function(ctx, run, validation = FALSE) {
  run <- normalizePath(run, winslash = "/", mustWork = TRUE)
  target <- file.path(run, "sealed_run.rds")
  if (file.exists(target)) stop("Run already sealed.")
  bundle <- readRDS(file.path(run, "simulation_bundle.rds"))
  tasks <- utils::read.csv(file.path(run, "tasks.csv"), stringsAsFactors = FALSE)
  ctx$api$paper_sim_validate_tasks(bundle, tasks)
  if (!isTRUE(validation)) {
    canonical <- ctx$api$paper_sim_read_config(ctx$root)
    canonical$replicates <- bundle$replicates
    canonical$bootstrap_size <- bundle$bootstrap_size
    if (!identical(canonical, bundle$config) ||
        !identical(ctx$api$paper_sim_result_catalog(canonical), bundle$catalog) ||
        !identical(ctx$api$paper_sim_table_map(bundle$catalog), bundle$table_map))
      stop("Prepared design/catalog/table map is not canonical.")
  }
  source <- bundle$source_manifest
  extra <- c("paper/lib/paper_external.R", "paper/lib/paper_worker.R",
             "paper/scripts/34_run_spectral_tasks.R", "paper/scripts/35_validate_spectral_workers.R",
             "paper/slurm/run_spectral_task.sbatch", "paper/slurm/submit_spectral_simulations.sh")
  source <- rbind(source, data.frame(path = extra,
    md5 = unname(tools::md5sum(file.path(ctx$root, extra)))))
  if (anyNA(source$md5) || !identical(source$md5,
      unname(tools::md5sum(file.path(ctx$root, source$path)))))
    stop("Prepared bundle source/input hashes do not match this checkout; prepare in this environment.")
  if (any(!grepl("^chunks/[A-Za-z0-9_.-]+\\.rds$", tasks$output_relative_path)))
    stop("Invalid task output path.")
  runtime <- paper_runtime(unique(bundle$catalog$method))
  seal <- list(schema = "spectral-run-1", validation_only = isTRUE(validation),
    preparation_provenance = bundle$provenance,
    bundle_md5 = unname(tools::md5sum(file.path(run, "simulation_bundle.rds"))),
    tasks_md5 = unname(tools::md5sum(file.path(run, "tasks.csv"))),
    sources = source, runtime = runtime, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  seal$run_id <- paper_object_md5(seal)
  paper_write_rds(seal, target)
  utils::write.csv(runtime$packages, file.path(run, "runtime_packages.csv"), row.names = FALSE)
  utils::write.csv(source, file.path(run, "run_sources.csv"), row.names = FALSE)
  cat("SPECTRAL RUN SEALED:", seal$run_id, "\n")
  invisible(seal)
}

paper_run_read <- function(ctx, run, check_runtime = TRUE) {
  run <- normalizePath(run, winslash = "/", mustWork = TRUE)
  seal <- readRDS(file.path(run, "sealed_run.rds"))
  stopifnot(identical(seal$schema, "spectral-run-1"))
  content <- seal; content$run_id <- NULL
  if (!identical(seal$run_id, paper_object_md5(content))) stop("Sealed metadata checksum mismatch.")
  if (!identical(seal$bundle_md5, unname(tools::md5sum(file.path(run, "simulation_bundle.rds")))) ||
      !identical(seal$tasks_md5, unname(tools::md5sum(file.path(run, "tasks.csv")))) ||
      !identical(seal$sources$md5, unname(tools::md5sum(file.path(ctx$root, seal$sources$path)))))
    stop("Sealed run input, task manifest or computation source changed.")
  bundle <- readRDS(file.path(run, "simulation_bundle.rds"))
  if (isTRUE(check_runtime) && !identical(seal$runtime, paper_runtime(unique(bundle$catalog$method))))
    stop("R/package runtime differs from the sealed run.")
  list(run = run, seal = seal, bundle = bundle,
       tasks = utils::read.csv(file.path(run, "tasks.csv"), stringsAsFactors = FALSE))
}

paper_fit_draw <- function(engine, job) {
  rows <- job$catalog
  fits <- stats::setNames(vector("list", nrow(rows)), rows$result_key)
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]; warning <- character(); value <- NULL
    error <- tryCatch({
      value <- withCallingHandlers(paper_sim_with_state(job$states[[row$result_key]], {
        if (row$method %in% c("ivw", "srivw", "mrdag")) {
          paper_external_fit(row$method, job$Y, job$X, job$Sigma_X, job$Sigma_Y, job$config)
        } else {
          method <- c(naive_mr_rr = "naive", mr_rr = "corrected",
                      regularized_mr_rr = "regularized", sparse_mr_rr = "sparse")[[row$method]]
          paper_engine_fit(engine, method, "simulation", job$Y, job$X, row$working_rank,
            job$Sigma_X, job$W, regularization_rate = job$config$regularization_rate,
            sparse_lambda = job$config$sparse_lambda, sparse_threshold = job$config$sparse_threshold,
            sparse_max_iter = job$config$sparse_max_iter, sparse_tol = job$config$sparse_tol,
            sparse_solver = job$config$sparse_solver)
        }
      }), warning = function(w) {
        warning <<- c(warning, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      stopifnot(is.matrix(value$AB), identical(dim(value$AB), c(ncol(job$Y), ncol(job$X))),
                all(is.finite(value$AB)))
      NULL
    }, error = function(e) conditionMessage(e))
    detail <- if (is.null(error) && job$bootstrap_id == 0L)
      value[intersect(names(value), c("A", "B", "B_raw", "AB_raw", "converged", "iter", "dist", "numerical_diagnostics"))]
      else NULL
    fits[[row$result_key]] <- list(ok = is.null(error),
      estimate = if (is.null(error)) value$AB else NULL, details = detail,
      warnings = warning, error = error)
  }
  list(bootstrap_id = job$bootstrap_id, fits = fits)
}

paper_worker_pool <- function(ctx, cores, methods) {
  cores <- paper_sim_integer(cores, "cores", maximum = 256L)
  if (cores == 1L) return(NULL)
  pool <- parallel::makePSOCKcluster(cores, rscript_args = "--vanilla", outfile = "")
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(pool), add = TRUE)
  initialize <- function(root, libraries, methods) {
    .libPaths(libraries)
    api <- new.env(parent = baseenv())
    for (f in c("paper_engine.R", "paper_simulation.R", "paper_external.R", "paper_worker.R"))
      sys.source(file.path(root, "paper/lib", f), envir = api)
    api$paper_external_require(methods)
    assign(".mrrr_worker_api", api, envir = .GlobalEnv)
    assign(".mrrr_worker_engine", api$paper_engine_load(root), envir = .GlobalEnv)
    TRUE
  }
  environment(initialize) <- baseenv()
  parallel::clusterCall(pool, initialize, ctx$root, .libPaths(), methods)
  ready <- TRUE
  pool
}

paper_worker_map <- function(ctx, jobs, pool = NULL) {
  if (is.null(pool)) return(lapply(jobs, function(job) paper_fit_draw(ctx$engine, job)))
  work <- function(job) get(".mrrr_worker_api", envir = .GlobalEnv)$paper_fit_draw(
    get(".mrrr_worker_engine", envir = .GlobalEnv), job)
  environment(work) <- baseenv()
  parallel::parLapply(pool, jobs, work)
}

paper_boot_summary <- function(draws, key, truth) {
  ok <- vapply(draws, function(d) isTRUE(d$fits[[key]]$ok), logical(1))
  issues <- lapply(which(!ok | vapply(draws, function(d) length(d$fits[[key]]$warnings) > 0L, logical(1))),
    function(i) list(draw = i, error = draws[[i]]$fits[[key]]$error, warnings = draws[[i]]$fits[[key]]$warnings))
  # Partial bootstrap draws are retained as failures, never silently used for
  # final SE/coverage. Resume retries the incomplete replicate with the same RNG.
  if (!all(ok)) return(list(complete = FALSE, successful_draws = sum(ok), issues = issues))
  mat <- vapply(draws, function(d) as.vector(d$fits[[key]]$estimate), numeric(length(truth)))
  ci <- t(apply(mat, 1L, stats::quantile, probs = c(.025, .975), names = FALSE, type = 7L))
  list(complete = TRUE, successful_draws = ncol(mat), se = apply(mat, 1L, stats::sd),
    lower = ci[, 1L], upper = ci[, 2L], coverage = as.numeric(truth >= ci[, 1L] & truth <= ci[, 2L]),
    issues = issues)
}

paper_worker_replicate <- function(ctx, state, task, replicate, pool = NULL) {
  b <- state$bundle
  data <- ctx$api$paper_sim_data(b, task$design, task$setting, replicate)
  keys <- strsplit(task$result_keys, ";", fixed = TRUE)[[1L]]
  rows <- b$catalog[match(keys, b$catalog$result_key), ]
  cfg <- b$config[b$config$design == task$design & b$config$setting == task$setting, ]
  p <- b$parameters[[paste(task$design, task$setting, sep = "/")]]
  B <- if (task$phase == "bootstrap") b$bootstrap_size else 0L
  seeds <- stats::setNames(lapply(seq_len(nrow(rows)), function(i)
    ctx$api$paper_sim_method_states(b, task$design, task$setting, replicate,
      rows$method[i], rows$working_rank[i], B = max(1L, B))), keys)
  indices <- if (B) ctx$api$paper_sim_resamples(b, task$design, task$setting, replicate) else NULL
  jobs <- lapply(if (B) seq_len(B) else 0L, function(draw) {
    idx <- if (draw) indices[, draw] else seq_len(nrow(data$X))
    list(bootstrap_id = as.integer(draw), X = data$X[idx, , drop = FALSE], Y = data$Y[idx, , drop = FALSE],
      Sigma_X = p$Sigma_X, Sigma_Y = p$Sigma_Y, W = p$weight.matrix, config = cfg, catalog = rows,
      states = lapply(seeds, function(s) s[, draw + 1L]))
  })
  draws <- paper_worker_map(ctx, jobs, pool)
  fits <- if (!B) draws[[1L]]$fits else stats::setNames(lapply(keys, function(k)
    paper_boot_summary(draws, k, as.vector(p$C))), keys)
  ok <- if (!B) all(vapply(fits, function(f) isTRUE(f$ok), logical(1))) else
    all(vapply(fits, function(f) isTRUE(f$complete), logical(1)))
  list(schema = "spectral-replicate-1", run_id = state$seal$run_id,
    task_id = task$task_id, phase = task$phase, design = task$design, setting = task$setting,
    replicate = as.integer(replicate), data_md5 = paper_object_md5(list(X = data$X, Y = data$Y)),
    resample_md5 = if (B) paper_object_md5(indices) else NA_character_,
    bootstrap_size = B, result_keys = keys, complete = ok, fits = fits)
}

paper_validate_replicate <- function(record, state, task, id, require_complete = TRUE) {
  keys <- strsplit(task$result_keys, ";", fixed = TRUE)[[1L]]
  stopifnot(identical(record$schema, "spectral-replicate-1"),
    identical(record$run_id, state$seal$run_id), identical(record$task_id, task$task_id),
    identical(record$replicate, as.integer(id)), identical(record$phase, task$phase),
    identical(record$design, task$design), identical(record$setting, task$setting),
    identical(record$result_keys, keys), identical(names(record$fits), keys),
    is.character(record$data_md5), length(record$data_md5) == 1L, nchar(record$data_md5) == 32L,
    identical(record$bootstrap_size, task$bootstrap_size))
  if (require_complete && !isTRUE(record$complete)) stop("Incomplete replicate ", id)
  if (isTRUE(record$complete)) for (key in keys) {
    fit <- record$fits[[key]]
    if (task$phase == "point") {
      stopifnot(isTRUE(fit$ok), identical(dim(fit$estimate), c(3L, 9L)), all(is.finite(fit$estimate)))
    } else {
      stopifnot(isTRUE(fit$complete), fit$successful_draws == task$bootstrap_size,
        is.character(record$resample_md5), length(record$resample_md5) == 1L, nchar(record$resample_md5) == 32L)
      for (field in c("se", "lower", "upper", "coverage"))
        stopifnot(is.numeric(fit[[field]]), length(fit[[field]]) == 27L, all(is.finite(fit[[field]])))
      truth <- as.vector(state$bundle$truths[[task$design]]$C)
      stopifnot(all(fit$se >= 0), all(fit$lower <= fit$upper),
        identical(fit$coverage, as.numeric(truth >= fit$lower & truth <= fit$upper)))
    }
  }
  invisible(TRUE)
}

paper_validate_chunk <- function(chunk, state, task) {
  ids <- seq.int(task$replicate_start, task$replicate_end)
  stopifnot(identical(chunk$schema, "spectral-chunk-1"), identical(chunk$run_id, state$seal$run_id),
    identical(chunk$task, task), identical(names(chunk$replicates), as.character(ids)))
  for (id in ids) paper_validate_replicate(chunk$replicates[[as.character(id)]], state, task, id)
  invisible(TRUE)
}

paper_run_task <- function(ctx, state, task_id, cores = 1L) {
  task_id <- paper_sim_integer(task_id, "task id", maximum = nrow(state$tasks))
  task <- state$tasks[state$tasks$task_id == task_id, , drop = FALSE]
  path <- file.path(state$run, task$output_relative_path)
  if (file.exists(path)) {
    chunk <- readRDS(path)
    paper_validate_chunk(chunk, state, task)
    cat("SKIP validated completed task:", task_id, "\n")
    return(invisible(chunk))
  }
  ids <- seq.int(task$replicate_start, task$replicate_end)
  folder <- file.path(state$run, "checkpoints", paste0("task-", task_id))
  records <- stats::setNames(vector("list", length(ids)), as.character(ids))
  methods <- state$bundle$catalog$method[match(strsplit(task$result_keys, ";", fixed = TRUE)[[1L]],
                                               state$bundle$catalog$result_key)]
  # Bootstrap parallelism uses a persistent socket pool; point jobs request one
  # CPU and avoid nested method parallelism. Explicit seeds make pool size irrelevant.
  cores <- paper_sim_integer(cores, "cores", maximum = 256L)
  effective <- if (task$phase == "bootstrap") min(cores, task$bootstrap_size) else 1L
  pool <- paper_worker_pool(ctx, effective, methods)
  if (!is.null(pool)) on.exit(parallel::stopCluster(pool), add = TRUE)
  for (id in ids) {
    saved <- file.path(folder, sprintf("rep-%04d.rds", id))
    record <- if (file.exists(saved)) readRDS(saved) else NULL
    if (!is.null(record)) paper_validate_replicate(record, state, task, id, require_complete = FALSE)
    if (is.null(record) || !isTRUE(record$complete)) {
      record <- paper_worker_replicate(ctx, state, task, id, pool)
      paper_write_rds(record, saved)
    }
    records[[as.character(id)]] <- record
    cat("Task", task_id, "replicate", id, if (isTRUE(record$complete)) "PASS" else "FAILED", "\n")
    utils::flush.console()
  }
  failed <- !vapply(records, function(r) isTRUE(r$complete), logical(1))
  if (any(failed)) {
    paper_write_rds(records[failed], file.path(folder, "failed_replicates.rds"))
    errors <- list()
    for (record in records[failed]) for (key in names(record$fits)) {
      fit <- record$fits[[key]]
      issues <- if (task$phase == "point" && !isTRUE(fit$ok))
        list(list(draw = 0L, error = fit$error)) else if (task$phase == "bootstrap") fit$issues else list()
      for (issue in issues) if (!is.null(issue$error)) errors[[length(errors) + 1L]] <-
        data.frame(replicate = record$replicate, result_key = key,
          bootstrap_draw = issue$draw, message = issue$error, stringsAsFactors = FALSE)
    }
    if (length(errors)) {
      errors <- do.call(rbind, errors)
      utils::write.csv(errors, file.path(folder, "failed_fits.csv"), row.names = FALSE)
      print(utils::head(errors, 20L), row.names = FALSE)
    }
    stop("Task has failed fits in replicate(s): ", paste(names(records)[failed], collapse = ", "),
         ". Checkpoints retained; no complete chunk written.")
  }
  paper_engine_check_sources(ctx$engine)
  if (!identical(state$seal$sources$md5,
      unname(tools::md5sum(file.path(ctx$root, state$seal$sources$path))))) stop("Source changed during task.")
  chunk <- list(schema = "spectral-chunk-1", run_id = state$seal$run_id, task = task, replicates = records)
  paper_validate_chunk(chunk, state, task)
  paper_write_rds(chunk, path)
  cat("SPECTRAL TASK: PASS", task_id, "\n")
  invisible(chunk)
}

paper_run_inventory <- function(state, validate = TRUE, write = TRUE) {
  out <- state$tasks
  out$status <- "missing"
  out$message <- ""
  for (i in seq_len(nrow(out))) {
    path <- file.path(state$run, out$output_relative_path[i])
    if (!file.exists(path)) next
    error <- if (isTRUE(validate)) tryCatch({
      paper_validate_chunk(readRDS(path), state, state$tasks[i, , drop = FALSE]); NULL
    }, error = function(e) conditionMessage(e)) else NULL
    out$status[i] <- if (is.null(error)) if (validate) "valid" else "present" else "invalid"
    if (!is.null(error)) out$message[i] <- error
  }
  if (isTRUE(write)) {
    utils::write.csv(out, file.path(state$run, "inventory.csv"), row.names = FALSE)
    utils::write.csv(out[out$status == "missing", names(state$tasks)],
                     file.path(state$run, "pending_tasks.csv"), row.names = FALSE)
    groups <- expand.grid(phase = c("point", "bootstrap"), resource_class = c("standard", "mrdag"),
                          stringsAsFactors = FALSE)
    groups$indices <- vapply(seq_len(nrow(groups)), function(i) {
      selected <- out$array_index[out$status == "missing" & out$phase == groups$phase[i] &
                                    out$resource_class == groups$resource_class[i]]
      if (length(selected)) paste(selected, collapse = ",") else "-"
    }, character(1))
    utils::write.table(groups, file.path(state$run, "pending_arrays.tsv"),
                       sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  print(table(out$phase, out$resource_class, out$status))
  if (any(out$status == "invalid")) stop("Invalid chunks detected; inspect inventory.csv. Existing files retained.")
  invisible(out)
}

paper_run_merge <- function(ctx, state) {
  inventory <- paper_run_inventory(state)
  if (any(inventory$status != "valid")) stop("Merge requires every expected chunk; see pending_tasks.csv.")
  b <- state$bundle; n <- b$replicates
  if (n < 2L) stop("At least two replicates are needed for empirical SD.")
  output <- file.path(state$run, "merged")
  if (dir.exists(output)) stop("Merged output exists; use the retained results or a new run directory.")
  blank <- function() matrix(NA_real_, 27L, n)
  results <- stats::setNames(lapply(seq_len(nrow(b$catalog)), function(i) {
    r <- list(point = blank(), point_diagnostics = vector("list", n))
    if (b$catalog$bootstrap[i]) {
      for (name in c("se", "lower", "upper", "coverage")) r[[name]] <- blank()
      r$bootstrap_issues <- vector("list", n)
    }
    r
  }), b$catalog$result_key)
  cells <- paste(b$config$design, b$config$setting, sep = "/")
  data_hash <- stats::setNames(lapply(cells, function(x) rep(NA_character_, n)), cells)
  resample_hash <- data_hash
  for (i in seq_len(nrow(state$tasks))) {
    task <- state$tasks[i, , drop = FALSE]
    chunk <- readRDS(file.path(state$run, task$output_relative_path))
    # Inventory already validated every chunk, but validate again on reading to
    # reject replacement between inventory and merge.
    paper_validate_chunk(chunk, state, task)
    cell <- paste(task$design, task$setting, sep = "/")
    for (r in chunk$replicates) {
      id <- r$replicate
      old <- data_hash[[cell]][id]
      if (!is.na(old) && !identical(old, r$data_md5)) stop("Point/bootstrap or resource classes used different data.")
      data_hash[[cell]][id] <- r$data_md5
      if (task$phase == "bootstrap") {
        old <- resample_hash[[cell]][id]
        if (!is.na(old) && !identical(old, r$resample_md5)) stop("Resource classes used different SNP resamples.")
        resample_hash[[cell]][id] <- r$resample_md5
      }
      for (key in r$result_keys) {
        fit <- r$fits[[key]]
        if (task$phase == "point") {
          if (any(is.finite(results[[key]]$point[, id]))) stop("Overlapping point results.")
          results[[key]]$point[, id] <- as.vector(fit$estimate)
          results[[key]]$point_diagnostics[[id]] <- list(details = fit$details, warnings = fit$warnings)
        } else {
          if (any(is.finite(results[[key]]$se[, id]))) stop("Overlapping bootstrap results.")
          for (name in c("se", "lower", "upper", "coverage")) results[[key]][[name]][, id] <- fit[[name]]
          results[[key]]$bootstrap_issues[[id]] <- fit$issues
        }
      }
    }
  }
  summaries <- list(); entrywise <- list()
  for (i in seq_len(nrow(b$catalog))) {
    row <- b$catalog[i, ]; key <- row$result_key; r <- results[[key]]
    truth <- as.vector(b$truths[[row$design]]$C)
    stopifnot(all(is.finite(r$point)))
    e <- data.frame(result_key = key, outcome = rep(1:3, 9), exposure = rep(1:9, each = 3),
      bias = abs(rowMeans(r$point) - truth), sd = apply(r$point, 1L, stats::sd),
      rmse = sqrt(rowMeans((r$point - truth)^2)), se = NA_real_, cp_percent = NA_real_)
    if (row$bootstrap) {
      for (name in c("se", "lower", "upper", "coverage")) stopifnot(all(is.finite(r[[name]])))
      e$se <- rowMeans(r$se); e$cp_percent <- 100 * rowMeans(r$coverage)
    }
    summary <- row
    for (metric in c("bias", "sd", "se", "cp_percent", "rmse")) {
      values <- if (all(is.na(e[[metric]]))) rep(NA_real_, 3) else
        stats::quantile(e[[metric]], probs = c(.5, .25, .75), names = FALSE, type = 7L)
      for (j in 1:3) summary[[paste0(metric, c("_median", "_q25", "_q75")[j])]] <- values[j]
    }
    summary$nonconverged_point_fits <- sum(vapply(r$point_diagnostics,
      function(d) identical(d$details$converged, FALSE), logical(1)))
    summary$projected_point_fits <- sum(vapply(r$point_diagnostics,
      function(d) isTRUE(d$details$numerical_diagnostics$corrected_covariance_projected), logical(1)))
    summary$point_warning_count <- sum(vapply(r$point_diagnostics, function(d) length(d$warnings), integer(1)))
    summary$bootstrap_draws_with_issues <- if (row$bootstrap)
      sum(vapply(r$bootstrap_issues, length, integer(1))) else 0L
    entrywise[[i]] <- e; summaries[[i]] <- summary
  }
  summary <- do.call(rbind, summaries)
  table_rows <- cbind(b$table_map[, c("table_id", "point_result_key", "bootstrap_result_key")],
                      summary[match(b$table_map$result_key, summary$result_key), ])
  shared <- intersect(table_rows$result_key[table_rows$table_id == "main_generic"],
                      table_rows$result_key[table_rows$table_id == "rank_misspecification"])
  for (key in shared) {
    selected <- table_rows[table_rows$result_key == key, names(summary), drop = FALSE]
    stopifnot(nrow(selected) == 2L, identical(unname(as.list(selected[1L, ])), unname(as.list(selected[2L, ]))))
  }
  paper_engine_check_sources(ctx$engine)
  if (!identical(state$seal$sources$md5,
      unname(tools::md5sum(file.path(ctx$root, state$seal$sources$path))))) stop("Source changed during merge.")
  merged <- list(schema = "spectral-merged-1", run_id = state$seal$run_id,
    seal = state$seal, config = b$config, truths = b$truths, parameters = b$parameters,
    prediction_exposure = b$prediction_exposure, catalog = b$catalog, table_map = b$table_map,
    results = results, data_md5 = data_hash, resample_md5 = resample_hash, summary = summary, table_rows = table_rows)
  staging <- tempfile("merge_build_", tmpdir = state$run)
  dir.create(staging)
  paper_write_rds(merged, file.path(staging, "spectral_simulation_results.rds"))
  utils::write.csv(summary, file.path(staging, "summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, entrywise), file.path(staging, "entrywise.csv"), row.names = FALSE)
  utils::write.csv(table_rows, file.path(staging, "table_rows.csv"), row.names = FALSE)
  writeLines(c("SPECTRAL SIMULATION MERGE: PASS", "All expected fits and bootstrap draws present and finite.",
    "Point/bootstrap datasets agree; SNP resamples agree among the present resource classes.",
    paste("Shared main/rank-sensitivity rows identical:", length(shared)),
    paste("Sparse point fits reported as nonconverged:", sum(summary$nonconverged_point_fits)),
    paste("Designs:", paste(unique(b$config$design), collapse = ", ")),
    paste("Validation fixture only:", state$seal$validation_only)),
    file.path(staging, "STATUS.txt"))
  if (dir.exists(output) || !file.rename(staging, output))
    stop("Could not finalize merge; prepared output retained at ", staging)
  cat(readLines(file.path(output, "STATUS.txt")), sep = "\n")
  invisible(merged)
}
