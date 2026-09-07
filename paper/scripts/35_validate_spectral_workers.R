#!/usr/bin/env Rscript
# Small execution/merge validation, including native external comparators.
paper_validate_spectral_workers <- function(root = getwd(), output = NULL, base_only = FALSE, cores = 2L) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (is.null(output)) output <- file.path(root, "paper/output/spectral_rebuild/worker_validation")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  writeLines("SPECTRAL WORKER VALIDATION: RUNNING", file.path(output, "STATUS.txt"))
  api <- new.env(parent = baseenv())
  sys.source(file.path(root, "paper/lib/paper_worker.R"), envir = api)
  ctx <- api$paper_worker_load(root); a <- ctx$api
  cores <- a$paper_sim_integer(cores, "validation cores", maximum = 4L)
  if (!base_only && cores < 2L) stop("Full native validation requires at least two workers.")
  bundle <- a$paper_sim_build_bundle(root, ctx$engine, replicates = 2L, bootstrap_size = 3L)
  bundle$config <- bundle$config[bundle$config$design == "generic" & bundle$config$setting == 4L, , drop = FALSE]
  bundle$catalog <- bundle$catalog[bundle$catalog$design == "generic" & bundle$catalog$setting == 4L, , drop = FALSE]
  if (base_only) bundle$catalog <- bundle$catalog[bundle$catalog$method %in%
    c("naive_mr_rr", "mr_rr", "regularized_mr_rr"), , drop = FALSE]
  rownames(bundle$catalog) <- NULL; rownames(bundle$config) <- NULL
  bundle$table_map <- a$paper_sim_table_map(bundle$catalog)
  bundle$parameters <- bundle$parameters["generic/4"]
  tasks <- a$paper_sim_tasks(bundle, 1L, 1L, 1L, 1L)
  run <- tempfile("run_", tmpdir = output); dir.create(run)
  saveRDS(bundle, file.path(run, "simulation_bundle.rds"), version = 2L)
  utils::write.csv(tasks, file.path(run, "tasks.csv"), row.names = FALSE)
  a$paper_run_seal(ctx, run, validation = TRUE)
  state <- a$paper_run_read(ctx, run)
  rows <- list()
  check <- function(label, expression) {
    error <- tryCatch({force(expression); NULL}, error = function(e) conditionMessage(e))
    rows[[length(rows) + 1L]] <<- data.frame(check = label,
      status = if (is.null(error)) "PASS" else "FAIL", detail = if (is.null(error)) "" else error)
    cat(if (is.null(error)) "PASS: " else "FAIL: ", label,
        if (is.null(error)) "" else paste0(" - ", error), "\n", sep = "")
  }
  throws <- function(expression) stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
  check("Comparator SE columns match their respective traits", {
    se <- a$paper_external_se(diag(c(1, 4, 9)), 5L)
    stopifnot(identical(se, matrix(rep(c(1, 2, 3), each = 5), 5, 3)))
    old <- matrix(rep(c(1, 2, 3), each = 5), nrow = 5, byrow = TRUE)
    stopifnot(!identical(se, old), all(se[, 1] == 1), all(se[, 2] == 2), all(se[, 3] == 3))
  })
  check("Bootstrap SD and percentile intervals use all draws", {
    draws <- lapply(1:3, function(b) list(fits = list(k = list(ok = TRUE,
      estimate = matrix(as.numeric(b), 3, 9), warnings = character(), error = NULL))))
    summary <- a$paper_boot_summary(draws, "k", rep(2, 27))
    stopifnot(summary$complete, all(summary$se == 1), all(abs(summary$lower - 1.05) < 1e-12),
              all(abs(summary$upper - 2.95) < 1e-12), all(summary$coverage == 1))
    draws[[2]]$fits$k$ok <- FALSE; draws[[2]]$fits$k$error <- "fixture failure"
    partial <- a$paper_boot_summary(draws, "k", rep(2, 27))
    stopifnot(!partial$complete, partial$successful_draws == 2L, is.null(partial$se))
  })
  check("Incomplete run cannot be merged", {
    throws(a$paper_run_merge(ctx, state))
    stopifnot(!dir.exists(file.path(run, "merged")))
  })
  check("Point and bootstrap tasks execute with the recorded method streams", {
    for (id in tasks$task_id) a$paper_run_task(ctx, state, id, cores = 1L)
    inventory <- a$paper_run_inventory(state)
    stopifnot(all(inventory$status == "valid"))
  })
  check("Completed tasks are reused without changing saved results", {
    paths <- file.path(run, tasks$output_relative_path)
    before <- unname(tools::md5sum(paths))
    for (id in tasks$task_id) a$paper_run_task(ctx, state, id, cores = 1L)
    stopifnot(identical(before, unname(tools::md5sum(paths))))
  })
  check("Missing chunk is reconstructed from completed checkpoints", {
    task <- state$tasks[1L, , drop = FALSE]
    path <- file.path(run, task$output_relative_path)
    original <- readRDS(path)
    checkpoint <- file.path(run, "checkpoints", paste0("task-", task$task_id), "rep-0001.rds")
    before <- unname(tools::md5sum(checkpoint))
    stopifnot(file.rename(path, paste0(path, ".validation_backup")))
    recovered <- a$paper_run_task(ctx, state, task$task_id, cores = 1L)
    stopifnot(identical(original, recovered), identical(before, unname(tools::md5sum(checkpoint))))
  })
  if (cores >= 2L) check("Serial and socket-parallel bootstrap results are identical", {
    for (resource in unique(tasks$resource_class)) {
      task <- state$tasks[state$tasks$phase == "bootstrap" & state$tasks$resource_class == resource, , drop = FALSE][1L, , drop = FALSE]
      expected <- readRDS(file.path(run, task$output_relative_path))$replicates[[1L]]
      keys <- strsplit(task$result_keys, ";", fixed = TRUE)[[1L]]
      methods <- bundle$catalog$method[match(keys, bundle$catalog$result_key)]
      pool <- a$paper_worker_pool(ctx, cores, methods)
      actual <- tryCatch(a$paper_worker_replicate(ctx, state, task, 1L, pool),
                         finally = parallel::stopCluster(pool))
      stopifnot(identical(expected, actual))
    }
  })
  check("Runtime/source/task identity mismatches are rejected", {
    task <- state$tasks[1L, , drop = FALSE]
    chunk <- readRDS(file.path(run, task$output_relative_path))
    chunk$run_id <- "other-run"
    throws(a$paper_validate_chunk(chunk, state, task))
    file <- file.path(run, "tasks.csv")
    bytes <- readBin(file, what = "raw", n = file.info(file)$size)
    tryCatch({
      cat("\n", file = file, append = TRUE)
      throws(a$paper_run_read(ctx, run))
    }, finally = writeBin(bytes, file))
  })
  check("Non-finite estimates cannot pass chunk validation", {
    task <- state$tasks[1L, , drop = FALSE]
    chunk <- readRDS(file.path(run, task$output_relative_path))
    chunk$replicates[[1L]]$fits[[1L]]$estimate[1L, 1L] <- Inf
    throws(a$paper_validate_chunk(chunk, state, task))
  })
  if (!base_only) check("IVW adapter agrees with constant-weight least squares", {
    data <- a$paper_sim_data(bundle, "generic", 4L, 1L)
    key <- bundle$catalog$result_key[bundle$catalog$method == "ivw"]
    task <- state$tasks[state$tasks$phase == "point" & state$tasks$resource_class == "standard", , drop = FALSE][1L, , drop = FALSE]
    actual <- readRDS(file.path(run, task$output_relative_path))$replicates[[1L]]$fits[[key]]$estimate
    expected <- t(solve(crossprod(data$X), crossprod(data$X, data$Y)))
    stopifnot(max(abs(actual - expected)) < 1e-8)
  })
  check("Merge validates common data/resamples and identical shared table rows", {
    merged <- a$paper_run_merge(ctx, state)
    stopifnot(nrow(merged$summary) == nrow(bundle$catalog),
      identical(readRDS(file.path(run, "merged/spectral_simulation_results.rds")), merged))
    key <- bundle$catalog$result_key[bundle$catalog$method == "regularized_mr_rr" & bundle$catalog$working_rank == 2L]
    selected <- merged$table_rows[merged$table_rows$result_key == key, ]
    stopifnot(nrow(selected) == 2L, length(unique(selected$bias_median)) == 1L,
              length(unique(selected$se_median)) == 1L, length(unique(selected$cp_percent_median)) == 1L)
  })
  table <- do.call(rbind, rows)
  utils::write.csv(table, file.path(output, "checks.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
  mode <- if (base_only) "PARTIAL_PASS (external comparators and Sparse excluded)" else "PASS"
  ok <- all(table$status == "PASS")
  status <- paste("SPECTRAL WORKER VALIDATION:", if (ok) mode else "FAIL")
  writeLines(c(status, paste("Socket workers tested:", cores >= 2L), paste("Validation run:", run)),
             file.path(output, "STATUS.txt"))
  if (!ok) stop("Worker validation failed; inspect checks.csv under ", output)
  cat(status, "\nReport:", normalizePath(output, winslash = "/"), "\n")
  invisible(table)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  base <- "--base-only" %in% args
  values <- args[args != "--base-only"]
  if (any(!grepl("^--(output|cores)=.+$", values))) stop("Unknown argument.")
  get <- function(key, default = NULL) {
    a <- values[startsWith(values, paste0("--", key, "="))]
    if (length(a) > 1L) stop("Duplicate argument.")
    if (length(a)) sub(paste0("^--", key, "="), "", a) else default
  }
  paper_validate_spectral_workers(output = get("output"), base_only = base,
                                  cores = as.numeric(get("cores", "2")))
}
