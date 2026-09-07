#!/usr/bin/env Rscript
# Actions: seal, verify, task, inventory, merge. Run from the repository root.
paper_spectral_command <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- c("action", "run-dir", "task-id", "phase", "resource-class", "array-index", "cores", "fast")
  if (any(!grepl("^--[a-z-]+=.+$", args))) stop("Use --name=value arguments.")
  names_arg <- sub("^--([^=]+)=.*$", "\\1", args)
  if (anyDuplicated(names_arg) || any(!names_arg %in% allowed)) stop("Unknown or duplicate argument.")
  values <- stats::setNames(sub("^--[^=]+=", "", args), names_arg)
  get <- function(name, default = NULL) if (name %in% names(values)) values[[name]] else default
  if (is.null(get("run-dir"))) stop("Supply --run-dir=PATH.")
  root <- normalizePath(getwd(), winslash = "/")
  api <- new.env(parent = baseenv())
  sys.source(file.path(root, "paper/lib/paper_worker.R"), envir = api)
  ctx <- api$paper_worker_load(root)
  action <- match.arg(get("action", "task"), c("seal", "verify", "task", "inventory", "merge"))
  if (action == "seal") return(ctx$api$paper_run_seal(ctx, get("run-dir")))
  state <- ctx$api$paper_run_read(ctx, get("run-dir"), check_runtime = action != "inventory")
  if (action == "verify") {
    if (isTRUE(state$seal$validation_only)) stop("Validation fixtures cannot be submitted as production runs.")
    cat("SEALED SOURCE AND RUNTIME: PASS\n")
    return(invisible(TRUE))
  }
  if (action == "inventory") return(ctx$api$paper_run_inventory(state, validate = get("fast", "false") != "true"))
  if (action == "merge") return(ctx$api$paper_run_merge(ctx, state))
  id <- get("task-id")
  if (is.null(id)) {
    row <- state$tasks[state$tasks$phase == get("phase", "") &
      state$tasks$resource_class == get("resource-class", "") &
      state$tasks$array_index == suppressWarnings(as.numeric(get("array-index", "NA"))), , drop = FALSE]
    if (nrow(row) != 1L || anyNA(row$task_id)) stop("Task selection must identify exactly one row.")
    id <- row$task_id
  }
  cores <- suppressWarnings(as.numeric(get("cores", "1")))
  allocated <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.finite(allocated) && cores > allocated) stop("Requested workers exceed allocated CPUs.")
  ctx$api$paper_run_task(ctx, state, as.numeric(id), cores)
}
if (sys.nframe() == 0L) paper_spectral_command()
