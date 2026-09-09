#!/usr/bin/env Rscript
# CLI for the remaining-paper computation. This script never reads frozen fits.
remaining_main <- function(args = commandArgs(TRUE)) {
  opts <- list(action = "", root = getwd(), run = "", profile = "full",
    reference = "", task = "", cores = "1")
  for (a in args) {
    if (!grepl("^--[a-z-]+=", a)) stop("Use --name=value arguments: ", a)
    name <- sub("^--([^=]+)=.*$", "\\1", a)
    if (!name %in% names(opts)) stop("Unknown option: ", name)
    opts[[name]] <- sub("^[^=]+=", "", a)
  }
  if (!opts$action %in% c("prepare", "verify", "task", "inventory", "merge"))
    stop("--action must be prepare, verify, task, inventory, or merge.")
  if (!nzchar(opts$run)) stop("Supply --run=PATH.")
  a <- new.env(parent = baseenv())
  sys.source(file.path(opts$root, "paper/lib/paper_remaining.R"), a)
  ctx <- a$rem_load(opts$root)
  if (opts$action == "prepare") return(invisible(ctx$api$rem_prepare(ctx, opts$run, opts$profile, opts$reference)))
  state <- ctx$api$rem_open(ctx, opts$run, runtime = opts$action != "inventory")
  if (opts$action == "verify") cat("SPECTRAL REMAINING RUN VERIFICATION: PASS\n")
  if (opts$action == "task") {
    task <- suppressWarnings(as.integer(opts$task)); cores <- suppressWarnings(as.integer(opts$cores))
    stopifnot(length(task) == 1L, !is.na(task), length(cores) == 1L, !is.na(cores), cores >= 1L)
    ctx$api$rem_run_task(ctx, state, task, cores)
  }
  if (opts$action == "inventory") {
    inv <- ctx$api$rem_inventory(state); print(table(inv$kind, inv$resource, inv$status))
    if (any(inv$status == "invalid")) stop("Invalid saved chunks require inspection; they will not be silently replaced.")
    groups <- list(rank = inv$kind == "rank", eta = inv$kind %in% c("sim_eta", "real_eta"),
      bootstrap_standard = inv$kind == "real_bootstrap" & inv$resource == "standard",
      bootstrap_mrdag = inv$kind == "real_bootstrap" & inv$resource == "mrdag")
    lines <- vapply(names(groups), function(g) {
      ids <- inv$task_id[groups[[g]] & inv$status == "missing"]
      paste(g, if (length(ids)) paste(ids, collapse = ",") else "-", sep = "\t")
    }, character(1))
    writeLines(lines, file.path(state$run, "pending_arrays.tsv"))
  }
  if (opts$action == "merge") ctx$api$rem_merge(ctx, state)
  invisible(TRUE)
}
if (sys.nframe() == 0L) remaining_main()
