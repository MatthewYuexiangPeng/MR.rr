#!/usr/bin/env Rscript
# Build the data/configuration bundle and computation manifests; submit no jobs.
paper_prepare_spectral_simulations <- function(root = getwd(), output = NULL,
                                               replicates = NULL, bootstrap_size = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (is.null(output)) output <- file.path(root, "paper/output/spectral_rebuild/simulations_v1")
  if (dir.exists(output) || file.exists(output)) stop("Output exists; use a new --output directory.")
  api <- new.env(parent = baseenv())
  for (name in c("paper_engine.R", "paper_simulation.R"))
    sys.source(file.path(root, "paper/lib", name), envir = api)
  engine <- api$paper_engine_load(root)
  bundle <- api$paper_sim_build_bundle(root, engine, replicates, bootstrap_size)
  tasks <- api$paper_sim_tasks(bundle)
  api$paper_sim_validate_tasks(bundle, tasks)
  api$paper_engine_check_sources(engine)
  git_commit <- if (nzchar(Sys.which("git"))) tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
                              stdout = TRUE, stderr = FALSE)), error = function(e) "unavailable") else "unavailable"
  bundle$provenance <- list(git_commit = git_commit, R_version = R.version.string,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    mode = if (bundle$replicates == 1000L && bundle$bootstrap_size == 300L) "full_size" else "development")
  bundle$source_manifest <- rbind(bundle$source_manifest,
    data.frame(path = "paper/scripts/32_prepare_spectral_simulations.R",
      md5 = unname(tools::md5sum(file.path(root, "paper/scripts/32_prepare_spectral_simulations.R")))))
  dir.create(output, recursive = TRUE)
  saveRDS(bundle, file.path(output, "simulation_bundle.rds"), version = 2)
  for (name in c("config", "catalog", "table_map", "source_manifest"))
    utils::write.csv(bundle[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  utils::write.csv(tasks, file.path(output, "tasks.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "session_info.txt"))
  writeLines(c("SIMULATION DESIGN AND MANIFEST: PASS", "Scope: three simulation designs and rank sensitivity.",
    "No estimation workers or Slurm jobs were run.",
    "Rank selection, standalone support recovery, tuning paths and real-data tasks are added in later milestones."),
    file.path(output, "STATUS.txt"))
  print(as.data.frame(table(tasks$phase, tasks$resource_class)))
  cat("Unique result keys:", nrow(bundle$catalog), "\nTable rows:", nrow(bundle$table_map),
      "\nTasks:", nrow(tasks), "\nReplicates:", bundle$replicates,
      "\nBootstrap draws:", bundle$bootstrap_size, "\nSIMULATION DESIGN AND MANIFEST: PASS\n",
      "Output:", normalizePath(output, winslash = "/"), "\n")
  invisible(bundle)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(!grepl("^--(repo|output|replicates|bootstrap-size)=.+$", args))) stop("Unknown argument.")
  value <- function(key, default = NULL) {
    selected <- args[startsWith(args, paste0("--", key, "="))]
    if (length(selected) > 1L) stop("Duplicate argument: ", key)
    if (length(selected)) sub(paste0("^--", key, "="), "", selected) else default
  }
  numeric_value <- function(key) {
    x <- value(key)
    if (is.null(x)) NULL else suppressWarnings(as.numeric(x))
  }
  paper_prepare_spectral_simulations(value("repo", getwd()), value("output"),
                                      numeric_value("replicates"), numeric_value("bootstrap-size"))
}
