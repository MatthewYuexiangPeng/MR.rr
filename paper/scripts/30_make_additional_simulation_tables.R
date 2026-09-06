#!/usr/bin/env Rscript
# Generate manuscript tables from merged additional-simulation results.
# Base R only. Does not run estimators, simulate data, or reconstruct the truth.

additional_tables_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

additional_tables_close <- function(x, y, label, tolerance = 1e-10) {
  additional_tables_assert(length(x) == length(y) && all(is.finite(x)) &&
    all(is.finite(y)), paste(label, "has missing or nonfinite values."))
  difference <- max(abs(x - y))
  additional_tables_assert(difference <= tolerance,
    paste(label, "differs by", format(difference, scientific = TRUE)))
  invisible(difference)
}

additional_tables_options <- function(args) {
  opts <- list(
    `input-dir` = "paper/output/additional_simulations/cluster_full_b300_delta01",
    `output-dir` = "paper/output/additional_simulations/manuscript_tables",
    overwrite = "false"
  )
  if (identical(args, "--help")) {
    cat("Usage: Rscript --vanilla paper/scripts/30_make_additional_simulation_tables.R\n",
        "  [--input-dir=RUN_DIRECTORY] [--output-dir=TABLE_DIRECTORY]\n",
        "  [--overwrite=true|false]\n", sep = "")
    return(NULL)
  }
  seen <- character()
  for (arg in args) {
    additional_tables_assert(grepl("^--[^=]+=.+$", arg),
      paste("Use --name=value syntax:", arg))
    name <- sub("^--([^=]+)=.*$", "\\1", arg)
    additional_tables_assert(name %in% names(opts) && !name %in% seen,
      paste("Unknown or duplicate option:", name))
    opts[[name]] <- sub("^--[^=]+=", "", arg)
    seen <- c(seen, name)
  }
  additional_tables_assert(opts$overwrite %in% c("true", "false"),
    "--overwrite must be true or false.")
  opts$overwrite <- identical(opts$overwrite, "true")
  opts
}

additional_tables_matrix <- function(x, n, label) {
  additional_tables_assert(is.matrix(x) && identical(dim(x), c(27L, as.integer(n))) &&
    all(is.finite(x)), paste(label, "must be a finite 27-by-N matrix."))
}

additional_tables_group_keys <- function(groups) {
  keys <- vapply(groups, function(g) paste(g$analysis, g$phase, g$setting_index,
    sep = "__"), character(1))
  additional_tables_assert(!anyDuplicated(keys), "Duplicate merged groups.")
  keys
}

additional_tables_quantiles <- function(x) {
  c(med = stats::median(x), q1 = unname(stats::quantile(x, 0.25, type = 7)),
    q3 = unname(stats::quantile(x, 0.75, type = 7)))
}

additional_tables_summarize <- function(point, bootstrap) {
  labels <- c(ivw = "IVW", srivw = "SRIVW", naive_mr_rr = "Naive MR-rr",
    mr_rr = "MR-rr", regularized_mr_rr = "Reg. MR-rr",
    sparse_mr_rr = "Sparse MR-rr", mrdag = "MrDAG")
  n <- point$metadata$total_replicates_per_setting
  B <- bootstrap$metadata$bootstrap_size
  additional_tables_assert(identical(as.integer(n), 1000L) &&
    identical(as.integer(bootstrap$metadata$total_replicates_per_setting), 1000L) &&
    identical(as.integer(B), 300L), "Final tables require N=1000 and B=300.")
  delta <- point$metadata$approximate_delta
  additional_tables_assert(length(delta) == 1L && is.finite(delta) &&
    delta > 0 && delta < 1, "Invalid approximate-low-rank delta.")
  additional_tables_close(delta, bootstrap$metadata$approximate_delta, "Delta")
  for (object in list(point, bootstrap)) {
    additional_tables_assert(identical(object$metadata$approximate_configuration_status,
      "locked"), "Final tables require a locked configuration.")
    inv <- object$inventory
    additional_tables_assert(is.data.frame(inv) && nrow(inv) > 0L &&
      all(c("present", "validated", "complete") %in% names(inv)) &&
      all(inv$present) && all(inv$validated) && all(inv$complete),
      "A merged inventory contains incomplete or unvalidated chunks.")
  }
  expected_groups <- c(paste("rank_misspecification", "standard", 1:4, sep = "__"),
    paste("approximate_low_rank", "standard", 1:4, sep = "__"),
    paste("approximate_low_rank", "mrdag", 1:4, sep = "__"))
  pk <- additional_tables_group_keys(point$groups)
  bk <- additional_tables_group_keys(bootstrap$groups)
  additional_tables_assert(setequal(pk, expected_groups) && setequal(bk, expected_groups),
    "The merged objects do not contain the expected twelve groups.")
  point$groups <- stats::setNames(point$groups, pk)
  bootstrap$groups <- stats::setNames(bootstrap$groups, bk)
  rows <- entries <- list()
  index <- 0L
  largest_bias_difference <- largest_coverage_difference <- 0

  for (group_key in expected_groups) {
    g <- point$groups[[group_key]]
    b <- bootstrap$groups[[group_key]]
    additional_tables_assert(identical(as.integer(g$replicate_ids), seq_len(n)) &&
      identical(as.integer(b$replicate_ids), seq_len(n)), "Replicate IDs are incomplete.")
    additional_tables_assert(identical(as.integer(g$data_seeds), as.integer(b$data_seeds)),
      paste("Point/bootstrap data seeds differ:", group_key))
    additional_tables_close(g$truth$C, b$truth$C, paste(group_key, "truth"))
    additional_tables_assert(is.matrix(g$truth$C) && identical(dim(g$truth$C), c(3L, 9L)),
      "The causal-effect truth must be 3-by-9.")
    additional_tables_assert(identical(g$configuration, b$configuration) &&
      all(g$configuration$status == "locked") &&
      all(g$configuration$monte_carlo_replicates == n) &&
      all(g$configuration$bootstrap_size == B), "Point/bootstrap configurations differ.")
    cfg <- g$configuration
    s <- g$setting_index
    expected_delta <- if (g$analysis == "rank_misspecification") 0 else delta
    additional_tables_close(unique(cfg$third_singular_value), expected_delta, "Configuration delta")
    additional_tables_close(unique(cfg$measurement_error_weight), c(2.5, 1, 2.5, 1)[s], "ME multiplier")
    additional_tables_close(unique(cfg$genetic_effect_weight), c(0.25, 0.25, 1, 1)[s], "Genetic multiplier")
    additional_tables_close(svd(g$truth$C, nu = 0L, nv = 0L)$d,
      c(1, 1, expected_delta), "Truth singular values")
    expected_keys <- if (g$analysis == "rank_misspecification") {
      unlist(lapply(1:3, function(r) paste0("r", r, "__",
        c("regularized_mr_rr", "sparse_mr_rr"))))
    } else if (g$phase == "standard") names(labels)[1:6] else "mrdag"
    additional_tables_assert(setequal(names(g$estimates), expected_keys) &&
      setequal(names(g$biases), expected_keys) &&
      setequal(g$result_specification$result_key, expected_keys) &&
      !anyDuplicated(g$result_specification$result_key), "Unexpected point result keys.")
    infer_keys <- expected_keys[!grepl("sparse_mr_rr$", expected_keys)]
    for (component in c("standard_errors", "ci_lower", "ci_upper", "coverage")) {
      additional_tables_assert(setequal(names(b[[component]]), infer_keys),
        paste("Unexpected bootstrap methods in", component))
    }
    additional_tables_assert(is.data.frame(b$errors) && nrow(b$errors) == 0L &&
      is.matrix(b$successful_draws) && ncol(b$successful_draws) == n &&
      setequal(rownames(b$successful_draws), infer_keys) &&
      all(b$successful_draws == B) && b$bootstrap_size == B,
      "Bootstrap draws are missing or failed; tables were not generated.")

    truth <- as.vector(g$truth$C)
    for (key in expected_keys) {
      spec <- g$result_specification[g$result_specification$result_key == key, , drop = FALSE]
      method <- sub("^r[123]__", "", key)
      additional_tables_assert(nrow(spec) == 1L && identical(spec$method, method),
        "Method label does not match its result key.")
      estimates <- g$estimates[[key]]
      biases <- g$biases[[key]]
      additional_tables_matrix(estimates, n, paste(group_key, key, "estimates"))
      additional_tables_matrix(biases, n, paste(group_key, key, "biases"))
      difference <- additional_tables_close(biases, sweep(estimates, 1L, truth, "-"),
        paste(group_key, key, "stored bias versus full truth"), tolerance = 1e-8)
      largest_bias_difference <- max(largest_bias_difference, difference)
      signed_bias <- rowMeans(biases)
      abs_bias <- abs(signed_bias)
      empirical_sd <- apply(estimates, 1L, stats::sd)
      se <- cp <- rep(NA_real_, 27L)
      if (method != "sparse_mr_rr") {
        for (component in c("standard_errors", "ci_lower", "ci_upper", "coverage"))
          additional_tables_matrix(b[[component]][[key]], n, paste(key, component))
        additional_tables_assert(all(b$standard_errors[[key]] >= 0) &&
          all(b$ci_lower[[key]] <= b$ci_upper[[key]]), "Invalid SE or confidence interval.")
        expected_coverage <- (b$ci_lower[[key]] <= truth) & (b$ci_upper[[key]] >= truth)
        difference <- additional_tables_close(as.numeric(b$coverage[[key]]),
          as.numeric(expected_coverage), paste(group_key, key, "coverage from intervals"), 0)
        largest_coverage_difference <- max(largest_coverage_difference, difference)
        se <- rowMeans(b$standard_errors[[key]])
        cp <- 100 * rowMeans(b$coverage[[key]])
      }
      index <- index + 1L
      identity <- data.frame(analysis = g$analysis, setting_index = s, scenario = g$scenario,
        measurement_error_weight = unique(cfg$measurement_error_weight),
        genetic_effect_weight = unique(cfg$genetic_effect_weight),
        true_rank = unique(cfg$truth_rank), third_singular_value = expected_delta,
        phase = g$phase, result_key = key, method = method,
        estimator = unname(labels[[method]]), working_rank = spec$working_rank,
        monte_carlo_replicates = n, bootstrap_size = if (method == "sparse_mr_rr") NA_integer_ else B,
        stringsAsFactors = FALSE)
      values <- list(Bias = abs_bias, SD = empirical_sd, SE = se, CP = cp)
      metrics <- unlist(lapply(values, function(x) {
        if (all(is.na(x))) c(med = NA_real_, q1 = NA_real_, q3 = NA_real_)
        else additional_tables_quantiles(x)
      }))
      names(metrics) <- gsub(".", "_", names(metrics), fixed = TRUE)
      rows[[index]] <- cbind(identity, as.data.frame(as.list(metrics)))
      entries[[index]] <- cbind(identity[rep(1L, 27L), , drop = FALSE],
        data.frame(entry_index = 1:27, outcome_index = rep(1:3, 9),
          exposure_index = rep(1:9, each = 3), truth = truth,
          mean_estimate = rowMeans(estimates), signed_bias = signed_bias,
          absolute_bias = abs_bias, SD = empirical_sd, SE = se, CP_percent = cp))
    }
  }

  table <- do.call(rbind, rows)
  entrywise <- do.call(rbind, entries)
  rownames(table) <- rownames(entrywise) <- NULL
  ordering <- order(match(table$analysis, c("rank_misspecification", "approximate_low_rank")),
    table$setting_index, table$working_rank, match(table$method, names(labels)))
  table <- table[ordering, , drop = FALSE]
  rownames(table) <- NULL
  additional_tables_assert(nrow(table) == 52L && nrow(entrywise) == 1404L &&
    sum(table$analysis == "rank_misspecification") == 24L &&
    sum(table$analysis == "approximate_low_rank") == 28L,
    "Final table row counts are incorrect.")
  for (s in 1:4) {
    std <- paste("approximate_low_rank", "standard", s, sep = "__")
    dag <- paste("approximate_low_rank", "mrdag", s, sep = "__")
    additional_tables_assert(identical(point$groups[[std]]$data_seeds,
      point$groups[[dag]]$data_seeds) &&
      identical(bootstrap$groups[[std]]$resample_seeds, bootstrap$groups[[dag]]$resample_seeds),
      "Approximate standard/MrDAG streams differ.")
    additional_tables_close(point$groups[[std]]$truth$C, point$groups[[dag]]$truth$C,
      "Approximate standard/MrDAG truth")
  }
  if (isTRUE(point$metadata$numerical_recovery_applied)) {
    audit <- point$numerical_recovery
    additional_tables_assert(is.list(audit) &&
      isTRUE(audit$previously_successful_entries_identical) &&
      audit$changed_fit_count == point$metadata$recovered_fit_count &&
      length(audit$result_keys) == audit$changed_fit_count,
      "The point result is missing its numerical recovery record.")
  }
  list(table = table, entrywise = entrywise, delta = delta, n = n, B = B,
    max_bias_difference = largest_bias_difference,
    max_coverage_difference = largest_coverage_difference)
}

additional_tables_check_bootstrap_summary <- function(table, recorded, label) {
  ids <- function(x) paste(x$analysis, x$setting_index, x$result_key, sep = "__")
  rows <- table[table$method != "sparse_mr_rr", , drop = FALSE]
  additional_tables_assert(nrow(recorded) == nrow(rows) && !anyDuplicated(ids(recorded)) &&
    setequal(ids(rows), ids(recorded)), paste(label, "has missing or duplicate rows."))
  recorded <- recorded[match(ids(rows), ids(recorded)), , drop = FALSE]
  additional_tables_assert(all(recorded$failed_bootstrap_fits == 0), "Recorded bootstrap failures exist.")
  for (metric in c("SE", "CP")) {
    digits <- if (metric == "CP") 1L else 3L
    for (q in c("med", "q1", "q3")) {
      column <- paste(metric, q, sep = "_")
      additional_tables_close(round(rows[[column]], digits), recorded[[column]],
        paste(label, column), 1e-12)
    }
  }
}

additional_tables_latex <- function(table, analysis) {
  x <- table[table$analysis == analysis, , drop = FALSE]
  rank_table <- analysis == "rank_misspecification"
  columns <- if (rank_table) 6L else 5L
  caption <- if (rank_table) {
    "Sensitivity to working-rank misspecification when the true causal-effect matrix has rank two."
  } else {
    sprintf(paste0("Simulation results under approximate low rank, with true singular values ",
      "$(1,1,%s)$ and working rank two for the reduced-rank estimators."),
      format(unique(x$third_singular_value), trim = TRUE))
  }
  notes <- paste0("\\textit{Note:} Each cell reports the median (first line) and the 25th and 75th ",
    "percentiles (second line) across the 27 entries of the causal-effect matrix. For each entry, ",
    "Bias is the absolute mean estimation error over 1,000 Monte Carlo replicates; SD is the ",
    "empirical standard deviation; SE is the mean bootstrap standard error across replicates; ",
    "and CP is the coverage probability of the nominal 95\\% percentile bootstrap interval, ",
    "reported as a percentage. Each bootstrap uses 300 SNP resamples. ",
    if (rank_table) "The true singular values are $(1,1,0)$. " else
      "Bias and coverage use the full rank-three truth $C_\\delta$, including its third component. ",
    "Bootstrap inference is not computed for Sparse MR-rr and is denoted by \\textemdash.")
  cell <- function(row, metric) {
    v <- as.numeric(row[paste0(metric, c("_med", "_q1", "_q3"))])
    if (all(is.na(v))) return("\\textemdash")
    digits <- if (metric == "CP") 1L else 3L
    strings <- sprintf(paste0("%.", digits, "f"), round(v, digits))
    if (metric == "CP") strings <- paste0(strings, "\\%")
    paste0("\\shortstack{", strings[1L], "\\\\(", strings[2L], ", ", strings[3L], ")}")
  }
  lines <- c("% Generated by 30_make_additional_simulation_tables.R; do not edit numbers by hand.",
    "% Requires \\usepackage{booktabs}. Input this file inside the document body.",
    "\\begin{table}[p]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:additional-", gsub("_", "-", analysis), "}"),
    "\\begingroup", "\\footnotesize", "\\renewcommand{\\arraystretch}{1.08}",
    "\\setlength{\\tabcolsep}{3pt}",
    paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}",
      if (rank_table) "clcccc" else "lcccc", "@{}}"), "\\toprule",
    paste0(if (rank_table) "Working rank & " else "", "Estimator & Bias & SD & SE & CP \\\\"),
    "\\midrule")
  for (s in 1:4) {
    part <- x[x$setting_index == s, , drop = FALSE]
    me <- format(unique(part$measurement_error_weight), nsmall = 1, trim = TRUE)
    ge <- format(unique(part$genetic_effect_weight), nsmall = 2, trim = TRUE)
    if (s > 1L) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf(paste0("\\multicolumn{%d}{@{}l}{\\textbf{Setting %d:} ",
      "$\\Sigma_X$ multiplier $=%s$, $\\Sigma_{\\gamma\\gamma}$ multiplier $=%s$} \\\\"),
      columns, s, me, ge))
    for (i in seq_len(nrow(part))) {
      row <- part[i, , drop = FALSE]
      contents <- c(if (rank_table) as.character(row$working_rank), row$estimator,
        vapply(c("Bias", "SD", "SE", "CP"), function(m) cell(row, m), character(1)))
      lines <- c(lines, paste0(paste(contents, collapse = " & "), " \\\\"))
    }
  }
  c(lines, "\\bottomrule", "\\end{tabular*}", "\\par\\smallskip",
    "\\begin{minipage}{\\textwidth}", "\\footnotesize", notes,
    "\\end{minipage}", "\\endgroup", "\\end{table}")
}

additional_tables_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- additional_tables_options(args)
  if (is.null(opts)) return(invisible(NULL))
  input <- normalizePath(opts[["input-dir"]], winslash = "/", mustWork = TRUE)
  output <- normalizePath(opts[["output-dir"]], winslash = "/", mustWork = FALSE)
  additional_tables_assert(!grepl("(^|/)freeze(/|$)", tolower(output)),
    "Writing into a frozen snapshot is not allowed.")
  files <- c(point = file.path(input, "point_merged_recovered/additional_point_results.rds"),
    bootstrap = file.path(input, "bootstrap_merged/additional_bootstrap_results.rds"),
    bootstrap_csv = file.path(input, "bootstrap_merged/additional_bootstrap_summary.csv"))
  additional_tables_assert(all(file.exists(files)), "One or more merged input files are missing.")
  input_hashes <- tools::md5sum(files)
  point <- readRDS(files[["point"]])
  bootstrap <- readRDS(files[["bootstrap"]])
  result <- additional_tables_summarize(point, bootstrap)
  additional_tables_check_bootstrap_summary(result$table, bootstrap$table_summary, "Embedded bootstrap summary")
  additional_tables_check_bootstrap_summary(result$table,
    read.csv(files[["bootstrap_csv"]], stringsAsFactors = FALSE), "Saved bootstrap summary CSV")
  names_out <- c("table_rank_misspecification.tex", "table_approximate_low_rank.tex",
    "additional_tables_summary.csv", "additional_tables_entrywise.csv",
    "additional_tables_preview.tex", "table_validation.txt", "table_provenance.rds")
  paths_out <- file.path(output, names_out)
  additional_tables_assert(opts$overwrite || !any(file.exists(paths_out)),
    "Generated table files already exist. Use --overwrite=true to regenerate them.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  for (analysis in c("rank_misspecification", "approximate_low_rank")) {
    writeLines(additional_tables_latex(result$table, analysis),
      file.path(output, paste0("table_", analysis, ".tex")), useBytes = TRUE)
  }
  write.csv(result$table, file.path(output, "additional_tables_summary.csv"), row.names = FALSE, na = "")
  write.csv(result$entrywise, file.path(output, "additional_tables_entrywise.csv"), row.names = FALSE, na = "")
  preview <- c("\\documentclass[10pt,letterpaper]{article}",
    "\\usepackage[margin=0.7in]{geometry}", "\\usepackage{booktabs}",
    "\\usepackage[T1]{fontenc}", "\\begin{document}",
    "\\input{table_rank_misspecification.tex}", "\\clearpage",
    "\\input{table_approximate_low_rank.tex}", "\\end{document}")
  writeLines(preview, file.path(output, "additional_tables_preview.tex"), useBytes = TRUE)
  recovery_count <- if (isTRUE(point$metadata$numerical_recovery_applied))
    point$metadata$recovered_fit_count else 0L
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else NA_character_
  if (is.na(script_file)) {
    source_file <- getSrcFilename(additional_tables_main, full.names = TRUE)
    if (length(source_file) == 1L && nzchar(source_file)) script_file <- source_file
  }
  provenance <- list(schema_version = "1.0.0", generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    inputs = data.frame(component = names(files), filename = basename(files), md5 = unname(input_hashes)),
    point_metadata = point$metadata, bootstrap_metadata = bootstrap$metadata,
    numerical_recovery = point$numerical_recovery,
    table_definitions = c(Bias = "Across-entry quantiles of abs(rowMeans(estimation errors))",
      SD = "Across-entry quantiles of sample SD across Monte Carlo replicates",
      SE = "Across-entry quantiles of mean bootstrap SE across Monte Carlo replicates",
      CP = "Across-entry quantiles of 100 times empirical coverage; nominal percentile CI level 95%"),
    quantile_type = 7L, script_file = script_file,
    script_md5 = if (!is.na(script_file) && file.exists(script_file)) unname(tools::md5sum(script_file)) else NA_character_,
    session_info = capture.output(sessionInfo()))
  saveRDS(provenance, file.path(output, "table_provenance.rds"))
  additional_tables_assert(identical(tools::md5sum(files), input_hashes), "An input file changed during table generation.")
  checks <- c("ADDITIONAL MANUSCRIPT TABLES: PASS", "Monte Carlo replicates per setting: 1000",
    "Bootstrap draws per replicate: 300", paste("Approximate delta:", result$delta),
    "Point/bootstrap truth and data-seed alignment: PASS", "Complete point estimates and bootstrap draws: PASS",
    paste("Maximum stored-bias versus full-truth difference:", result$max_bias_difference),
    paste("Maximum coverage versus interval-derived coverage difference:", result$max_coverage_difference),
    "SE/CP versus existing bootstrap RDS summary and CSV (all 36 rows): PASS",
    "Rank-misspecification table: 24 rows", "Approximate-low-rank table: 28 rows",
    "Entrywise output: 1404 rows", "Sparse SE/CP: not computed, left missing",
    paste("Numerically recovered point fits retained:", recovery_count),
    "Input files unchanged: PASS", "No simulations or bootstrap fits were rerun.")
  writeLines(checks, file.path(output, "table_validation.txt"), useBytes = TRUE)
  cat(paste(checks, collapse = "\n"), "\nOutput:", normalizePath(output, winslash = "/"), "\n")
  invisible(result)
}

if (sys.nframe() == 0L) additional_tables_main()
