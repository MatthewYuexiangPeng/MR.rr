#!/usr/bin/env Rscript
# Export the unified spectral simulation: base R only; no fits or RNG calls.
# The completed merged RDS is the sole numerical input. Run with --help.

spout_source <- local({
  frames <- sys.frames()
  files <- unlist(lapply(frames, function(f) f$ofile), use.names = FALSE)
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  f <- if (length(files)) tail(files, 1L) else if (length(arg)) sub("^--file=", "", arg[1L]) else ""
  if (nzchar(f) && file.exists(f)) normalizePath(f, winslash = "/") else NA_character_
})
spout_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
spout_equal <- function(x, y, message, tolerance = 1e-10) {
  spout_assert(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)), message)
}
spout_labels <- c(ivw = "IVW", srivw = "SRIVW", naive_mr_rr = "Naive MR-rr",
  mr_rr = "MR-rr", regularized_mr_rr = "Reg. MR-rr", sparse_mr_rr = "Sparse MR-rr", mrdag = "MrDAG")
spout_colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
spout_key <- function(design, setting, method, rank = 2L) paste(design, paste0("setting-", setting),
  method, paste0("rank-", if (method %in% c("ivw", "srivw", "mrdag")) 0L else rank), sep = "__")
spout_options <- function(args) {
  opts <- list(input = "paper/output/spectral_rebuild/cluster_simulations_v1_blas_recovery/run/merged/spectral_simulation_results.rds",
    output = "paper/output/spectral_rebuild/manuscript_simulations_v1", `support-threshold` = "0.001")
  if (identical(args, "--help")) {
    cat("Usage: Rscript --vanilla paper/scripts/36_make_spectral_simulation_outputs.R\n",
      "  [--input=MERGED_RDS] [--output=NEW_OUTPUT_DIRECTORY] [--support-threshold=0.001]\n",
      "The output directory must not already exist. Requires base R only.\n")
    return(NULL)
  }
  seen <- character()
  for (arg in args) {
    spout_assert(grepl("^--[^=]+=.+$", arg), paste("Use --name=value:", arg))
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    spout_assert(key %in% names(opts) && !key %in% seen, paste("Unknown or duplicate option:", key))
    opts[[key]] <- sub("^--[^=]+=", "", arg); seen <- c(seen, key)
  }
  opts$`support-threshold` <- as.numeric(opts$`support-threshold`)
  spout_assert(length(opts$`support-threshold`) == 1L && is.finite(opts$`support-threshold`) &&
    opts$`support-threshold` > 0, "Support threshold must be positive and finite.")
  opts
}
spout_matrix <- function(x, nr, nc, label) {
  spout_assert(is.matrix(x) && identical(dim(x), as.integer(c(nr, nc))) && all(is.finite(x)),
    paste("Invalid matrix:", label))
}

spout_validate <- function(x) {
  spout_assert(identical(x$schema, "spectral-merged-1") && identical(x$seal$validation_only, FALSE),
    "A completed production spectral-merged-1 object is required.")
  spout_assert(identical(x$run_id, x$seal$run_id), "Run identities disagree.")
  spout_assert(nrow(x$config) == 12L && all(x$config$replicates == 1000L) &&
    all(x$config$bootstrap_size == 300L), "This paper export requires 12 cells, N=1000 and B=300.")
  expected <- character()
  for (design in c("generic", "sparse_loading", "approximate_low_rank")) {
    spout_matrix(x$truths[[design]]$C, 3L, 9L, paste(design, "truth"))
    for (s in 1:4) for (method in names(spout_labels)) {
      ranks <- if (design == "generic" && method %in% c("regularized_mr_rr", "sparse_mr_rr")) 1:3 else 2L
      for (rank in ranks) expected <- c(expected, spout_key(design, s, method, rank))
    }
  }
  spout_assert(nrow(x$catalog) == 100L && !anyDuplicated(x$catalog$result_key) &&
    setequal(expected, x$catalog$result_key) && setequal(expected, names(x$results)), "Result catalog is incomplete.")
  spout_assert(!anyDuplicated(paste(x$config$design, x$config$setting)), "Duplicate design cells.")
  cells <- unlist(lapply(c("generic", "sparse_loading", "approximate_low_rank"), function(d) paste(d, 1:4, sep = "/")))
  for (field in c("data_md5", "resample_md5")) {
    spout_assert(setequal(names(x[[field]]), cells) && all(vapply(x[[field]], function(h)
      length(h) == 1000L && all(grepl("^[0-9a-f]{32}$", h)), logical(1))), paste("Incomplete", field))
  }
  spout_assert(length(x$prediction_exposure) == 9L && all(is.finite(x$prediction_exposure)), "Missing prediction exposure profile.")
  spout_equal(svd(x$truths$generic$C, nu = 0, nv = 0)$d, c(1, 1, 0), "Generic truth singular values differ.")
  spout_equal(svd(x$truths$approximate_low_rank$C, nu = 0, nv = 0)$d, c(1, 1, .1), "Approximate truth singular values differ.")
  entries <- summaries <- list()
  for (i in seq_len(nrow(x$catalog))) {
    row <- x$catalog[i, ]; key <- row$result_key; r <- x$results[[key]]
    truth <- as.vector(x$truths[[row$design]]$C)
    spout_assert(identical(key, spout_key(row$design, row$setting, row$method, row$working_rank)) &&
      identical(row$bootstrap, row$method != "sparse_mr_rr"), "Invalid method/rank/bootstrap mapping.")
    spout_matrix(r$point, 27L, 1000L, key)
    spout_assert(length(r$point_diagnostics) == 1000L, paste("Missing point diagnostics:", key))
    e <- data.frame(result_key = key, outcome = rep(1:3, 9), exposure = rep(1:9, each = 3),
      truth = truth, signed_bias = rowMeans(r$point) - truth, bias = abs(rowMeans(r$point) - truth),
      sd = apply(r$point, 1, stats::sd), rmse = sqrt(rowMeans((r$point - truth)^2)), se = NA_real_, cp_percent = NA_real_)
    if (row$bootstrap) {
      for (field in c("se", "lower", "upper", "coverage")) spout_matrix(r[[field]], 27L, 1000L, paste(key, field))
      spout_assert(all(r$se >= 0) && all(r$lower <= r$upper), paste("Invalid bootstrap intervals:", key))
      spout_assert(identical(as.numeric(r$coverage), as.numeric(r$lower <= truth & r$upper >= truth)),
        paste("Coverage disagrees with percentile intervals:", key))
      spout_assert(length(r$bootstrap_issues) == 1000L, paste("Missing bootstrap diagnostics:", key))
      e$se <- rowMeans(r$se); e$cp_percent <- 100 * rowMeans(r$coverage)
    } else spout_assert(all(vapply(c("se", "lower", "upper", "coverage"), function(f) is.null(r[[f]]), logical(1))),
      "Sparse simulation bootstrap must not be fabricated.")
    sr <- row
    for (metric in c("bias", "sd", "se", "cp_percent", "rmse")) {
      vals <- if (all(is.na(e[[metric]]))) rep(NA_real_, 3) else
        stats::quantile(e[[metric]], c(.5, .25, .75), names = FALSE, type = 7L)
      for (j in 1:3) sr[[paste0(metric, c("_median", "_q25", "_q75")[j])]] <- vals[j]
    }
    sr$nonconverged_point_fits <- sum(vapply(r$point_diagnostics, function(d) identical(d$details$converged, FALSE), logical(1)))
    sr$projected_point_fits <- sum(vapply(r$point_diagnostics, function(d) isTRUE(d$details$numerical_diagnostics$corrected_covariance_projected), logical(1)))
    sr$point_warning_count <- sum(vapply(r$point_diagnostics, function(d) length(d$warnings), integer(1)))
    sr$bootstrap_draws_with_issues <- if (row$bootstrap) sum(vapply(r$bootstrap_issues, length, integer(1))) else 0L
    entries[[i]] <- e; summaries[[i]] <- sr
  }
  summary <- do.call(rbind, summaries); entrywise <- do.call(rbind, entries)
  spout_equal(summary, x$summary[match(summary$result_key, x$summary$result_key), names(summary)],
    "Recomputed summaries disagree with the original strict merger.")
  expected_map <- x$catalog[x$catalog$working_rank %in% c(0, 2), ]
  expected_map$table_id <- c(generic = "main_generic", sparse_loading = "main_sparse_loading",
    approximate_low_rank = "approximate_low_rank")[expected_map$design]
  sensitivity <- x$catalog[x$catalog$design == "generic" & x$catalog$method %in% c("regularized_mr_rr", "sparse_mr_rr"), ]
  sensitivity$table_id <- "rank_misspecification"
  expected_map <- rbind(expected_map, sensitivity)
  expected_map$point_result_key <- expected_map$result_key
  expected_map$bootstrap_result_key <- ifelse(expected_map$bootstrap, expected_map$result_key, NA_character_)
  map_columns <- names(expected_map)
  spout_assert(nrow(x$table_map) == 108L && !anyDuplicated(paste(x$table_map$table_id, x$table_map$result_key)), "Invalid table map.")
  spout_equal(expected_map, x$table_map[, map_columns], "The table map differs from the canonical shared-result mapping.", 0)
  tr <- cbind(expected_map[, c("table_id", "point_result_key", "bootstrap_result_key")],
    summary[match(expected_map$result_key, summary$result_key), ])
  spout_equal(tr, x$table_rows[, names(tr)], "Stored table rows disagree with the recomputed rows.")
  shared <- intersect(tr$result_key[tr$table_id == "main_generic"], tr$result_key[tr$table_id == "rank_misspecification"])
  spout_assert(length(shared) == 8L, "Expected eight shared rank-two rows.")
  for (key in shared) {
    a <- tr[tr$result_key == key, names(summary)]
    spout_assert(identical(unname(as.list(a[1, ])), unname(as.list(a[2, ]))), "Shared table rows are not exactly equal.")
  }
  if (!is.null(x$recovery_provenance)) {
    p <- x$recovery_provenance
    spout_assert(identical(p$schema, "spectral-blas-recovery-1") && identical(p$original_run_id, x$run_id), "Invalid recovery identity.")
    spout_assert(length(intersect(p$recomputed_task_ids, p$retained_task_ids)) == 0L &&
      setequal(c(p$recomputed_task_ids, p$retained_task_ids), 1:1800), "Recovery task coverage is incomplete.")
  }
  list(summary = summary, entrywise = entrywise, table_rows = tr, shared = shared)
}

spout_table <- function(rows, table_id, x, out) {
  meta <- list(
    main_generic = c("table_main_generic.tex", "tab:main_regular_C", "Finite-sample performance under the generic low-rank design across four levels of instrument strength."),
    main_sparse_loading = c("table_main_sparse_loading.tex", "tab:main_sparse_loading_C", "Finite-sample performance under the sparse-loading design across four levels of instrument strength."),
    rank_misspecification = c("table_rank_misspecification.tex", "tab:rank_misspecification", "Sensitivity to working-rank misspecification when the true causal-effect matrix has rank two."),
    approximate_low_rank = c("table_approximate_low_rank.tex", "tab:approximate_low_rank", "Performance under an approximately rank-two causal-effect matrix with singular values $(1,1,0.1)$."))[[table_id]]
  has_rank <- table_id == "rank_misspecification"; nc <- if (has_rank) 6L else 5L
  cell <- function(r, metric) {
    if (is.na(r[[paste0(metric, "_median")]])) return("---")
    v <- as.numeric(r[paste0(metric, c("_median", "_q25", "_q75"))])
    fmt <- if (metric == "cp_percent") "%.1f" else "%.3f"
    sprintf("\\shortstack{%s\\\\(%s, %s)}", sprintf(fmt, v[1]), sprintf(fmt, v[2]), sprintf(fmt, v[3]))
  }
  lines <- c("\\begin{table}[p]", "\\centering", paste0("\\caption{", meta[3], "}\\label{", meta[2], "}"),
    "\\begingroup\\fontsize{8.7}{9.6}\\selectfont", "\\setlength{\\tabcolsep}{7pt}\\renewcommand{\\arraystretch}{0.94}",
    paste0("\\begin{tabular}{", if (has_rank) "clrrrr" else "lrrrr", "}\\toprule"),
    paste0(if (has_rank) "Working rank & " else "", "Estimator & Bias & SD & SE & CP (\\%) \\\\ \\midrule"))
  rr <- rows[rows$table_id == table_id, ]
  for (s in 1:4) {
    group <- rr[rr$setting == s, ]; design <- unique(group$design)
    cfg <- x$config[x$config$design == design & x$config$setting == s, ]
    siv <- x$parameters[[paste(design, s, sep = "/")]]$iv_strength
    lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{Setting %d:} $S_{\\mathrm{IV}}=%.2f$; $\\Sigma_X\\times%.1f$, $\\Sigma_{\\gamma\\gamma}\\times%.2f$} \\\\[2pt]", nc, s, siv, cfg$measurement_error_weight, cfg$genetic_effect_weight))
    group <- group[if (has_rank) order(group$working_rank, match(group$method, names(spout_labels))) else
      order(match(group$method, names(spout_labels))), ]
    for (i in seq_len(nrow(group))) {
      r <- group[i, ]; values <- c(if (has_rank) r$working_rank, spout_labels[[r$method]],
        vapply(c("bias", "sd", "se", "cp_percent"), function(m) cell(r, m), character(1)))
      lines <- c(lines, paste0(paste(values, collapse = " & "), " \\\\[1.5pt]"))
    }
    lines <- c(lines, if (s < 4) "\\midrule" else "\\bottomrule")
  }
  note <- paste0("Each cell gives the median and interquartile range across the 27 entries of $C$. ",
    "For each entry, Bias is the absolute mean estimation error over 1,000 Monte Carlo replicates; SD is the empirical standard deviation; ",
    "SE is the mean bootstrap standard error; and CP is coverage of a nominal 95\\% percentile interval, expressed as a percentage. ",
    "Each bootstrap uses 300 SNP resamples. Sparse MR-rr has no bootstrap inference (---). ",
    if (has_rank) "The working-rank-two rows use exactly the same estimates and bootstrap summaries as Table~\\ref{tab:main_regular_C}." else "")
  lines <- c(lines, "\\end{tabular}\\endgroup", "\\par\\vspace{4pt}\\begin{minipage}{0.98\\textwidth}\\footnotesize",
    paste0("\\textit{Note:} ", note), "\\end{minipage}", "\\end{table}")
  writeLines(lines, file.path(out, meta[1]), useBytes = TRUE)
  meta[1]
}

spout_predict <- function(mat, exposure) {
  # Explicit exposure-order summation: avoids dispatching a new BLAS product.
  ans <- matrix(0, 3L, ncol(mat))
  for (j in 1:9) ans <- ans + mat[(3L * (j - 1L) + 1L):(3L * j), , drop = FALSE] * exposure[j]
  ans
}
spout_align <- function(A, B, truth) {
  spout_matrix(A, 3L, 2L, "sparse A"); spout_matrix(B, 2L, 9L, "sparse B")
  candidates <- expand.grid(s1 = c(1, -1), s2 = c(1, -1), swap = c(FALSE, TRUE))
  best <- Inf; result <- NULL
  for (i in seq_len(nrow(candidates))) {
    signs <- as.numeric(candidates[i, 1:2]); perm <- if (candidates$swap[i]) 2:1 else 1:2
    bb <- B[perm, , drop = FALSE] * signs
    loss <- sum((bb - truth)^2)
    if (loss < best) {
      best <- loss
      result <- list(A = sweep(A[, perm, drop = FALSE], 2, signs, "*"), B = bb,
        swap = candidates$swap[i], s1 = signs[1], s2 = signs[2], loss = loss)
    }
  }
  result
}

spout_support <- function(x, threshold) {
  truth <- x$truths$sparse_loading$B; spout_matrix(truth, 2L, 9L, "sparse truth B")
  truth_nonzero <- as.vector(truth != 0); fits <- summaries <- entries <- aligned <- list()
  max_difference <- 0
  for (s in 1:4) {
    key <- spout_key("sparse_loading", s, "sparse_mr_rr"); r <- x$results[[key]]
    bm <- matrix(NA_real_, 18L, 1000L); transformations <- vector("list", 1000L)
    for (i in 1:1000) {
      d <- r$point_diagnostics[[i]]$details; a <- spout_align(d$A, d$B, truth)
      old_c <- d$A %*% d$B; new_c <- a$A %*% a$B
      spout_equal(as.vector(old_c), r$point[, i], "Saved sparse factors disagree with point C.")
      spout_equal(old_c, new_c, "Permutation/sign alignment changed the fitted C.")
      max_difference <- max(max_difference, abs(old_c - new_c))
      bm[, i] <- as.vector(a$B)
      transformations[[i]] <- data.frame(setting = s, replicate = i, swap = a$swap,
        sign_1 = a$s1, sign_2 = a$s2, squared_distance = a$loss)
    }
    pred <- abs(bm) > threshold
    tp <- colSums(pred & truth_nonzero); fp <- colSums(pred & !truth_nonzero)
    tn <- colSums(!pred & !truth_nonzero); fn <- colSums(!pred & truth_nonzero)
    rep <- do.call(rbind, transformations)
    rep$TP <- tp; rep$FP <- fp; rep$TN <- tn; rep$FN <- fn
    rep$exact_support <- fp + fn == 0
    rep$sensitivity <- tp / (tp + fn); rep$specificity <- tn / (tn + fp)
    rep$precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)
    rep$FDR <- ifelse(tp + fp > 0, fp / (tp + fp), NA_real_)
    summaries[[s]] <- data.frame(setting = s, support_threshold = threshold, replicates = 1000,
      exact_count = sum(rep$exact_support), exact_percent = 100 * mean(rep$exact_support),
      sensitivity = sum(tp) / sum(tp + fn), specificity = sum(tn) / sum(tn + fp),
      precision = sum(tp) / sum(tp + fp), FDR = sum(fp) / sum(tp + fp))
    entries[[s]] <- data.frame(setting = s, pathway = rep(1:2, 9), exposure = rep(1:9, each = 2),
      truth = as.vector(truth), mean = rowMeans(bm), sd = apply(bm, 1, stats::sd),
      median = apply(bm, 1, stats::median), q25 = apply(bm, 1, stats::quantile, .25),
      q75 = apply(bm, 1, stats::quantile, .75), selected_proportion = rowMeans(pred))
    fits[[s]] <- rep; aligned[[paste0("setting-", s)]] <- bm
  }
  list(replicates = do.call(rbind, fits), summary = do.call(rbind, summaries),
    entrywise = do.call(rbind, entries), aligned = aligned, max_C_difference = max_difference)
}

spout_limits <- function(values, truth) {
  lim <- range(stats::quantile(values, c(.01, .99), names = FALSE), truth)
  if (diff(lim) == 0) lim <- lim + c(-1, 1)
  lim + c(-1, 1) * diff(lim) * .06
}
spout_box <- function(values, lim, truth, names, colors, title, yaxis = TRUE, xaxis = FALSE) {
  graphics::boxplot(values, names = rep("", length(values)), ylim = lim, outline = TRUE,
    col = colors, border = "#444444", axes = FALSE, main = title, cex.main = .83,
    pars = list(outpch = 16, outcex = .13, boxwex = .58, medlwd = 1.1), yaxs = "i")
  graphics::abline(h = truth, col = if (truth == 0 && length(values) == 1) "#888888" else "#B22222", lty = 2, lwd = 1)
  if (yaxis) graphics::axis(2, las = 1, cex.axis = .67, tck = -.025)
  if (xaxis) graphics::axis(1, at = seq_along(values), labels = names, cex.axis = .67, tck = -.015)
  graphics::box(col = "#CCCCCC")
}
spout_box_record <- function(values, lim, figure, panel, keys) do.call(rbind, lapply(seq_along(values), function(i) {
  a <- values[[i]]; b <- grDevices::boxplot.stats(a)
  data.frame(figure = figure, panel = panel, result_key = keys[i], n = length(a),
    minimum = min(a), lower_whisker = b$stats[1], lower_hinge = b$stats[2], median = b$stats[3],
    upper_hinge = b$stats[4], upper_whisker = b$stats[5], maximum = max(a),
    axis_lower = lim[1], axis_upper = lim[2], outside_axis = sum(a < lim[1] | a > lim[2]))
}))
spout_legend <- function() {
  graphics::par(oma = rep(0, 4), fig = c(0, 1, 0, .085), new = TRUE, mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend("center", legend = paste(1:7, spout_labels), fill = spout_colors,
    ncol = 4, bty = "n", cex = .88, border = NA)
}

spout_figures <- function(x, support, out) {
  stats <- pred_entries <- captions <- list(); index <- 0L
  for (design in c("generic", "sparse_loading")) for (s in c(4L, 1L)) {
    keys <- vapply(names(spout_labels), function(m) spout_key(design, s, m), character(1))
    points <- lapply(keys, function(k) x$results[[k]]$point)
    truth <- x$truths[[design]]$C
    title <- paste(if (design == "generic") "Generic low-rank design" else "Sparse-loading design",
      sprintf(" | Setting %d | SIV = %.2f", s, x$parameters[[paste(design, s, sep = "/")]]$iv_strength))
    file <- sprintf("figure_C_%s_setting%d.pdf", design, s)
    grDevices::pdf(file.path(out, file), width = 12.4, height = 6.6, pointsize = 11, useDingbats = FALSE)
    graphics::par(mfrow = c(3, 9), mar = c(1.4, 2.4, 1.65, .2), oma = c(3.4, 1, 3.2, .3), mgp = c(1.2, .3, 0), cex = 1)
    for (y in 1:3) {
      rows <- seq(y, 27, by = 3)
      lim <- spout_limits(unlist(lapply(points, function(p) p[rows, ])), truth[y, ])
      for (j in 1:9) {
        vals <- lapply(points, function(p) p[y + 3 * (j - 1), ])
        spout_box(vals, lim, truth[y, j], 1:7, spout_colors, paste0("Y", y, " / X", j), xaxis = y == 3)
        index <- index + 1L; stats[[index]] <- spout_box_record(vals, lim, file, paste(y, j, sep = "/"), keys)
      }
    }
    graphics::mtext(title, outer = TRUE, side = 3, line = 1.5, cex = 1.08)
    graphics::mtext("C estimates; red dashed lines = truth", outer = TRUE, side = 3, line = .1, cex = .86)
    spout_legend(); grDevices::dev.off()
    captions[[file]] <- paste(title, ". Entrywise estimates of C. Red dashed lines indicate truth. Display limits use pooled 1st/99th percentiles within each outcome, extended to include truth; all 1,000 replicates enter boxplot and table statistics.")
    pred <- lapply(points, spout_predict, exposure = x$prediction_exposure)
    target <- as.vector(spout_predict(matrix(as.vector(truth), 27, 1), x$prediction_exposure))
    for (m in seq_along(pred)) for (y in 1:3) {
      v <- pred[[m]][y, ]; err <- v - target[y]
      pred_entries[[length(pred_entries) + 1L]] <- data.frame(design = design, setting = s,
        result_key = keys[m], method = names(spout_labels)[m], outcome = y, truth = target[y],
        mean = mean(v), signed_bias = mean(err), absolute_bias = abs(mean(err)), sd = stats::sd(v),
        rmse = sqrt(mean(err^2)), median = stats::median(v), q25 = stats::quantile(v, .25), q75 = stats::quantile(v, .75))
    }
    file <- sprintf("figure_prediction_%s_setting%d.pdf", design, s)
    grDevices::pdf(file.path(out, file), width = 9.6, height = 4.4, pointsize = 11, useDingbats = FALSE)
    graphics::par(mfrow = c(1, 3), mar = c(2.5, 3.2, 2, .7), oma = c(3.4, 0, 3, 0), mgp = c(1.7, .5, 0), cex = 1)
    for (y in 1:3) {
      vals <- lapply(pred, function(p) p[y, ]); lim <- spout_limits(unlist(vals), target[y])
      spout_box(vals, lim, target[y], 1:7, spout_colors, paste("Outcome", y), xaxis = TRUE)
      index <- index + 1L; stats[[index]] <- spout_box_record(vals, lim, file, as.character(y), keys)
    }
    graphics::mtext(title, outer = TRUE, side = 3, line = 1.3, cex = 1.04)
    graphics::mtext("Predicted risk scores; red dashed lines = true Cx", outer = TRUE, side = 3, line = -.1, cex = .85)
    spout_legend(); grDevices::dev.off()
    captions[[file]] <- paste(title, ". Predicted risk scores using the stored fixed exposure vector. Red dashed lines indicate true Cx. Display limits use pooled 1st/99th percentiles within each outcome, extended to include truth; statistics use all 1,000 replicates.")
  }
  file <- "figure_B_sparse_loading_setting4.pdf"; bm <- support$aligned[["setting-4"]]; truth <- x$truths$sparse_loading$B
  grDevices::pdf(file.path(out, file), width = 12.4, height = 4.6, pointsize = 11, useDingbats = FALSE)
  graphics::par(mfrow = c(2, 9), mar = c(1, 2.6, 1.7, .2), oma = c(.5, .5, 3.7, 0), mgp = c(1.3, .3, 0), cex = 1)
  for (k in 1:2) {
    rows <- seq(k, 18, by = 2); lim <- spout_limits(bm[rows, ], truth[k, ])
    for (j in 1:9) {
      vals <- list(bm[k + 2 * (j - 1), ])
      spout_box(vals, lim, truth[k, j], "", "#AA3377", paste0("B", k, " / X", j))
      index <- index + 1L; stats[[index]] <- spout_box_record(vals, lim, file, paste(k, j, sep = "/"),
        spout_key("sparse_loading", 4, "sparse_mr_rr"))
    }
  }
  graphics::mtext("Sparse loading estimates | Setting 4 | 1,000 replicates", outer = TRUE, side = 3, line = 2, cex = 1.12)
  graphics::mtext("Row permutations and signs aligned to truth; A transformed jointly to preserve C", outer = TRUE, side = 3, line = .8, cex = .88)
  graphics::mtext("Dashed lines: red = nonzero truth; gray = zero truth", outer = TRUE, side = 3, line = -.3, cex = .83)
  grDevices::dev.off()
  captions[[file]] <- "Sparse-loading design, Setting 4. Estimated B is aligned to the stored truth over all row permutations and sign changes, minimizing squared Frobenius distance. The same transformation is applied to A. Dashed lines indicate true loadings (red: nonzero; gray: zero). Display limits use pooled 1st/99th percentiles within each pathway, extended to include truth; statistics use all replicates."
  list(boxplots = do.call(rbind, stats), prediction = do.call(rbind, pred_entries), captions = captions)
}

spout_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- spout_options(args); if (is.null(opts)) return(invisible(NULL))
  spout_assert(file.exists(opts$input), paste("Missing merged RDS:", opts$input))
  spout_assert(!file.exists(opts$output) && !dir.exists(opts$output), "Output already exists; choose a new directory.")
  input <- normalizePath(opts$input, winslash = "/"); input_md5 <- unname(tools::md5sum(input))
  x <- readRDS(input); v <- spout_validate(x)
  cat("PASS: production result completeness, intervals and independently recomputed summaries\n")
  cat("PASS: eight main/rank-sensitivity rows share exactly the same result keys and values\n")
  support <- spout_support(x, opts$`support-threshold`)
  cat("PASS: sparse loading alignment preserves every stored C\n")
  parent <- dirname(opts$output); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile("spectral_export_build_", tmpdir = parent); dir.create(stage)
  # Keep an incomplete staging directory on error for diagnosis; never label it PASS.
  for (d in c("tables", "figures", "statistics", "provenance")) dir.create(file.path(stage, d))
  tables <- vapply(c("main_generic", "main_sparse_loading", "rank_misspecification", "approximate_low_rank"),
    function(id) spout_table(v$table_rows, id, x, file.path(stage, "tables")), character(1))
  plots <- spout_figures(x, support, file.path(stage, "figures"))
  write_csv <- function(object, name) utils::write.csv(object, file.path(stage, "statistics", name), row.names = FALSE, na = "NA")
  write_csv(v$summary, "simulation_summary.csv"); write_csv(v$table_rows, "simulation_table_rows.csv")
  write_csv(v$entrywise, "simulation_entrywise.csv"); write_csv(plots$prediction, "prediction_summary.csv")
  write_csv(plots$boxplots, "figure_boxplot_statistics.csv"); write_csv(support$summary, "sparse_support_summary.csv")
  write_csv(support$replicates, "sparse_support_replicates.csv"); write_csv(support$entrywise, "sparse_B_entrywise.csv")
  write_csv(data.frame(exposure = 1:9, value = x$prediction_exposure), "prediction_exposure.csv")
  write_csv(data.frame(result_key = v$shared, main_table = "main_generic", sensitivity_table = "rank_misspecification",
    same_point_object = TRUE, same_bootstrap_object = !grepl("sparse_mr_rr", v$shared), exact_summary_match = TRUE), "shared_rank2_rows.csv")
  saveRDS(list(truth = x$truths$sparse_loading, aligned_B = support$aligned,
    alignment = "Minimum squared Frobenius distance over signed row permutations", support_threshold = opts$`support-threshold`),
    file.path(stage, "statistics", "sparse_B_aligned.rds"), version = 2L)
  table_preview <- c("\\documentclass[10pt,a4paper]{article}", "\\usepackage[margin=15mm]{geometry}",
    "\\usepackage{booktabs,amsmath,lmodern}", "\\usepackage[T1]{fontenc}", "\\begin{document}",
    unlist(lapply(tables, function(f) c(paste0("\\input{", f, "}"), "\\clearpage"))), "\\end{document}")
  writeLines(table_preview, file.path(stage, "tables", "tables_preview.tex"))
  esc <- function(s) gsub("_", "\\_", s, fixed = TRUE)
  fig_preview <- c("\\documentclass[10pt,a4paper]{article}", "\\usepackage[margin=16mm]{geometry}",
    "\\usepackage{graphicx,lmodern}", "\\usepackage[T1]{fontenc}", "\\setlength{\\parindent}{0pt}", "\\begin{document}")
  for (f in names(plots$captions)) fig_preview <- c(fig_preview, paste0("\\noindent\\textbf{", esc(f), "}\\par\\medskip"),
    paste0("\\includegraphics[width=\\textwidth]{", f, "}\\par\\medskip"), plots$captions[[f]], "\\clearpage")
  writeLines(c(fig_preview, "\\end{document}"), file.path(stage, "figures", "figures_preview.tex"))
  writeLines(unlist(lapply(names(plots$captions), function(f) c(f, plots$captions[[f]], ""))), file.path(stage, "figures", "captions.txt"))
  utils::write.csv(x$config, file.path(stage, "provenance", "simulation_config.csv"), row.names = FALSE)
  utils::write.csv(x$seal$sources, file.path(stage, "provenance", "computation_sources.csv"), row.names = FALSE)
  utils::write.csv(x$seal$runtime$packages, file.path(stage, "provenance", "computation_packages.csv"), row.names = FALSE)
  script_md5 <- if (!is.na(spout_source)) unname(tools::md5sum(spout_source)) else NA_character_
  saveRDS(list(schema = "spectral-output-1", input_run_id = x$run_id, input_md5 = input_md5,
    exporter_md5 = script_md5, exporter_R = R.version.string, export_options = opts,
    computation_seal = x$seal, recovery_provenance = x$recovery_provenance), file.path(stage, "provenance", "output_provenance.rds"), version = 2L)
  writeLines(c(paste("Input run ID:", x$run_id), paste("Input MD5:", input_md5), paste("Exporter MD5:", script_md5),
    paste("Exporter R:", R.version.string), paste("Computation R:", x$seal$runtime$R),
    paste("Recovery kernel:", if (is.null(x$recovery_provenance)) "none recorded" else x$recovery_provenance$kernel),
    "Reproduction route: completed merged RDS -> this exporter -> LaTeX, PDF, CSV.",
    "This export validates saved numerical results; it does not rerun the estimators or validate a new cluster environment."),
    file.path(stage, "provenance", "provenance.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(stage, "provenance", "export_session_info.txt"))
  spout_assert(identical(unname(tools::md5sum(input)), input_md5), "Input RDS changed during export.")
  report <- c("SPECTRAL SIMULATION OUTPUT VALIDATION: PASS", paste("Input run ID:", x$run_id),
    "Production: 12 design/setting cells; 100 result keys; 108 table rows; N=1000; B=300.",
    "Stored point/interval dimensions and finite values: PASS", "Coverage recomputed from interval bounds: PASS",
    "Summary recomputed from point estimates and saved bootstrap summaries: PASS (tolerance 1e-10)",
    "Main/rank-sensitivity shared rows identical: 8", "Sparse simulation SE/CP remain unavailable: PASS",
    paste("Sparse alignment maximum change in C:", format(support$max_C_difference, scientific = TRUE)),
    paste("Sparse support definition: abs(aligned B) >", opts$`support-threshold`),
    paste("Recorded nonconverged point fits:", sum(v$summary$nonconverged_point_fits)),
    paste("Recorded projected sparse point fits (across result keys):", sum(v$summary$projected_point_fits)),
    paste("Recorded point warnings:", sum(v$summary$point_warning_count)),
    paste("Recorded bootstrap draws with issues:", sum(v$summary$bootstrap_draws_with_issues)),
    "Figure limits affect the displayed viewport only; all observations enter statistics.",
    "Original input RDS unchanged: PASS", "Full computation/recovery verification is inherited from the recorded strict merge.",
    "Scope excludes rank-selection experiments, tuning paths and real-data refits/bootstrap.")
  writeLines(report, file.path(stage, "VALIDATION.txt"))
  files <- list.files(stage, recursive = TRUE, full.names = FALSE)
  utils::write.csv(data.frame(path = files, md5 = unname(tools::md5sum(file.path(stage, files)))),
    file.path(stage, "OUTPUT_MD5.csv"), row.names = FALSE)
  spout_assert(file.rename(stage, opts$output), "Could not install output directory.")
  cat(paste(report, collapse = "\n"), "\nOutput:", opts$output, "\n")
  invisible(list(output = opts$output, summary = v$summary, support = support$summary))
}
if (sys.nframe() == 0L) spout_main()
