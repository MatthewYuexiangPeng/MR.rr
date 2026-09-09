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
    output = "paper/output/spectral_rebuild/manuscript_simulations_style_v2", `support-threshold` = "0.001",
    formats = "pdf,png", `png-dpi` = "600", `figure-palette` = "paper", `figure-range` = "central",
    `compare-with` = "")
  if (identical(args, "--help")) {
    cat("Usage: Rscript --vanilla paper/scripts/36_make_spectral_simulation_outputs.R\n",
      "  [--input=MERGED_RDS] [--output=NEW_OUTPUT_DIRECTORY] [--support-threshold=0.001]\n",
      "  [--formats=pdf,png|pdf] [--png-dpi=600] [--figure-palette=paper|legacy]\n",
      "  [--figure-range=central|full] [--compare-with=PREVIOUS_EXPORT_DIRECTORY]\n",
      "The output directory must not already exist. Requires only components included with R.\n")
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
  spout_assert(opts$formats %in% c("pdf", "pdf,png"), "Formats must be pdf or pdf,png.")
  opts$formats <- strsplit(opts$formats, ",", fixed = TRUE)[[1]]
  opts$`png-dpi` <- as.numeric(opts$`png-dpi`)
  spout_assert(opts$`png-dpi` %in% c(300, 600, 1200), "PNG DPI must be 300, 600 or 1200.")
  spout_assert(opts$`figure-palette` %in% c("paper", "legacy"), "Palette must be paper or legacy.")
  spout_assert(opts$`figure-range` %in% c("central", "full"), "Figure range must be central or full.")
  if ("png" %in% opts$formats) spout_assert(isTRUE(capabilities("png")), "This R build needs a PNG device, or use --formats=pdf.")
  if (nzchar(opts$`compare-with`)) spout_assert(dir.exists(opts$`compare-with`), "Previous export directory not found.")
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

# Figure presentation only. Table, prediction and sparse-alignment calculations
# above this section retain the validated v1 definitions.
spout_limits <- function(values, truth, mode = "central") {
  limits <- if (mode == "full") range(values, truth) else
    range(stats::quantile(values, c(.01, .99), names = FALSE), truth)
  if (diff(limits) == 0) limits <- limits + c(-1, 1)
  limits + c(-1, 1) * diff(limits) * .06
}
spout_box_record <- function(values, lim, figure, panel, keys) do.call(rbind, lapply(seq_along(values), function(i) {
  a <- values[[i]]; b <- grDevices::boxplot.stats(a)
  data.frame(figure = figure, panel = panel, result_key = keys[i], n = length(a),
    minimum = min(a), lower_whisker = b$stats[1], lower_hinge = b$stats[2], median = b$stats[3],
    upper_hinge = b$stats[4], upper_whisker = b$stats[5], maximum = max(a),
    axis_lower = lim[1], axis_upper = lim[2], outside_axis = sum(a < lim[1] | a > lim[2]))
}))

spout_display_B <- function(x, support) {
  # One fixed display convention, determined exclusively by TRUE B, as in freeze:
  # make the largest-magnitude entry in each true pathway positive, then order
  # pathways by that entry's exposure. All replicates already underwent the
  # validated signed-permutation alignment to the stored truth in spout_support.
  target <- x$truths$sparse_loading$B
  anchors <- apply(abs(target), 1L, which.max)
  signs <- sign(target[cbind(seq_len(nrow(target)), anchors)])
  signs[signs == 0] <- 1
  perm <- order(anchors, -apply(abs(target), 1L, max))
  transform <- function(b) (b * signs)[perm, , drop = FALSE]
  b <- support$aligned[["setting-4"]]
  shown <- matrix(NA_real_, nrow(b), ncol(b))
  for (i in seq_len(ncol(b))) shown[, i] <- as.vector(transform(matrix(b[, i], 2L, 9L)))
  # The fixed signed permutation has inverse equal to its transpose. Applying
  # it jointly to A preserves C for arbitrary factors, not just the truth.
  transform_matrix <- diag(signs)[perm, , drop = FALSE]
  spout_equal(crossprod(transform_matrix), diag(2), "Display transformation is not orthogonal.", 0)
  transformed_A <- x$truths$sparse_loading$A %*% t(transform_matrix)
  max_delta <- max(abs(transformed_A %*% transform(target) - x$truths$sparse_loading$C))
  spout_assert(max_delta < 1e-10, "B display convention changed true C.")
  spout_assert(identical(abs(shown), abs(b[as.vector(matrix(1:18, 2, 9)[perm, ]), , drop = FALSE])),
    "B display convention changed loading magnitudes.")
  list(values = shown, truth = transform(target), map = data.frame(display_pathway = 1:2,
    stored_pathway = perm, sign = signs[perm], anchor_exposure = anchors[perm]),
    max_true_C_difference = max_delta)
}

# grid is distributed with R. Explicit facet geometry avoids automatic base
# graphics text shrinking and retains the same Tukey box statistics as v1.
spout_draw_facets <- function(panels, rows, cols, row_labels, col_labels,
                              ylab, colors, xlabels = NULL, legend = TRUE,
                              letters = TRUE) {
  grid::grid.newpage()
  size <- grDevices::dev.size("in") * 72
  W <- size[1]; H <- size[2]
  is_dense <- cols == 9L
  left <- 38; right <- if (rows > 1L) 19 else 5
  bottom <- if (is.null(xlabels)) 14 else 85
  top <- 23; gap_x <- if (is_dense) 3.5 else 9; gap_y <- 13
  pw <- (W - left - right - (cols - 1L) * gap_x) / cols
  ph <- (H - top - bottom - (rows - 1L) * gap_y) / rows
  txt <- function(label, x, y, fontsize = 8, rot = 0, just = "centre", fontface = "plain", col = "#222222")
    grid::grid.text(label, x / W, y / H, rot = rot, just = just,
      gp = grid::gpar(fontsize = fontsize, fontfamily = "sans", fontface = fontface, col = col))
  rect <- function(x, y, w, h, fill, col = "#AAAAAA")
    grid::grid.rect(x = x / W, y = y / H, width = w / W, height = h / H,
      just = c("left", "bottom"), gp = grid::gpar(fill = fill, col = col, lwd = .55))
  line <- function(x, y, col = "#555555", lwd = .55)
    grid::grid.lines(x / W, y / H, gp = grid::gpar(col = col, lwd = lwd))
  for (r in seq_len(rows)) for (cc in seq_len(cols)) {
    k <- (r - 1L) * cols + cc; p <- panels[[k]]
    xx <- left + (cc - 1L) * (pw + gap_x)
    yy <- H - top - r * ph - (r - 1L) * gap_y
    n <- length(p$values); xr <- if (n == 1L) c(.4, 1.6) else c(.4, n + .6)
    ticks <- pretty(p$limits, n = 5L)
    ticks <- ticks[ticks >= p$limits[1] & ticks <= p$limits[2]]
    grid::pushViewport(grid::viewport(x = xx / W, y = yy / H, width = pw / W, height = ph / H,
      just = c("left", "bottom"), xscale = xr, yscale = p$limits, clip = "on"))
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    native_lines <- function(x, y, col, lwd = .6, lty = 1)
      grid::grid.lines(grid::unit(x, "native"), grid::unit(y, "native"),
        gp = grid::gpar(col = col, lwd = lwd, lty = lty))
    if (length(ticks) > 1L) {
      minor <- head(ticks, -1L) + diff(ticks) / 2
      for (t in minor) native_lines(xr, rep(t, 2L), "#F3F3F3", .45)
    }
    for (t in ticks) native_lines(xr, rep(t, 2L), "#E6E6E6", .55)
    if (n > 1L) for (j in seq_len(n)) native_lines(rep(j, 2L), p$limits, "#F0F0F0", .45)
    for (j in seq_len(n)) {
      b <- grDevices::boxplot.stats(p$values[[j]]); q <- b$stats
      box_width <- if (n == 1L) .72 else .7
      native_lines(rep(j, 2L), q[c(1, 5)], "#454545", .65)
      for (t in q[c(1, 5)]) native_lines(j + c(-.18, .18), rep(t, 2L), "#454545", .65)
      if (length(b$out)) grid::grid.points(x = grid::unit(rep(j, length(b$out)), "native"),
        y = grid::unit(b$out, "native"), pch = 16, size = grid::unit(.35, "mm"),
        gp = grid::gpar(col = "#777777"))
      grid::grid.rect(x = grid::unit(j, "native"), y = grid::unit(mean(q[c(2, 4)]), "native"),
        width = grid::unit(box_width, "native"), height = grid::unit(q[4] - q[2], "native"),
        gp = grid::gpar(fill = colors[j], col = "#353535", lwd = .65))
      native_lines(j + c(-box_width, box_width) / 2, rep(q[3], 2L), "#222222", 1)
    }
    native_lines(xr, rep(p$truth, 2L), if (n == 1L && p$truth == 0) "#737373" else "#C22E32", .95, 2)
    grid::popViewport()
    rect(xx, yy, pw, ph, NA)
    if (cc == 1L) {
      pos <- yy + (ticks - p$limits[1]) / diff(p$limits) * ph
      labels <- format(signif(ticks, 4), trim = TRUE, scientific = max(abs(ticks)) > 1e4)
      for (i in seq_along(ticks)) {
        line(c(xx - 2.3, xx), rep(pos[i], 2L))
        txt(labels[i], xx - 4.3, pos[i], 7.3, just = "right")
      }
    }
    if (r == 1L) {
      rect(xx, yy + ph, pw, 16, "#EEEEEE")
      txt(col_labels[cc], xx + pw / 2, yy + ph + 8, if (is_dense) 7 else 8.5)
    }
    if (cc == cols && rows > 1L) {
      rect(xx + pw, yy, 15, ph, "#EEEEEE")
      txt(row_labels[r], xx + pw + 7.5, yy + ph / 2, 8, rot = 270)
    }
    if (letters) {
      label <- if (k <= 26L) LETTERS[k] else paste0("A", LETTERS[k - 26L])
      txt(label, xx + 3, yy + ph - 3, 6.3, just = c("left", "top"), col = "#555555")
    }
    if (r == rows && !is.null(xlabels)) {
      pos <- xx + (seq_len(n) - xr[1]) / diff(xr) * pw
      for (j in seq_len(n)) {
        line(rep(pos[j], 2L), c(yy, yy - 2))
        txt(xlabels[j], pos[j], yy - 4, if (is_dense) 6.5 else 7.5,
          rot = 90, just = c("right", "centre"))
      }
    }
  }
  txt(ylab, 8, bottom + (H - top - bottom) / 2, 9, rot = 90)
  if (legend) {
    # One compact legend, leaving all nine exposure columns available to data.
    widths <- vapply(unname(spout_labels), function(s)
      grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = grid::gpar(fontsize = 7.2, fontfamily = "sans"))),
        "inches", valueOnly = TRUE) * 72 + 18, numeric(1))
    cursor <- (W - sum(widths)) / 2
    for (j in seq_along(widths)) {
      rect(cursor, 7, 7, 6, colors[j], "#444444")
      txt(spout_labels[j], cursor + 10, 10, 7.2, just = "left")
      cursor <- cursor + widths[j]
    }
  }
}

spout_figure_device <- function(out, stem, width, height, formats, dpi, draw) {
  for (fmt in formats) {
    file <- file.path(out, paste0(stem, ".", fmt))
    if (fmt == "pdf") {
      if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(file, width, height, family = "sans", bg = "white") else
        grDevices::pdf(file, width, height, family = "Helvetica", bg = "white", useDingbats = FALSE)
    } else grDevices::png(file, width = width, height = height, units = "in", res = dpi,
      type = if (isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType"), bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
  }
}

spout_figures <- function(x, support, out, options) {
  records <- pred_entries <- captions <- alt <- list(); index <- 0L
  palette <- if (options$`figure-palette` == "legacy")
    c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7") else spout_colors
  formats <- options$formats; mode <- options$`figure-range`; dpi <- options$`png-dpi`
  record_panels <- function(panels, stem) {
    for (p in panels) {
      index <<- index + 1L
      records[[index]] <<- spout_box_record(p$values, p$limits, paste0(stem, ".pdf"), p$id, p$keys)
    }
  }
  window_note <- if (mode == "central") paste("The display window uses pooled 1st and 99th percentiles,",
    "includes the true values and adds a 6% margin; observations outside this window are clipped only in the display.") else
    "The display window includes the full range of estimates and true values, with a 6% margin."
  box_note <- paste("Boxes show Tukey hinges and medians; whiskers extend to observations within 1.5 interquartile ranges,",
    "and points show outliers. All 1,000 Monte Carlo replicates enter every boxplot and numerical summary.")
  for (design in c("generic", "sparse_loading")) for (s in c(4L, 1L)) {
    keys <- vapply(names(spout_labels), function(m) spout_key(design, s, m), character(1))
    points <- lapply(keys, function(k) x$results[[k]]$point); truth <- x$truths[[design]]$C
    setting <- sprintf("%s, Setting %d (SIV = %.2f).", if (design == "generic") "Generic low-rank design" else
      "Sparse-loading design", s, x$parameters[[paste(design, s, sep = "/")]]$iv_strength)
    stem <- sprintf("figure_C_%s_setting%d", design, s); panels <- list()
    for (y in 1:3) {
      rows <- seq(y, 27, by = 3)
      lim <- spout_limits(unlist(lapply(points, function(p) p[rows, ])), truth[y, ], mode)
      for (j in 1:9) panels[[length(panels) + 1L]] <- list(values = lapply(points, function(p) p[y + 3 * (j - 1), ]),
        limits = lim, truth = truth[y, j], id = paste(y, j, sep = "/"), keys = keys)
    }
    record_panels(panels, stem)
    spout_figure_device(out, stem, 7.2, 5.35, formats, dpi, function()
      spout_draw_facets(panels, 3, 9, paste("Outcome", 1:3), paste("Exposure", 1:9),
        "Estimated C entry", palette, unname(spout_labels)))
    captions[[paste0(stem, ".pdf")]] <- paste(setting, "Entrywise estimates of C. Rows represent outcomes and columns exposures.",
      "Seven estimators appear in the same left-to-right order in every facet. Red dashed lines indicate the true entries.",
      "A common vertical scale is used across exposures within each outcome.", window_note, box_note)
    alt[[paste0(stem, ".pdf")]] <- paste("A three-by-nine grid of boxplots for three outcomes and nine exposures,",
      "with seven labelled estimators per facet and a dashed line marking the corresponding true effect.", setting)
    pred <- lapply(points, spout_predict, exposure = x$prediction_exposure)
    target <- as.vector(spout_predict(matrix(as.vector(truth), 27, 1), x$prediction_exposure))
    for (m in seq_along(pred)) for (y in 1:3) {
      v <- pred[[m]][y, ]; err <- v - target[y]
      pred_entries[[length(pred_entries) + 1L]] <- data.frame(design = design, setting = s,
        result_key = keys[m], method = names(spout_labels)[m], outcome = y, truth = target[y],
        mean = mean(v), signed_bias = mean(err), absolute_bias = abs(mean(err)), sd = stats::sd(v),
        rmse = sqrt(mean(err^2)), median = stats::median(v), q25 = stats::quantile(v, .25), q75 = stats::quantile(v, .75))
    }
    stem <- sprintf("figure_prediction_%s_setting%d", design, s)
    lim <- spout_limits(unlist(pred), target, mode)
    panels <- lapply(1:3, function(y) list(values = lapply(pred, function(p) p[y, ]),
      limits = lim, truth = target[y], id = as.character(y), keys = keys))
    record_panels(panels, stem)
    spout_figure_device(out, stem, 7.2, 3.9, formats, dpi, function()
      spout_draw_facets(panels, 1, 3, NULL, paste("Outcome", 1:3), "Estimated risk score",
        palette, unname(spout_labels)))
    captions[[paste0(stem, ".pdf")]] <- paste(setting, "Predicted risk scores from the stored fixed exposure vector.",
      "Red dashed lines indicate the true Cx. A common vertical scale is used for the three outcomes within this figure.", window_note, box_note)
    alt[[paste0(stem, ".pdf")]] <- paste("Three side-by-side boxplot panels compare seven labelled estimators' risk scores",
      "on a shared vertical scale; dashed lines mark the true scores for the fixed exposure profile.", setting)
  }
  display <- spout_display_B(x, support); truth <- display$truth; bm <- display$values
  stem <- "figure_B_sparse_loading_setting4"; panels <- list()
  lim <- spout_limits(bm, truth, mode)
  for (k in 1:2) for (j in 1:9) panels[[length(panels) + 1L]] <- list(
    values = list(bm[k + 2 * (j - 1), ]), limits = lim, truth = truth[k, j],
    id = paste(k, j, sep = "/"), keys = spout_key("sparse_loading", 4, "sparse_mr_rr"))
  record_panels(panels, stem)
  spout_figure_device(out, stem, 7.2, 3.3, formats, dpi, function()
    spout_draw_facets(panels, 2, 9, paste("Pathway", 1:2), paste("Exposure", 1:9), "Estimated B entry",
      if (options$`figure-palette` == "legacy") "#90EE90" else palette[6], legend = FALSE))
  captions[[paste0(stem, ".pdf")]] <- paste("Sparse-loading design, Setting 4. Estimated B is aligned to the stored truth over all signed row permutations,",
    "with A transformed jointly to preserve C. For display, a single fixed signed row permutation makes the largest absolute true loading",
    "in each pathway positive and orders pathways by its exposure, matching the convention in the archived plotting code.",
    "Red dashed lines mark nonzero true loadings; gray dashed lines mark zero true loadings. Both pathways share one vertical scale.", window_note, box_note)
  alt[[paste0(stem, ".pdf")]] <- paste("A two-by-nine grid shows aligned sparse loading estimates for two pathways and nine exposures.",
    "Dashed red or gray lines mark nonzero or zero true loadings, respectively. Pathway orientation is fixed by the true loading matrix.")
  manuscript_order <- c("figure_C_generic_setting4.pdf", "figure_C_generic_setting1.pdf",
    "figure_C_sparse_loading_setting4.pdf", "figure_C_sparse_loading_setting1.pdf",
    "figure_B_sparse_loading_setting4.pdf", "figure_prediction_generic_setting4.pdf",
    "figure_prediction_generic_setting1.pdf", "figure_prediction_sparse_loading_setting4.pdf",
    "figure_prediction_sparse_loading_setting1.pdf")
  list(boxplots = do.call(rbind, records), prediction = do.call(rbind, pred_entries),
    captions = captions[manuscript_order], alt = alt[manuscript_order], display_B = display, palette = palette)
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
  plots <- spout_figures(x, support, file.path(stage, "figures"), opts)
  write_csv <- function(object, name) utils::write.csv(object, file.path(stage, "statistics", name), row.names = FALSE, na = "NA")
  write_csv(v$summary, "simulation_summary.csv"); write_csv(v$table_rows, "simulation_table_rows.csv")
  write_csv(v$entrywise, "simulation_entrywise.csv"); write_csv(plots$prediction, "prediction_summary.csv")
  write_csv(plots$display_B$map, "figure_B_display_map.csv")
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
  esc <- function(s) gsub("%", "\\%", gsub("_", "\\_", s, fixed = TRUE), fixed = TRUE)
  fig_preview <- c("\\documentclass[10pt,a4paper]{article}", "\\usepackage[margin=16mm]{geometry}",
    "\\usepackage{graphicx,lmodern}", "\\usepackage[T1]{fontenc}", "\\setlength{\\parindent}{0pt}", "\\begin{document}")
  for (f in names(plots$captions)) fig_preview <- c(fig_preview, paste0("\\noindent\\textbf{", esc(f), "}\\par\\medskip"),
    paste0("\\includegraphics[width=\\textwidth]{", f, "}\\par\\medskip"), esc(plots$captions[[f]]), "\\clearpage")
  writeLines(c(fig_preview, "\\end{document}"), file.path(stage, "figures", "figures_preview.tex"))
  writeLines(unlist(lapply(names(plots$captions), function(f) c(f, plots$captions[[f]], ""))), file.path(stage, "figures", "captions.txt"))
  writeLines(unlist(lapply(names(plots$alt), function(f) c(f, paste("Alt text:", plots$alt[[f]]), ""))),
    file.path(stage, "figures", "alt_text.txt"))
  utils::write.csv(data.frame(method = names(spout_labels), label = unname(spout_labels), color = plots$palette),
    file.path(stage, "figures", "method_colors.csv"), row.names = FALSE)
  utils::write.csv(data.frame(figure = names(plots$captions), width_inches = 7.2,
    height_inches = ifelse(grepl("figure_C_", names(plots$captions)), 5.35,
      ifelse(grepl("figure_B_", names(plots$captions)), 3.3, 3.9)),
    pdf_device = if (isTRUE(capabilities("cairo"))) "cairo_pdf (embedded fonts)" else "pdf (standard Helvetica)",
    png_dpi = if ("png" %in% opts$formats) opts$`png-dpi` else NA, palette = opts$`figure-palette`,
    viewport = opts$`figure-range`), file.path(stage, "figures", "figure_manifest.csv"), row.names = FALSE)
  comparison <- character()
  if (nzchar(opts$`compare-with`)) {
    old <- opts$`compare-with`
    csv <- c("simulation_summary.csv", "simulation_table_rows.csv", "simulation_entrywise.csv",
      "prediction_summary.csv", "sparse_support_summary.csv", "sparse_support_replicates.csv",
      "sparse_B_entrywise.csv", "prediction_exposure.csv", "shared_rank2_rows.csv")
    previous_provenance <- readRDS(file.path(old, "provenance", "output_provenance.rds"))
    spout_assert(identical(previous_provenance$input_md5, input_md5), "Previous export used a different input RDS.")
    for (f in csv) spout_equal(utils::read.csv(file.path(stage, "statistics", f)),
      utils::read.csv(file.path(old, "statistics", f)), paste("Restyling changed numerical output:", f))
    for (f in tables) spout_assert(identical(readLines(file.path(stage, "tables", f)),
      readLines(file.path(old, "tables", f))), paste("Restyling changed table:", f))
    spout_equal(readRDS(file.path(stage, "statistics", "sparse_B_aligned.rds")),
      readRDS(file.path(old, "statistics", "sparse_B_aligned.rds")), "Stored sparse alignment changed.", 0)
    comparison <- c("Compared with previous export: PASS (same input MD5)",
      "All four LaTeX tables unchanged: PASS", "Nine numerical CSV outputs unchanged: PASS (tolerance 1e-10)",
      "Stored aligned sparse B unchanged: PASS (tolerance 0)",
      "B figure uses a fixed recorded signed permutation; stored loading summaries retain their original orientation.")
    writeLines(comparison, file.path(stage, "STYLE_COMPARISON.txt"))
  }
  utils::write.csv(x$config, file.path(stage, "provenance", "simulation_config.csv"), row.names = FALSE)
  utils::write.csv(x$seal$sources, file.path(stage, "provenance", "computation_sources.csv"), row.names = FALSE)
  utils::write.csv(x$seal$runtime$packages, file.path(stage, "provenance", "computation_packages.csv"), row.names = FALSE)
  script_md5 <- if (!is.na(spout_source)) unname(tools::md5sum(spout_source)) else NA_character_
  saveRDS(list(schema = "spectral-output-2", presentation_version = "facets-v2", input_run_id = x$run_id, input_md5 = input_md5,
    exporter_md5 = script_md5, exporter_R = R.version.string, export_options = opts, B_display_map = plots$display_B$map,
    computation_seal = x$seal, recovery_provenance = x$recovery_provenance), file.path(stage, "provenance", "output_provenance.rds"), version = 2L)
  writeLines(c(paste("Input run ID:", x$run_id), paste("Input MD5:", input_md5), paste("Exporter MD5:", script_md5),
    paste("Exporter R:", R.version.string), paste("Computation R:", x$seal$runtime$R),
    paste("Recovery kernel:", if (is.null(x$recovery_provenance)) "none recorded" else x$recovery_provenance$kernel),
    "Reproduction route: completed merged RDS -> this exporter -> LaTeX, PDF/PNG, CSV.",
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
    paste("Figure formats:", paste(opts$formats, collapse = ", ")),
    paste("Palette:", opts$`figure-palette`, "; viewport:", opts$`figure-range`),
    "B display orientation follows the archived truth-based sign/order convention; stored fits and alignment unchanged.",
    comparison,
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
