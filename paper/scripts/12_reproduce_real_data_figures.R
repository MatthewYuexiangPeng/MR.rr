#!/usr/bin/env Rscript

# Reproduce manuscript Figure 2 and supplementary Figures S2 and S12 from
# the frozen real-data inputs and archived bootstrap results. Generated files
# are written to paper/output/figures/. The frozen project is read-only.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    freeze_root <- file.path(
      current,
      "freeze",
      "current_analysis_20260823",
      "project"
    )

    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(freeze_root) &&
        dir.exists(file.path(current, "paper", "scripts"))) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      stop(
        "Could not locate the MR.rr repository root.",
        call. = FALSE
      )
    }

    current <- parent
  }
}


require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]

  if (length(missing) > 0L) {
    stop(
      paste(
        "Required files are missing:",
        paste(missing, collapse = "\n  "),
        sep = "\n  "
      ),
      call. = FALSE
    )
  }

  invisible(paths)
}


require_packages <- function(packages) {
  installed <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  if (any(!installed)) {
    stop(
      paste0(
        "Missing packages: ",
        paste(packages[!installed], collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  invisible(packages)
}


reconstruct_real_data <- function(data_path, correlation_path) {
  real_data <- utils::read.csv(data_path)
  correlation_data <- utils::read.csv(correlation_path)

  required_columns <- c(
    "ImpMAF",
    paste0("gamma_exp", 1:9),
    paste0("se_exp", 1:9),
    paste0("gamma_out", 1:4),
    paste0("se_out", 1:4)
  )
  missing_columns <- setdiff(required_columns, names(real_data))

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "The frozen real-data file is missing columns: ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  variant_variance <-
    2 * real_data$ImpMAF * (1 - real_data$ImpMAF)
  scale_factor <- sqrt(variant_variance)

  gamma_exp <- sweep(
    as.matrix(real_data[, paste0("gamma_exp", 1:9)]),
    1L,
    scale_factor,
    "*"
  )
  se_exp <- sweep(
    as.matrix(real_data[, paste0("se_exp", 1:9)]),
    1L,
    scale_factor,
    "*"
  )
  gamma_out <- sweep(
    as.matrix(real_data[, paste0("gamma_out", 1:4)]),
    1L,
    scale_factor,
    "*"
  )
  se_out <- sweep(
    as.matrix(real_data[, paste0("se_out", 1:4)]),
    1L,
    scale_factor,
    "*"
  )

  correlation_exp <- as.matrix(
    correlation_data[1:9, 1:9]
  )

  # The primary paper analysis removes the first outcome.
  gamma_out <- gamma_out[, 2:4, drop = FALSE]
  se_out <- se_out[, 2:4, drop = FALSE]
  correlation_out <- as.matrix(
    correlation_data[11:13, 11:13]
  )

  if (!identical(
    c(nrow(gamma_exp), ncol(gamma_exp), ncol(gamma_out)),
    c(177L, 9L, 3L)
  )) {
    stop(
      "The frozen real-data inputs have unexpected dimensions.",
      call. = FALSE
    )
  }

  n_instruments <- nrow(gamma_exp)
  n_exposures <- ncol(gamma_exp)
  n_outcomes <- ncol(gamma_out)
  Sigma_X_sum <- matrix(0, n_exposures, n_exposures)
  Sigma_Y_sum <- matrix(0, n_outcomes, n_outcomes)

  for (instrument in seq_len(n_instruments)) {
    D_exp <- diag(se_exp[instrument, ])
    D_out <- diag(se_out[instrument, ])

    Sigma_X_sum <-
      Sigma_X_sum +
      D_exp %*% correlation_exp %*% D_exp
    Sigma_Y_sum <-
      Sigma_Y_sum +
      D_out %*% correlation_out %*% D_out
  }

  Sigma_X <- Sigma_X_sum / n_instruments
  Sigma_Y <- Sigma_Y_sum / n_instruments

  list(
    gamma_exp = gamma_exp,
    gamma_out = gamma_out,
    Sigma_X = Sigma_X,
    Sigma_Y = Sigma_Y,
    W = solve(Sigma_Y)
  )
}


compute_eta_path <- function(
    rank,
    eta_grid,
    inputs,
    legacy,
    zero_tolerance = 1e-2) {
  path_rows <- vector("list", length(eta_grid))
  W_sqrt <- legacy$.sqrt_matrix(inputs$W)

  for (eta_index in seq_along(eta_grid)) {
    eta <- eta_grid[[eta_index]]

    fit <- legacy$mr_rr_sparse(
      GAMMA_hat = inputs$gamma_out,
      gamma_hat = inputs$gamma_exp,
      W = inputs$W,
      Sigma_X = inputs$Sigma_X,
      lambda = rep(eta, ncol(inputs$gamma_exp)),
      r = rank,
      max_iter = 100L,
      tol = 1e-2
    )

    B_hat <- if (!is.null(fit$B_raw)) {
      fit$B_raw
    } else {
      fit$B
    }
    C_hat <- if (!is.null(fit$AB_raw)) {
      fit$AB_raw
    } else {
      fit$A %*% B_hat
    }

    residual_term <- norm(
      (
        inputs$gamma_out -
          inputs$gamma_exp %*% t(C_hat)
      ) %*% W_sqrt,
      type = "F"
    )^2 / nrow(inputs$gamma_exp)

    debias_term <- sum(diag(
      W_sqrt %*%
        C_hat %*%
        inputs$Sigma_X %*%
        t(C_hat) %*%
        W_sqrt
    ))

    path_rows[[eta_index]] <- data.frame(
      Rank = rank,
      Eta = eta,
      Debiased_objective = residual_term - debias_term,
      Zero_proportion = mean(abs(B_hat) < zero_tolerance),
      Nonzero_count = sum(abs(B_hat) >= zero_tolerance),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  path <- do.call(rbind, path_rows)
  rownames(path) <- NULL

  if (nrow(path) != length(eta_grid) ||
      any(!is.finite(as.matrix(path)))) {
    stop(
      "A sparse eta path is incomplete or non-finite.",
      call. = FALSE
    )
  }

  path
}


eta_label <- function(x) {
  label <- formatC(x, format = "e", digits = 1)
  label <- sub("\\.0e", "e", label)
  sub("e-0", "e-", label)
}


make_eta_plot <- function(
    path,
    rank,
    selected_eta,
    eta_grid) {
  objective_plot <- ggplot2::ggplot(
    path,
    ggplot2::aes(x = Eta, y = Debiased_objective)
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
      xintercept = selected_eta,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::scale_x_log10(
      breaks = eta_grid,
      labels = NULL
    ) +
    ggplot2::labs(
      title = paste0("Rank ", rank),
      x = NULL,
      y = "Debiased objective"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  sparsity_plot <- ggplot2::ggplot(
    path,
    ggplot2::aes(x = Eta, y = Zero_proportion)
  ) +
    ggplot2::geom_step(
      direction = "hv",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
      xintercept = selected_eta,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::scale_x_log10(
      breaks = eta_grid,
      labels = eta_label
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent()
    ) +
    ggplot2::labs(
      x = expression(eta),
      y = "Zero proportion in B"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      )
    )

  patchwork::wrap_plots(
    objective_plot,
    sparsity_plot,
    ncol = 1,
    heights = c(1, 1)
  )
}


make_halfeye_data <- function(result, bootstrap) {
  method_keys <- c(
    "ivw",
    "adivw",
    "naive_MRrr",
    "MRrr",
    "MRrr_regularized",
    "MRrr_sparse",
    "Mr_DAG"
  )
  method_labels <- c(
    "IVW",
    "SRIVW",
    "Naive MR-rr",
    "MR-rr",
    "Reg. MR-rr",
    "Sparse MR-rr",
    "MrDAG"
  )
  names(method_labels) <- method_keys

  estimates <- list(
    ivw = result$C_ivw,
    adivw = result$C_adivw,
    naive_MRrr = result$C_naive_MRrr,
    MRrr = result$C_MRrr,
    MRrr_regularized = result$C_MRrr_regularized,
    MRrr_sparse = result$C_MRrr_sparse,
    Mr_DAG = result$C_Mr_DAG
  )

  outcome_names <- c("LAS", "CES", "SVS")
  exposure_names <- c(
    "MMP12",
    "CNTN1",
    "FGL1",
    "MXRA8",
    "CNTFR",
    "SCG3",
    "HTRA1",
    "CLEC3B",
    "ANTXR2"
  )
  long_rows <- list()
  interval_rows <- list()
  long_index <- 1L
  interval_index <- 1L

  for (outcome_index in seq_along(outcome_names)) {
    for (exposure_index in seq_along(exposure_names)) {
      for (method in method_keys) {
        values <- bootstrap[[method]][
          outcome_index,
          exposure_index,
        ]
        values <- values[is.finite(values)]

        if (length(values) == 0L) {
          stop(
            "A real-data bootstrap cell has no finite values.",
            call. = FALSE
          )
        }

        limits <- stats::quantile(
          values,
          probs = c(0.025, 0.975),
          names = FALSE,
          type = 7
        )

        long_rows[[long_index]] <- data.frame(
          Outcome = outcome_names[[outcome_index]],
          Exposure = exposure_names[[exposure_index]],
          Estimator = unname(method_labels[[method]]),
          Value = values,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        interval_rows[[interval_index]] <- data.frame(
          Outcome = outcome_names[[outcome_index]],
          Exposure = exposure_names[[exposure_index]],
          Estimator = unname(method_labels[[method]]),
          Point_estimate = estimates[[method]][
            outcome_index,
            exposure_index
          ],
          Lower = limits[[1L]],
          Upper = limits[[2L]],
          Significant = limits[[1L]] > 0 || limits[[2L]] < 0,
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        long_index <- long_index + 1L
        interval_index <- interval_index + 1L
      }
    }
  }

  long <- do.call(rbind, long_rows)
  intervals <- do.call(rbind, interval_rows)

  long$Outcome <- factor(long$Outcome, levels = outcome_names)
  long$Exposure <- factor(long$Exposure, levels = exposure_names)
  long$Estimator <- factor(
    long$Estimator,
    levels = unname(method_labels)
  )
  intervals$Outcome <- factor(
    intervals$Outcome,
    levels = outcome_names
  )
  intervals$Exposure <- factor(
    intervals$Exposure,
    levels = exposure_names
  )
  intervals$Estimator <- factor(
    intervals$Estimator,
    levels = unname(method_labels)
  )

  list(
    long = long,
    intervals = intervals
  )
}


make_halfeye_plot <- function(plot_data) {
  ggplot2::ggplot(
    plot_data$long,
    ggplot2::aes(
      x = Estimator,
      y = Value,
      fill = Estimator
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35
    ) +
    ggdist::stat_halfeye(
      adjust = 1,
      .width = 0.95,
      justification = -0.2,
      scale = 0.8,
      point_colour = NA,
      interval_colour = "black",
      interval_size = 0.75,
      slab_alpha = 0.7,
      normalize = "xy"
    ) +
    ggplot2::geom_point(
      data = plot_data$intervals,
      inherit.aes = FALSE,
      ggplot2::aes(
        x = Estimator,
        y = Point_estimate,
        color = Significant
      ),
      shape = 18,
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "black", `TRUE` = "red"),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.5, 0.5)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Outcome),
      cols = ggplot2::vars(Exposure)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "",
      y = "Estimated C Value",
      fill = "Estimator"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
}


write_png <- function(plot, path, width, height) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop(
      "A figure file was not created: ",
      path,
      call. = FALSE
    )
  }

  cat("Wrote:", path, "\n")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}


repo_root <- locate_repo_root()
setwd(repo_root)

require_packages(c(
  "ggplot2",
  "ggdist",
  "patchwork",
  "scales"
))

scripts_root <- file.path(repo_root, "paper", "scripts")
validation_script <- file.path(
  scripts_root,
  "09_smoke_test_real_data_manuscript.R"
)
freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)
estimator_path <- file.path(
  freeze_root,
  "scripts",
  "MR_rr_estimators.R"
)
data_path <- file.path(
  freeze_root,
  "data",
  "dat_1e-4.csv"
)
correlation_path <- file.path(
  freeze_root,
  "data",
  "rho_mat_1e-4.csv"
)

require_files(c(
  validation_script,
  estimator_path,
  data_path,
  correlation_path
))

cat("Loading and validating archived real-data results.\n")
verified <- new.env(parent = globalenv())
sys.source(
  validation_script,
  envir = verified
)

cat("\nBuilding manuscript Figure 2.\n")
rank1_plot_data <- make_halfeye_data(
  verified$rank1_result,
  verified$rank1_bootstrap
)
figure_02 <- make_halfeye_plot(rank1_plot_data)

cat("Building supplementary Figure S12.\n")
rank2_plot_data <- make_halfeye_data(
  verified$rank2_result,
  verified$rank2_bootstrap
)
figure_S12 <- make_halfeye_plot(rank2_plot_data)

cat("Recomputing sparse eta paths for supplementary Figure S2.\n")
legacy <- new.env(parent = globalenv())
sys.source(
  estimator_path,
  envir = legacy
)
real_inputs <- reconstruct_real_data(
  data_path,
  correlation_path
)
eta_grid <- c(
  1e-4,
  3e-4,
  5e-4,
  7e-4,
  1e-3,
  1.2e-3,
  1.5e-3,
  2e-3,
  2.5e-3,
  3e-3,
  4e-3,
  5e-3,
  1e-2
)

set.seed(123)
rank1_eta_path <- compute_eta_path(
  rank = 1L,
  eta_grid = eta_grid,
  inputs = real_inputs,
  legacy = legacy
)
rank2_eta_path <- compute_eta_path(
  rank = 2L,
  eta_grid = eta_grid,
  inputs = real_inputs,
  legacy = legacy
)

if (abs(verified$rank1_result$sparse_eta - 1.2e-3) > 1e-15 ||
    abs(verified$rank2_result$sparse_eta - 1e-3) > 1e-15) {
  stop(
    "The selected sparse eta values differ from the manuscript.",
    call. = FALSE
  )
}

figure_S02 <- patchwork::wrap_plots(
  make_eta_plot(
    rank1_eta_path,
    rank = 1L,
    selected_eta = verified$rank1_result$sparse_eta,
    eta_grid = eta_grid
  ),
  make_eta_plot(
    rank2_eta_path,
    rank = 2L,
    selected_eta = verified$rank2_result$sparse_eta,
    eta_grid = eta_grid
  ),
  ncol = 2
)

output_root <- file.path(
  repo_root,
  "paper",
  "output",
  "figures"
)
dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_paths <- c(
  write_png(
    figure_02,
    file.path(output_root, "Figure_02_real_data_rank1.png"),
    width = 13,
    height = 6
  ),
  write_png(
    figure_S02,
    file.path(output_root, "Figure_S02_sparse_eta_paths.png"),
    width = 14,
    height = 7.5
  ),
  write_png(
    figure_S12,
    file.path(output_root, "Figure_S12_real_data_rank2.png"),
    width = 13,
    height = 6
  )
)

eta_path_file <- file.path(
  output_root,
  "Figure_S02_sparse_eta_paths.csv"
)
utils::write.csv(
  rbind(rank1_eta_path, rank2_eta_path),
  eta_path_file,
  row.names = FALSE
)
cat("Wrote:", eta_path_file, "\n")

generated_paths <- c(
  figure_paths,
  normalizePath(
    eta_path_file,
    winslash = "/",
    mustWork = TRUE
  )
)
manifest <- data.frame(
  File = basename(generated_paths),
  Bytes = as.numeric(file.info(generated_paths)$size),
  MD5 = unname(tools::md5sum(generated_paths)),
  stringsAsFactors = FALSE,
  row.names = NULL
)
manifest_path <- file.path(output_root, "manifest.csv")
utils::write.csv(
  manifest,
  manifest_path,
  row.names = FALSE
)
cat("Wrote:", manifest_path, "\n")

if (nrow(rank1_plot_data$intervals) != 189L ||
    nrow(rank2_plot_data$intervals) != 189L ||
    nrow(rank1_eta_path) != 13L ||
    nrow(rank2_eta_path) != 13L ||
    any(manifest$Bytes <= 0) ||
    any(!nzchar(manifest$MD5))) {
  stop(
    "The generated real-data figure inventory is incomplete.",
    call. = FALSE
  )
}

cat("\nGenerated real-data figures:", length(figure_paths), "\n")
cat("Generated sparse eta path rows: 26\n")
cat("Manuscript Figure 2 and supplementary Figures S2/S12: PASS\n")
