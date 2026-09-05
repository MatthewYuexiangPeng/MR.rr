#!/usr/bin/env Rscript

# Reproduce supplementary simulation Figures S1 and S3-S11 from the frozen
# Monte Carlo result objects. Generated files are written to
# paper/output/figures/. The frozen project is read-only.

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


load_required_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)

  if (!object_name %in% loaded) {
    stop(
      basename(path),
      " does not contain `",
      object_name,
      "`.",
      call. = FALSE
    )
  }

  environment[[object_name]]
}


validate_matrix <- function(
    matrix,
    expected_rows,
    expected_columns,
    label,
    allow_na = TRUE) {
  expected_dimension <- as.integer(c(
    expected_rows,
    expected_columns
  ))

  if (!identical(dim(matrix), expected_dimension)) {
    stop(
      label,
      " has unexpected dimensions.",
      call. = FALSE
    )
  }

  if (any(is.infinite(matrix)) ||
      (!allow_na && anyNA(matrix))) {
    stop(
      label,
      " contains an invalid value.",
      call. = FALSE
    )
  }

  invisible(matrix)
}


validate_prediction_reconstruction <- function(
    main_result,
    sparse_result,
    fixed_exposure,
    label,
    tolerance = 1e-10) {
  component_pairs <- data.frame(
    effect = c(
      "C_ivw_list",
      "C_adivw_list",
      "AB_list",
      "AB_d_list",
      "AB_d_r_list",
      "MrDAG_list"
    ),
    prediction = c(
      "Y_pred_C_ivw_list",
      "Y_pred_C_adivw_list",
      "Y_pred_AB_list",
      "Y_pred_AB_d_list",
      "Y_pred_AB_d_r_list",
      "Y_pred_MrDAG_list"
    ),
    stringsAsFactors = FALSE
  )
  scenario_names <- names(main_result$AB_list)
  reconstruct_prediction <- function(effect_matrix) {
    do.call(
      rbind,
      lapply(
        seq_len(3L),
        function(outcome) {
          effect_rows <- seq(
            outcome,
            nrow(effect_matrix),
            by = 3L
          )
          as.vector(crossprod(
            fixed_exposure,
            effect_matrix[effect_rows, , drop = FALSE]
          ))
        }
      )
    )
  }
  maximum_difference <- 0

  for (scenario in scenario_names) {
    for (pair_index in seq_len(nrow(component_pairs))) {
      effect_component <- component_pairs$effect[[pair_index]]
      prediction_component <-
        component_pairs$prediction[[pair_index]]
      effect_matrix <- main_result[[effect_component]][[scenario]]
      archived_prediction <-
        main_result[[prediction_component]][[scenario]]

      validate_matrix(
        effect_matrix,
        expected_rows = 27L,
        expected_columns = 1000L,
        label = paste(
          label,
          scenario,
          effect_component
        )
      )
      validate_matrix(
        archived_prediction,
        expected_rows = 3L,
        expected_columns = 1000L,
        label = paste(
          label,
          scenario,
          prediction_component
        )
      )

      reconstructed_prediction <-
        reconstruct_prediction(effect_matrix)
      finite_mask <-
        is.finite(reconstructed_prediction) &
        is.finite(archived_prediction)

      if (!identical(
        is.finite(reconstructed_prediction),
        is.finite(archived_prediction)
      ) || !any(finite_mask)) {
        stop(
          label,
          " prediction finite-value masks differ for ",
          scenario,
          ".",
          call. = FALSE
        )
      }

      difference <- max(abs(
        reconstructed_prediction[finite_mask] -
          archived_prediction[finite_mask]
      ))
      maximum_difference <- max(maximum_difference, difference)
    }

    sparse_effect <- sparse_result$C_sparse_list[[scenario]]
    sparse_prediction <-
      sparse_result$Y_pred_C_sparse_list[[scenario]]

    validate_matrix(
      sparse_effect,
      expected_rows = 27L,
      expected_columns = 1000L,
      label = paste(label, scenario, "C_sparse_list")
    )
    validate_matrix(
      sparse_prediction,
      expected_rows = 3L,
      expected_columns = 1000L,
      label = paste(label, scenario, "Y_pred_C_sparse_list")
    )

    reconstructed_sparse <-
      reconstruct_prediction(sparse_effect)
    finite_mask <-
      is.finite(reconstructed_sparse) &
      is.finite(sparse_prediction)

    if (!identical(
      is.finite(reconstructed_sparse),
      is.finite(sparse_prediction)
    ) || !any(finite_mask)) {
      stop(
        label,
        " sparse prediction finite-value masks differ for ",
        scenario,
        ".",
        call. = FALSE
      )
    }

    difference <- max(abs(
      reconstructed_sparse[finite_mask] -
        sparse_prediction[finite_mask]
    ))
    maximum_difference <- max(maximum_difference, difference)
  }

  if (!is.finite(maximum_difference) ||
      maximum_difference > tolerance) {
    stop(
      label,
      " archived risk scores differ from C-hat times X by ",
      format(maximum_difference, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  maximum_difference
}


make_eta_figure <- function(eta_summary) {
  required_columns <- c(
    "sim_name",
    "eta",
    "objective",
    "zero_prop"
  )
  missing_columns <- setdiff(
    required_columns,
    names(eta_summary)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "The archived eta summary is missing columns: ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  scenario_levels <- c(
    "me_1_effect_0.25",
    "me_1_effect_1",
    "me_2.5_effect_0.25",
    "me_2.5_effect_1"
  )
  eta_summary$sim_name <- factor(
    eta_summary$sim_name,
    levels = scenario_levels
  )

  if (nrow(eta_summary) != 124L ||
      anyNA(eta_summary$sim_name) ||
      any(!is.finite(eta_summary$eta)) ||
      any(!is.finite(eta_summary$objective)) ||
      any(!is.finite(eta_summary$zero_prop))) {
    stop(
      "The archived eta path is incomplete or non-finite.",
      call. = FALSE
    )
  }

  selected_eta <- c(
    me_1_effect_0.25 = 1e-3,
    me_1_effect_1 = 1e-3,
    me_2.5_effect_0.25 = 1e-4,
    me_2.5_effect_1 = 1e-3
  )

  # The archived logarithmic grid explicitly included 1e-3 but not 1e-4.
  # Setting 1 used 1e-4 as a visually selected value on the continuous path,
  # so the selected value must lie within the displayed range but need not be
  # an exact grid point.
  for (scenario in names(selected_eta)) {
    scenario_path <- eta_summary[
      eta_summary$sim_name == scenario,
      ,
      drop = FALSE
    ]

    eta_range <- range(scenario_path$eta)

    if (selected_eta[[scenario]] < eta_range[[1L]] ||
        selected_eta[[scenario]] > eta_range[[2L]]) {
      stop(
        "The selected eta lies outside the path for ",
        scenario,
        ".",
        call. = FALSE
      )
    }
  }

  objective_plot <- ggplot2::ggplot(
    eta_summary,
    ggplot2::aes(x = eta, y = objective)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(
      ggplot2::vars(sim_name),
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Debiased objective",
      title = "Sparse MR-rr objective against eta"
    ) +
    ggplot2::theme_bw()

  sparsity_plot <- ggplot2::ggplot(
    eta_summary,
    ggplot2::aes(x = eta, y = zero_prop)
  ) +
    ggplot2::geom_step(direction = "hv") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(sim_name),
      ncol = 2
    ) +
    ggplot2::labs(
      x = expression(eta),
      y = "Zero proportion in B"
    ) +
    ggplot2::theme_bw()

  patchwork::wrap_plots(
    objective_plot,
    sparsity_plot,
    ncol = 1
  )
}


make_entrywise_data <- function(
    main_result,
    sparse_result,
    scenario) {
  component_labels <- c(
    C_ivw_list = "IVW",
    C_adivw_list = "SRIVW",
    AB_list = "Naive MR-rr",
    AB_d_list = "MR-rr",
    AB_d_r_list = "Reg. MR-rr",
    MrDAG_list = "MrDAG"
  )
  estimator_levels <- c(
    "IVW",
    "SRIVW",
    "Naive MR-rr",
    "MR-rr",
    "Reg. MR-rr",
    "Sparse MR-rr",
    "MrDAG"
  )
  position_grid <- expand.grid(
    Row = 0:2,
    Col = 0:8
  )
  rows <- list()
  row_index <- 1L

  for (component in names(component_labels)) {
    result_matrix <- main_result[[component]][[scenario]]
    validate_matrix(
      result_matrix,
      expected_rows = 27L,
      expected_columns = 1000L,
      label = paste(scenario, component)
    )
    finite_values <- !is.na(result_matrix)
    row_ids <- rep(position_grid$Row, each = ncol(result_matrix))
    column_ids <- rep(
      position_grid$Col,
      each = ncol(result_matrix)
    )
    values <- as.vector(t(result_matrix))

    rows[[row_index]] <- data.frame(
      Estimator = unname(component_labels[[component]]),
      Row = row_ids[as.vector(t(finite_values))],
      Col = column_ids[as.vector(t(finite_values))],
      Value = values[as.vector(t(finite_values))],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    row_index <- row_index + 1L
  }

  sparse_matrix <- sparse_result$C_sparse_list[[scenario]]
  validate_matrix(
    sparse_matrix,
    expected_rows = 27L,
    expected_columns = 1000L,
    label = paste(scenario, "C_sparse_list")
  )
  finite_values <- !is.na(sparse_matrix)
  row_ids <- rep(position_grid$Row, each = ncol(sparse_matrix))
  column_ids <- rep(
    position_grid$Col,
    each = ncol(sparse_matrix)
  )
  values <- as.vector(t(sparse_matrix))
  rows[[row_index]] <- data.frame(
    Estimator = "Sparse MR-rr",
    Row = row_ids[as.vector(t(finite_values))],
    Col = column_ids[as.vector(t(finite_values))],
    Value = values[as.vector(t(finite_values))],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  output <- do.call(rbind, rows)
  output$Estimator <- factor(
    output$Estimator,
    levels = estimator_levels
  )
  rownames(output) <- NULL
  output
}


make_entrywise_figure <- function(
    plot_data,
    true_effect,
    winsor_probability = c(0.99, 0.99, 0.99)) {
  if (length(winsor_probability) != 3L ||
      any(winsor_probability <= 0.5) ||
      any(winsor_probability > 1)) {
    stop(
      "The entrywise winsor probabilities are invalid.",
      call. = FALSE
    )
  }

  outcome_labels <- setNames(
    paste0("Outcome ", 1:3),
    as.character(0:2)
  )
  exposure_labels <- setNames(
    paste0("Exposure ", 1:9),
    as.character(0:8)
  )
  true_data <- expand.grid(
    Row = 0:2,
    Col = 0:8
  )
  true_data$TrueValue <- as.vector(true_effect)
  row_limits <- lapply(
    0:2,
    function(row) {
      row_values <- plot_data$Value[plot_data$Row == row]
      c(
        lower = unname(stats::quantile(
          row_values,
          1 - winsor_probability[[row + 1L]],
          na.rm = TRUE
        )),
        upper = unname(stats::quantile(
          row_values,
          winsor_probability[[row + 1L]],
          na.rm = TRUE
        ))
      )
    }
  )
  row_plots <- vector("list", 3L)

  for (row in 0:2) {
    row_data <- plot_data[
      plot_data$Row == row,
      ,
      drop = FALSE
    ]
    row_data$Outcome <- outcome_labels[[as.character(row)]]
    row_truth <- true_data[
      true_data$Row == row,
      ,
      drop = FALSE
    ]
    row_truth$Outcome <- outcome_labels[[as.character(row)]]
    show_x_text <- row == 2L
    show_top_strip <- row == 0L

    row_plots[[row + 1L]] <- ggplot2::ggplot(
      row_data,
      ggplot2::aes(
        x = Estimator,
        y = Value,
        fill = Estimator
      )
    ) +
      ggplot2::geom_boxplot(outlier.size = 0.3) +
      ggplot2::geom_hline(
        data = row_truth,
        ggplot2::aes(yintercept = TrueValue),
        color = "red",
        linetype = "dashed",
        linewidth = 0.5,
        inherit.aes = FALSE
      ) +
      ggplot2::coord_cartesian(ylim = row_limits[[row + 1L]]) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(Outcome),
        cols = ggplot2::vars(Col),
        labeller = ggplot2::labeller(Col = exposure_labels)
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.x = if (show_x_text) {
          ggplot2::element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
          )
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (show_x_text) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },
        strip.text.x = if (show_top_strip) {
          ggplot2::element_text(size = 9)
        } else {
          ggplot2::element_blank()
        },
        strip.background.x = if (show_top_strip) {
          ggplot2::element_rect()
        } else {
          ggplot2::element_blank()
        },
        strip.text.y = ggplot2::element_text(size = 9),
        legend.position = "right"
      ) +
      ggplot2::labs(x = "", y = "")
  }

  combined <- patchwork::wrap_plots(
    row_plots,
    ncol = 1,
    guides = "collect"
  )

  cowplot::ggdraw() +
    cowplot::draw_plot(
      combined,
      x = 0.025,
      y = 0,
      width = 0.975,
      height = 1
    ) +
    cowplot::draw_label(
      "Estimated C Value",
      x = 0.012,
      y = 0.5,
      angle = 90,
      size = 12
    )
}


canonicalize_loading <- function(loading) {
  for (pathway in seq_len(nrow(loading))) {
    largest_index <- which.max(abs(loading[pathway, ]))

    if (loading[pathway, largest_index] < 0) {
      loading[pathway, ] <- -loading[pathway, ]
    }
  }

  loading
}


make_loading_figure <- function(support_result) {
  true_loading <- canonicalize_loading(
    support_result$B_sparse
  )
  estimated_loading <- support_result$result_B_sparse
  rank <- nrow(true_loading)
  n_exposures <- ncol(true_loading)
  simulation_count <- ncol(estimated_loading)

  validate_matrix(
    estimated_loading,
    expected_rows = rank * n_exposures,
    expected_columns = 1000L,
    label = "Sparse loading estimates",
    allow_na = FALSE
  )

  key_index <- apply(abs(true_loading), 1L, which.max)
  key_value <- apply(abs(true_loading), 1L, max)
  pathway_order <- order(key_index, -key_value)
  true_loading <- true_loading[
    pathway_order,
    ,
    drop = FALSE
  ]
  aligned_loading <- matrix(
    NA_real_,
    nrow = rank * n_exposures,
    ncol = simulation_count
  )

  for (simulation in seq_len(simulation_count)) {
    loading <- matrix(
      estimated_loading[, simulation],
      nrow = rank,
      ncol = n_exposures
    )
    loading <- canonicalize_loading(loading)
    loading <- loading[pathway_order, , drop = FALSE]
    aligned_loading[, simulation] <- as.vector(loading)
  }

  entry_grid <- expand.grid(
    Row = 0:(rank - 1L),
    Col = 0:(n_exposures - 1L)
  )
  entry_grid <- entry_grid[
    order(entry_grid$Col, entry_grid$Row),
    ,
    drop = FALSE
  ]
  plot_data <- data.frame(
    Row = rep(entry_grid$Row, each = simulation_count),
    Col = rep(entry_grid$Col, each = simulation_count),
    Value = as.vector(t(aligned_loading)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  true_data <- data.frame(
    Row = entry_grid$Row,
    Col = entry_grid$Col,
    TrueValue = as.vector(true_loading),
    Is_zero = as.vector(true_loading == 0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  pathway_labels <- setNames(
    paste0("Pathway ", seq_len(rank)),
    as.character(0:(rank - 1L))
  )
  exposure_labels <- setNames(
    paste0("Exposure ", seq_len(n_exposures)),
    as.character(0:(n_exposures - 1L))
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = "", y = Value)
  ) +
    ggplot2::geom_boxplot(
      fill = "lightgreen",
      outlier.size = 0.3
    ) +
    ggplot2::geom_hline(
      data = true_data,
      ggplot2::aes(
        yintercept = TrueValue,
        color = Is_zero
      ),
      linetype = "dashed",
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "gray60", `FALSE` = "red")
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Row),
      cols = ggplot2::vars(Col),
      labeller = ggplot2::labeller(
        Row = pathway_labels,
        Col = exposure_labels
      ),
      scales = "free"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(size = 8),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = "",
      y = "Estimated B Entry",
      title = paste(
        "Simulation Boxplots of Estimated B Entries",
        "(Red = Nonzero True B, Gray = Zero True B)",
        sep = "\n"
      )
    )
}


make_prediction_data <- function(
    main_result,
    sparse_result,
    scenario) {
  component_labels <- c(
    Y_pred_C_ivw_list = "IVW",
    Y_pred_C_adivw_list = "SRIVW",
    Y_pred_AB_list = "Naive MR-rr",
    Y_pred_AB_d_list = "MR-rr",
    Y_pred_AB_d_r_list = "Reg. MR-rr",
    Y_pred_MrDAG_list = "MrDAG"
  )
  estimator_levels <- c(
    "IVW",
    "SRIVW",
    "Naive MR-rr",
    "MR-rr",
    "Reg. MR-rr",
    "Sparse MR-rr",
    "MrDAG"
  )
  rows <- list()
  row_index <- 1L

  for (component in names(component_labels)) {
    result_matrix <- main_result[[component]][[scenario]]
    validate_matrix(
      result_matrix,
      expected_rows = 3L,
      expected_columns = 1000L,
      label = paste(scenario, component)
    )
    finite_values <- !is.na(result_matrix)
    outcome_ids <- rep(0:2, each = ncol(result_matrix))
    values <- as.vector(t(result_matrix))

    rows[[row_index]] <- data.frame(
      Estimator = unname(component_labels[[component]]),
      Outcome = outcome_ids[as.vector(t(finite_values))],
      Value = values[as.vector(t(finite_values))],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    row_index <- row_index + 1L
  }

  sparse_matrix <-
    sparse_result$Y_pred_C_sparse_list[[scenario]]
  validate_matrix(
    sparse_matrix,
    expected_rows = 3L,
    expected_columns = 1000L,
    label = paste(scenario, "Y_pred_C_sparse_list")
  )
  finite_values <- !is.na(sparse_matrix)
  outcome_ids <- rep(0:2, each = ncol(sparse_matrix))
  values <- as.vector(t(sparse_matrix))
  rows[[row_index]] <- data.frame(
    Estimator = "Sparse MR-rr",
    Outcome = outcome_ids[as.vector(t(finite_values))],
    Value = values[as.vector(t(finite_values))],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  output <- do.call(rbind, rows)
  output$Estimator <- factor(
    output$Estimator,
    levels = estimator_levels
  )
  rownames(output) <- NULL
  output
}


make_prediction_figure <- function(plot_data, true_risk_score) {
  outcome_labels <- setNames(
    paste0("Outcome ", 1:3),
    as.character(0:2)
  )
  true_data <- data.frame(
    Outcome = 0:2,
    TrueValue = as.vector(true_risk_score),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Estimator)
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(y = Value, fill = Estimator),
      outlier.size = 0.3
    ) +
    ggplot2::geom_hline(
      data = true_data,
      ggplot2::aes(yintercept = TrueValue),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(Outcome),
      ncol = 3,
      labeller = ggplot2::labeller(Outcome = outcome_labels)
    ) +
    ggplot2::coord_cartesian(ylim = c(-2, 2)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(size = 10),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = "",
      y = "Estimated Risk Score"
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
  "patchwork",
  "cowplot",
  "scales"
))

scripts_root <- file.path(repo_root, "paper", "scripts")
validation_script <- file.path(
  scripts_root,
  "08_smoke_test_simulation_tables.R"
)
freeze_results_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project",
  "results"
)
eta_result_path <- file.path(
  freeze_results_root,
  "eta_selection_result_260729_regular_C.RData"
)

require_files(c(
  validation_script,
  eta_result_path
))

cat("Loading and validating archived simulation results.\n")
verified <- new.env(parent = globalenv())
sys.source(
  validation_script,
  envir = verified
)
eta_result <- load_required_object(
  eta_result_path,
  "eta_selection_result"
)

set.seed(123)
fixed_exposure <- stats::rnorm(9L)
fixed_exposure_target <- c(
  -0.560,
  -0.230,
   1.559,
   0.071,
   0.129,
   1.715,
   0.461,
  -1.265,
  -0.687
)

if (!identical(
  round(fixed_exposure, 3L),
  fixed_exposure_target
)) {
  stop(
    "The fixed prediction exposure profile differs from Table S2.",
    call. = FALSE
  )
}

generic_prediction_difference <-
  validate_prediction_reconstruction(
    verified$generic_main,
    verified$generic_sparse,
    fixed_exposure,
    "Generic-design"
  )
sparse_prediction_difference <-
  validate_prediction_reconstruction(
    verified$sparse_main,
    verified$sparse_sparse,
    fixed_exposure,
    "Sparse-loading-design"
  )

cat(
  "Archived prediction matrix reconstruction: PASS",
  " (maximum difference = ",
  format(
    max(
      generic_prediction_difference,
      sparse_prediction_difference
    ),
    scientific = TRUE
  ),
  ")\n",
  sep = ""
)

cat("Building supplementary Figure S1.\n")
figure_S01 <- make_eta_figure(eta_result$summary)

cat("Building supplementary Figures S3-S6.\n")
generic_strong_data <- make_entrywise_data(
  verified$generic_main,
  verified$generic_sparse,
  "me_1_effect_1"
)
generic_weak_data <- make_entrywise_data(
  verified$generic_main,
  verified$generic_sparse,
  "me_2.5_effect_0.25"
)
sparse_strong_data <- make_entrywise_data(
  verified$sparse_main,
  verified$sparse_sparse,
  "me_1_effect_1"
)
sparse_weak_data <- make_entrywise_data(
  verified$sparse_main,
  verified$sparse_sparse,
  "me_2.5_effect_0.25"
)

figure_S03 <- make_entrywise_figure(
  generic_strong_data,
  verified$generic_true_effect
)
figure_S04 <- make_entrywise_figure(
  generic_weak_data,
  verified$generic_true_effect,
  winsor_probability = c(0.98, 0.9925, 0.985)
)
figure_S05 <- make_entrywise_figure(
  sparse_strong_data,
  verified$sparse_true_effect
)
figure_S06 <- make_entrywise_figure(
  sparse_weak_data,
  verified$sparse_true_effect,
  winsor_probability = c(0.98, 0.9925, 0.985)
)

cat("Building supplementary Figure S7.\n")
figure_S07 <- make_loading_figure(
  verified$support_result
)

cat("Building supplementary Figures S8-S11.\n")
generic_true_risk <-
  verified$generic_true_effect %*% fixed_exposure
sparse_true_risk <-
  verified$sparse_true_effect %*% fixed_exposure

figure_S08 <- make_prediction_figure(
  make_prediction_data(
    verified$generic_main,
    verified$generic_sparse,
    "me_1_effect_1"
  ),
  generic_true_risk
)
figure_S09 <- make_prediction_figure(
  make_prediction_data(
    verified$generic_main,
    verified$generic_sparse,
    "me_2.5_effect_0.25"
  ),
  generic_true_risk
)
figure_S10 <- make_prediction_figure(
  make_prediction_data(
    verified$sparse_main,
    verified$sparse_sparse,
    "me_1_effect_1"
  ),
  sparse_true_risk
)
figure_S11 <- make_prediction_figure(
  make_prediction_data(
    verified$sparse_main,
    verified$sparse_sparse,
    "me_2.5_effect_0.25"
  ),
  sparse_true_risk
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

figure_labels <- c(
  "S1",
  "S3",
  "S4",
  "S5",
  "S6",
  "S7",
  "S8",
  "S9",
  "S10",
  "S11"
)
figure_paths <- c(
  write_png(
    figure_S01,
    file.path(output_root, "Figure_S01_simulation_eta_paths.png"),
    width = 10,
    height = 10
  ),
  write_png(
    figure_S03,
    file.path(
      output_root,
      "Figure_S03_generic_strong_entrywise.png"
    ),
    width = 10,
    height = 5
  ),
  write_png(
    figure_S04,
    file.path(
      output_root,
      "Figure_S04_generic_weak_entrywise.png"
    ),
    width = 10,
    height = 5
  ),
  write_png(
    figure_S05,
    file.path(
      output_root,
      "Figure_S05_sparse_strong_entrywise.png"
    ),
    width = 10,
    height = 5
  ),
  write_png(
    figure_S06,
    file.path(
      output_root,
      "Figure_S06_sparse_weak_entrywise.png"
    ),
    width = 10,
    height = 5
  ),
  write_png(
    figure_S07,
    file.path(output_root, "Figure_S07_sparse_loading_B.png"),
    width = 12,
    height = 4
  ),
  write_png(
    figure_S08,
    file.path(
      output_root,
      "Figure_S08_generic_strong_prediction.png"
    ),
    width = 12,
    height = 6
  ),
  write_png(
    figure_S09,
    file.path(
      output_root,
      "Figure_S09_generic_weak_prediction.png"
    ),
    width = 12,
    height = 6
  ),
  write_png(
    figure_S10,
    file.path(
      output_root,
      "Figure_S10_sparse_strong_prediction.png"
    ),
    width = 12,
    height = 6
  ),
  write_png(
    figure_S11,
    file.path(
      output_root,
      "Figure_S11_sparse_weak_prediction.png"
    ),
    width = 12,
    height = 6
  )
)

eta_summary_path <- file.path(
  output_root,
  "Figure_S01_simulation_eta_paths.csv"
)
utils::write.csv(
  as.data.frame(eta_result$summary),
  eta_summary_path,
  row.names = FALSE
)
cat("Wrote:", eta_summary_path, "\n")

artifact_paths <- c(
  figure_paths,
  normalizePath(
    eta_summary_path,
    winslash = "/",
    mustWork = TRUE
  )
)
artifact_labels <- c(
  paste("Figure", figure_labels),
  "Figure S1 source data"
)
manifest <- data.frame(
  Artifact = artifact_labels,
  File = basename(artifact_paths),
  Bytes = as.numeric(file.info(artifact_paths)$size),
  MD5 = unname(tools::md5sum(artifact_paths)),
  stringsAsFactors = FALSE,
  row.names = NULL
)
manifest_path <- file.path(
  output_root,
  "simulation_figure_manifest.csv"
)
utils::write.csv(
  manifest,
  manifest_path,
  row.names = FALSE
)
cat("Wrote:", manifest_path, "\n")

if (length(figure_paths) != 10L ||
    any(manifest$Bytes <= 0) ||
    any(!nzchar(manifest$MD5))) {
  stop(
    "The generated simulation figure inventory is incomplete.",
    call. = FALSE
  )
}

cat("\nGenerated simulation figures:", length(figure_paths), "\n")
cat("Supplementary Figures S1 and S3-S11: PASS\n")
