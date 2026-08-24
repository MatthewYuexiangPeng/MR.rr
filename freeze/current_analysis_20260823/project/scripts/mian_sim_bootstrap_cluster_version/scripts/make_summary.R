# Combine all bootstrap CI results into summary values
suppressPackageStartupMessages(library(dplyr))

# sim setting 1
# sim_tag <- "regular_C"

# sim setting 2
sim_tag <- "sparseB_C"

files <- list.files(
  "results",
  pattern = paste0(
    "^CI_(no_mrdag|mrdag_only)_idx.*_",
    sim_tag,
    "\\.RData$"
  ),
  full.names = TRUE
)

cat(sprintf("Found %d result files\n", length(files)))

# ------------------------------------------------------------ #
# Helpers
# ------------------------------------------------------------ #

rename_estimator <- function(est) {
  dplyr::recode(
    est,
    "adIVW" = "SRIVW",
    "Naive" = "Naive MR-rr",
    "MR" = "MR-rr",
    "MR_r" = "Reg. MR-rr",
    .default = est
  )
}

est_order <- c(
  "IVW",
  "SRIVW",
  "Naive MR-rr",
  "MR-rr",
  "Reg. MR-rr",
  "MrDAG"
)

summarize_med_iqr <- function(x, scale = 1, digits = 3) {
  x <- scale * x
  
  c(
    med = round(median(x, na.rm = TRUE), digits),
    q1  = round(unname(quantile(x, 0.25, na.rm = TRUE)), digits),
    q3  = round(unname(quantile(x, 0.75, na.rm = TRUE)), digits)
  )
}

compute_avg_se_by_entry <- function(estimates_list) {
  # estimates_list: list over Monte Carlo simulations
  # each element: bootstrap_size x d matrix
  se_mat <- sapply(estimates_list, function(bt_mat) {
    apply(bt_mat, 2, sd, na.rm = TRUE)
  })
  
  if (is.vector(se_mat)) {
    se_mat <- matrix(se_mat, nrow = length(se_mat))
  }
  
  # average bootstrap SE across simulations, entrywise
  rowMeans(se_mat, na.rm = TRUE)
}

compute_cp_by_entry <- function(coverage_mat) {
  # coverage_mat: d x Monte Carlo simulations
  rowMeans(coverage_mat, na.rm = TRUE)
}

# ------------------------------------------------------------ #
# Main
# ------------------------------------------------------------ #

rows <- list()

for (f in files) {
  load(f)  # loads object `result`
  
  s <- result$setting
  setting_tag <- sprintf("me=%s, effect=%s", s$me_weight, s$effect_weight)
  
  for (est in names(result$coverage_matrix)) {
    
    est_label <- rename_estimator(est)
    
    # CP: entrywise coverage probability across simulations
    cp_entry <- compute_cp_by_entry(result$coverage_matrix[[est]])
    cp_sum <- summarize_med_iqr(cp_entry, scale = 100, digits = 1)
    
    # SE: average bootstrap SE across simulations, then summarize across entries
    if (!is.null(result$estimates) && !is.null(result$estimates[[est]])) {
      avg_se_entry <- compute_avg_se_by_entry(result$estimates[[est]])
      se_sum <- summarize_med_iqr(avg_se_entry, scale = 1, digits = 3)
    } else {
      se_sum <- c(med = NA_real_, q1 = NA_real_, q3 = NA_real_)
    }
    
    rows[[length(rows) + 1]] <- data.frame(
      Setting = setting_tag,
      Estimator = est_label,
      
      SE_med = se_sum["med"],
      SE_q1  = se_sum["q1"],
      SE_q3  = se_sum["q3"],
      
      CP_med = cp_sum["med"],
      CP_q1  = cp_sum["q1"],
      CP_q3  = cp_sum["q3"],
      
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
}

setting_order <- c(
  "me=2.5, effect=0.25",
  "me=1, effect=0.25",
  "me=2.5, effect=1",
  "me=1, effect=1"
)

summary_table <- do.call(rbind, rows) |>
  mutate(
    Setting = factor(Setting, levels = setting_order),
    Estimator = factor(Estimator, levels = est_order)
  ) |>
  arrange(Setting, Estimator)

print(summary_table)

out_csv <- paste0(
  "results/bootstrap_summary_SE_CP_med_IQR_",
  sim_tag,
  ".csv"
)

write.csv(
  summary_table,
  out_csv,
  row.names = FALSE
)

cat("Saved to", out_csv, "\n")


# # Combine all bootstrap CI results into one summary table
# suppressPackageStartupMessages(library(dplyr))
# 
# files <- list.files("results",
#                     pattern = "^CI_(no_mrdag|mrdag_only)_idx.*\\.RData$",
#                     full.names = TRUE)
# cat(sprintf("Found %d result files\n", length(files)))
# 
# rows <- list()
# for (f in files) {
#   load(f)  # `result`
#   s <- result$setting
#   setting_tag <- sprintf("me=%s, effect=%s", s$me_weight, s$effect_weight)
#   for (est in names(result$med_coverage)) {
#     rows[[length(rows) + 1]] <- data.frame(
#       Setting          = setting_tag,
#       Estimator        = est,
#       Coverage         = round(result$med_coverage[[est]], 4),
#       CI_Median_Length = round(result$ci_length_summary[[est]]$median, 3),
#       CI_Mean_Length   = round(result$ci_length_summary[[est]]$mean,   3),
#       stringsAsFactors = FALSE)
#   }
# }
# 
# summary_table <- do.call(rbind, rows) |>
#   arrange(Setting,
#           factor(Estimator, levels = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG")))
# print(summary_table)
# write.csv(summary_table, "results/summary_table_all.csv", row.names = FALSE)
# cat("Saved to results/summary_table_all.csv\n")