rank_selection_test_data <- function() {
  set.seed(20260827)

  W <- matrix(
    c(
      1.4, 0.2, 0.1,
      0.2, 1.2, 0.15,
      0.1, 0.15, 1.0
    ),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )

  list(
    Y = matrix(stats::rnorm(180), nrow = 60, ncol = 3),
    X = matrix(stats::rnorm(240), nrow = 60, ncol = 4),
    Sigma_X = diag(rep(1e-4, 4)),
    W = W
  )
}


reference_rank_selection <- function(
    Y,
    X,
    Sigma_X,
    W,
    alpha,
    min_rank) {
  n <- nrow(Y)
  px <- ncol(X)
  py <- ncol(Y)
  full_rank <- min(px, py)

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)
  W_sqrt <- .sqrt_matrix(W)

  V <-
    W_sqrt %*%
    t(Sigma_xy) %*%
    debiased_Sigma_xx_inv %*%
    Sigma_xy %*%
    W_sqrt

  eigenvalues <- eigen(V)$values[seq_len(full_rank)]
  candidate_ranks <- seq.int(min_rank, full_rank - 1L)

  statistic <- vapply(candidate_ranks, function(rank) {
    tail_indices <- seq.int(rank + 1L, full_rank)
    (n - (px + py + 1) / 2) *
      sum(log(1 + eigenvalues[tail_indices]))
  }, numeric(1))

  df <- (py - candidate_ranks) * (px - candidate_ranks)
  p_value <- 1 - stats::pchisq(statistic, df = df)
  accepted <- which(p_value >= alpha)

  selected_rank <- if (length(accepted) == 0L) {
    full_rank
  } else {
    candidate_ranks[accepted[1L]]
  }

  list(
    selected_rank = selected_rank,
    candidate_ranks = candidate_ranks,
    statistic = statistic,
    df = df,
    p_value = p_value,
    eigenvalues = eigenvalues
  )
}


test_that("select_mr_rr_rank matches the paper diagnostic", {
  dat <- rank_selection_test_data()

  fit <- select_mr_rr_rank(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = dat$W,
    alpha = 0.05,
    min_rank = 1
  )

  reference <- reference_rank_selection(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = dat$W,
    alpha = 0.05,
    min_rank = 1
  )

  expect_s3_class(fit, "mr_rr_rank_selection")
  expect_identical(fit$selected_rank, as.integer(reference$selected_rank))
  expect_identical(
    fit$candidate_ranks,
    as.integer(reference$candidate_ranks)
  )
  expect_equal(fit$statistic, reference$statistic, tolerance = 1e-10)
  expect_equal(fit$df, reference$df, tolerance = 0)
  expect_equal(fit$p_value, reference$p_value, tolerance = 1e-10)
  expect_equal(fit$eigenvalues, reference$eigenvalues, tolerance = 1e-10)
})


test_that("select_mr_rr_rank handles the full-rank boundary", {
  dat <- rank_selection_test_data()

  fit <- select_mr_rr_rank(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = dat$W,
    min_rank = 3
  )

  expect_identical(fit$selected_rank, 3L)
  expect_length(fit$candidate_ranks, 0L)
  expect_length(fit$statistic, 0L)
  expect_length(fit$df, 0L)
  expect_length(fit$p_value, 0L)
})


test_that("identity W is used when W is omitted", {
  dat <- rank_selection_test_data()

  implicit <- select_mr_rr_rank(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X
  )

  explicit <- select_mr_rr_rank(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = diag(ncol(dat$Y))
  )

  expect_identical(implicit$selected_rank, explicit$selected_rank)
  expect_equal(implicit$statistic, explicit$statistic, tolerance = 1e-12)
  expect_equal(implicit$p_value, explicit$p_value, tolerance = 1e-12)
})


test_that("select_mr_rr_rank validates its controls", {
  dat <- rank_selection_test_data()

  common_args <- list(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = dat$W
  )

  expect_error(
    do.call(select_mr_rr_rank, c(common_args, list(alpha = 0))),
    "`alpha`"
  )
  expect_error(
    do.call(select_mr_rr_rank, c(common_args, list(alpha = 1))),
    "`alpha`"
  )
  expect_error(
    do.call(select_mr_rr_rank, c(common_args, list(min_rank = -1))),
    "`min_rank`"
  )
  expect_error(
    do.call(select_mr_rr_rank, c(common_args, list(min_rank = 4))),
    "`min_rank`"
  )
})
