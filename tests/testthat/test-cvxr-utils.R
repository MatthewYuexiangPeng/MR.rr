test_that("CVXR compatibility layer solves an OSQP matrix problem", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("osqp")
  skip_if(
    !"OSQP" %in% CVXR::installed_solvers(),
    "OSQP is not available"
  )

  target <- matrix(c(1, -2), nrow = 1)

  beta <- .cvxr_matrix_variable(
    n_rows = 1,
    n_cols = 2
  )

  objective <- CVXR::sum_squares(beta - target)

  problem <- CVXR::Problem(
    CVXR::Minimize(objective)
  )

  observed <- .solve_cvxr_variable(
    problem = problem,
    variable = beta,
    solver = "OSQP"
  )

  expect_equal(dim(observed), c(1L, 2L))
  expect_equal(observed, target, tolerance = 1e-5)
})
