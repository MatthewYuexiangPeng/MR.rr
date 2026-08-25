.cvxr_matrix_variable <- function(n_rows, n_cols) {
  if (utils::packageVersion("CVXR") >= "1.8.0") {
    CVXR::Variable(c(n_rows, n_cols))
  } else {
    CVXR::Variable(rows = n_rows, cols = n_cols)
  }
}

.solve_cvxr_variable <- function(problem, variable, solver = "OSQP") {
  if (!solver %in% CVXR::installed_solvers()) {
    stop(
      sprintf("CVXR solver '%s' is not installed.", solver),
      call. = FALSE
    )
  }

  if (utils::packageVersion("CVXR") >= "1.8.0") {
    CVXR::psolve(problem, solver = solver)
    solution <- CVXR::value(variable)
  } else {
    result <- CVXR::solve(problem, solver = solver)
    solution <- result$getValue(variable)
  }

  if (is.null(solution) || any(!is.finite(solution))) {
    stop("CVXR did not return a finite solution.", call. = FALSE)
  }

  as.matrix(solution)
}
