# solution without deficits ----
local({
  expected_results <- data.table::data.table(
        region = c("Canada", "Mexico", "USA"),
        total = c(-0.06, 1.31, 0.08),
        tot = c(-0.11, -0.41, 0.04),
        vot = c(0.04, 1.72, 0.04),
        realwage = c(0.32, 1.72, 0.11),
        tech = c(0, 0, 0)
    )

    results <- run_cp2015(
        data = cp2015nafta,
        zero_aggregate_deficit = TRUE,
        tol = 1e-7,
        verbose = TRUE
    )

    nafta <- c("Canada", "Mexico", "USA")
    
    expect_true(all(vapply(results, data.table::is.data.table, logical(1))))
    
    results <- results$welfare[region %in% nafta]
    num_cols <- names(results)[vapply(results, is.numeric, logical(1))]
    results[, (num_cols) := lapply(.SD, round, 2), .SDcols = num_cols]
    
    expect_equal(results, expected_results)
})

# typed solver boundary ----
local({
  prepare <- getFromNamespace("prepare_cp2015", "cp2015")
  slices <- getFromNamespace("array_to_matrices", "cp2015")
  data <- prepare(cp2015nafta, zero_aggregate_deficit = TRUE)
  before <- serialize(data, NULL)
  solution <- suppressWarnings(cp2015:::solve_model(data, maxiter = 2L, trace = FALSE))

  expect_identical(serialize(data, NULL), before)
  expect_identical(
    suppressWarnings(cp2015:::solve_model(data, maxiter = 2L, trace = FALSE)),
    solution
  )
  for (name in names(data$variables)) {
    expect_identical(attributes(solution$variables[[name]]),
                     attributes(data$variables[[name]]))
  }
  expect_identical(dimnames(solution$variables$P_n_hat),
                   dimnames(data$variables$w_n_hat))
  expect_true(any(solution$variables$X_nj1 != data$variables$X_nj1))

  x <- array(as.numeric(seq_len(6L)), c(1L, 2L, 3L))
  matrices <- slices(x)
  expect_identical(lapply(matrices, dim), rep(list(c(1L, 2L)), 3L))
  expect_identical(unlist(matrices, use.names = FALSE), as.numeric(x))
  expect_error(slices(matrix(1, 2, 2)), "three-dimensional")

  invalid <- data
  invalid$parameters$theta_j <- 1
  expect_error(cp2015:::solve_model(invalid, trace = FALSE), "theta_j.*length")
  invalid <- data
  invalid$parameters$pi_nij0 <- array(1, c(1L, 1L, 40L))
  expect_error(cp2015:::solve_model(invalid, trace = FALSE), "pi_nij0.*dimensions")
  expect_error(cp2015:::solve_model(data, maxiter = 0L), "positive integer")
  expect_error(cp2015:::solve_model(data, triter = 0L), "positive integer")

  data$sets <- list(regions = "region", sectors = "sector")
  for (name in names(data$parameters)) {
    x <- data$parameters[[name]]
    data$parameters[[name]] <- array(
      1, dim = rep(1L, length(dim(x))),
      dimnames = lapply(dimnames(x), head, 1L)
    )
  }
  for (name in names(data$variables)) {
    x <- data$variables[[name]]
    data$variables[[name]] <- array(
      1, dim = rep(1L, length(dim(x))),
      dimnames = lapply(dimnames(x), head, 1L)
    )
  }
  data$parameters$gamma_nj[] <- 0.5
  data$parameters$gamma_nkj[] <- 0.5
  data$parameters$tau_nij0[] <- 0
  data$variables$tau_nij1[] <- 0
  data$variables$D_n[] <- 0
  data$variables$X_nj1[] <- 2
  solution <- cp2015:::solve_model(data, trace = FALSE)
  expect_identical(solution$message, "Successful convergence")
  expect_equal(solution$convergence_criteria, 0)
  expect_identical(dim(solution$variables$pi_nij1), c(1L, 1L, 1L))
  expect_equal(as.numeric(solution$variables$P_n_hat), 1)
})

# solution with deficits ----

local({
  expected_results <- data.table::data.table(
        region = c("Canada", "Mexico", "USA"),
        total = c(-0.04, 1.17, 0.08),
        tot = c(-0.08, -0.41, 0.05),
        vot = c(0.04, 1.59, 0.04),
        realwage = c(0.33, 1.64, 0.12),
        tech = c(0, 0, 0)
        )
    
    results <- run_cp2015(
        data = cp2015nafta,
        zero_aggregate_deficit = FALSE,
        tol = 1e-7,
        verbose = TRUE
    )

    nafta <- c("Canada", "Mexico", "USA")

    expect_true(all(vapply(results, data.table::is.data.table, logical(1))))

    results <- results$welfare[region %in% nafta]
    num_cols <- names(results)[vapply(results, is.numeric, logical(1))]
    results[, (num_cols) := lapply(.SD, round, 2), .SDcols = num_cols]

    expect_equal(results, expected_results)
})
