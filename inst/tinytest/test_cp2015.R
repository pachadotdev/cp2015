# solution without deficits ----
local({
  expected_results <- data.table::data.table(
        region = c("Canada", "Mexico", "USA"),
        tot = c(-0.11, -0.41, 0.04),
        vot = c(0.04, 1.72, 0.04),
        tech = c(0, 0, 0),
        welfare = c(-0.06, 1.31, 0.08),
        realwage = c(0.32, 1.72, 0.11)
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

# solution with deficits ----

local({
  expected_results <- data.table::data.table(
        region = c("Canada", "Mexico", "USA"),
        tot = c(-0.08, -0.41, 0.05),
        vot = c(0.04, 1.59, 0.04),
        tech = c(0, 0, 0),
        welfare = c(-0.04, 1.17, 0.08),
        realwage = c(0.33, 1.64, 0.12)
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
