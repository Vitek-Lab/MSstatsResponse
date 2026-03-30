test_that("build_concentration_ladders returns correct subset sizes", {
  concs <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
  result <- MSstatsResponse:::.build_concentration_ladders(concs, c(2, 5))

  expect_equal(length(result), 4)
  expect_equal(names(result), c("2", "3", "4", "5"))

  # All subsets should contain control and highest dose
  for (k in names(result)) {
    expect_true(0 %in% result[[k]])
    expect_true(3000 %in% result[[k]])
    expect_equal(length(result[[k]]), as.integer(k))
  }
})

test_that("build_concentration_ladders always includes control and max", {
  concs <- c(0, 5, 50, 200, 1000)
  result <- MSstatsResponse:::.build_concentration_ladders(concs, c(2, 5))

  for (k in names(result)) {
    expect_true(0 %in% result[[k]])
    expect_true(1000 %in% result[[k]])
  }
})

test_that("build_concentration_ladders errors on fewer than 2 concentrations", {
  expect_error(
    MSstatsResponse:::.build_concentration_ladders(c(0), c(2, 3)),
    "At least 2 unique"
  )
})

test_that("build_concentration_ladders handles arbitrary user concentrations", {
  concs <- c(0, 2, 4, 20, 50, 100, 500, 1000)
  result <- MSstatsResponse:::.build_concentration_ladders(concs, c(2, 6))

  expect_equal(length(result), 5)
  for (k in names(result)) {
    subset <- result[[k]]
    expect_true(all(subset %in% concs))
  }
})

test_that("run_tpr_simulation validates rep_range", {
  expect_error(
    run_tpr_simulation(c(3, 1), c(0, 100, 1000), c(2, 3)),
    "min <= max"
  )
  expect_error(
    run_tpr_simulation("bad", c(0, 100, 1000), c(2, 3)),
    "numeric vector"
  )
})

test_that("run_tpr_simulation validates dose_range", {
  expect_error(
    run_tpr_simulation(c(1, 2), c(0, 100, 1000), c(5, 2)),
    "min <= max"
  )
  expect_error(
    run_tpr_simulation(c(1, 2), c(0, 100, 1000), c(1, 3)),
    "at least 2"
  )
})

test_that("run_tpr_simulation returns correct structure", {
  results <- run_tpr_simulation(
    rep_range = c(1, 2),
    concentrations = c(0, 10, 100, 1000),
    dose_range = c(2, 3),
    n_proteins = 50
  )

  expect_true(is.data.frame(results))
  expect_true(all(c("Interaction", "TPR", "N_rep", "NumConcs") %in% colnames(results)))
  expect_true(all(results$Interaction %in% c("Strong", "Weak")))
  expect_true(all(results$TPR >= 0 & results$TPR <= 100))
})

test_that("plot_tpr_power_curve returns a plotly object", {
  results <- run_tpr_simulation(
    rep_range = c(1, 2),
    concentrations = c(0, 10, 100, 1000),
    dose_range = c(2, 3),
    n_proteins = 50
  )

  p <- plot_tpr_power_curve(results)
  expect_true(inherits(p, "plotly"))
})

test_that("run_tpr_simulation accepts data and protein parameters", {
  # Verify the function signature accepts these params without error
  # (actual template extraction is tested via integration/manual testing)
  expect_error(
    run_tpr_simulation(
      rep_range = c(1, 1),
      concentrations = c(0, 100, 1000),
      dose_range = c(2, 3),
      n_proteins = 10,
      data = NULL,
      protein = NULL
    ),
    NA  # NA means "expect NO error"
  )
})