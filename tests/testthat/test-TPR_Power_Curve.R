test_that("build_concentration_ladders returns correct subset sizes", {
  concs <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
  result <- MSstatsResponse:::.build_concentration_ladders(concs, c(2, 5))

  expect_equal(length(result), 4)
  expect_equal(names(result), c("2", "3", "4", "5"))

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

test_that("build_concentration_ladders errors without control dose", {
  expect_error(
    MSstatsResponse:::.build_concentration_ladders(c(1, 10, 100), c(2, 3)),
    "must include 0"
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
  mock_data <- data.frame(protein = "P1", drug = "D1", dose = 0, response = 20)
  expect_error(
    run_tpr_simulation(c(3, 1), c(0, 100, 1000), c(2, 3), data = mock_data, protein = "P1"),
    "min <= max"
  )
  expect_error(
    run_tpr_simulation("bad", c(0, 100, 1000), c(2, 3), data = mock_data, protein = "P1"),
    "numeric vector"
  )
})

test_that("run_tpr_simulation validates dose_range", {
  mock_data <- data.frame(protein = "P1", drug = "D1", dose = 0, response = 20)
  expect_error(
    run_tpr_simulation(c(1, 2), c(0, 100, 1000), c(5, 2), data = mock_data, protein = "P1"),
    "min <= max"
  )
  expect_error(
    run_tpr_simulation(c(1, 2), c(0, 100, 1000), c(1, 3), data = mock_data, protein = "P1"),
    "at least 2"
  )
})

test_that("run_tpr_simulation returns correct structure", {
  mock_data <- data.frame(
    protein = rep("P1", 6),
    drug = c("DMSO", "DMSO", "Drug1", "Drug1", "Drug1", "Drug1"),
    dose = c(0, 0, 10, 100, 1000, 1000),
    response = c(20, 20, 18, 15, 12, 12)
  )

  results <- run_tpr_simulation(
    rep_range = c(1, 2),
    concentrations = c(0, 10, 100, 1000),
    dose_range = c(2, 3),
    data = mock_data,
    protein = "P1",
    n_proteins = 50
  )

  expect_true(is.data.frame(results))
  expect_true(all(c("Interaction", "TPR", "N_rep", "NumConcs") %in% colnames(results)))
  expect_true(all(results$Interaction == "Strong"))
  expect_true(all(results$TPR >= 0 & results$TPR <= 100))

  # TPR should generally increase with more replicates for a given concentration count
  for (nc in unique(results$NumConcs)) {
    subset <- results[results$NumConcs == nc, ]
    subset <- subset[order(subset$N_rep), ]
    if (nrow(subset) >= 2) {
      # Allow for some noise, but overall trend should be non-decreasing
      expect_true(subset$TPR[nrow(subset)] >= subset$TPR[1],
                  info = paste("TPR should increase with replicates at NumConcs =", nc))
    }
  }

  # TPR should generally increase with more concentrations for a given replicate count
  for (nr in unique(results$N_rep)) {
    subset <- results[results$N_rep == nr, ]
    subset <- subset[order(subset$NumConcs), ]
    if (nrow(subset) >= 2) {
      expect_true(subset$TPR[nrow(subset)] >= subset$TPR[1],
                  info = paste("TPR should increase with concentrations at N_rep =", nr))
    }
  }
})

test_that("plot_tpr_power_curve returns a plotly object", {
  mock_data <- data.frame(
    protein = rep("P1", 6),
    drug = c("DMSO", "DMSO", "Drug1", "Drug1", "Drug1", "Drug1"),
    dose = c(0, 0, 10, 100, 1000, 1000),
    response = c(20, 20, 18, 15, 12, 12)
  )

  results <- run_tpr_simulation(
    rep_range = c(1, 2),
    concentrations = c(0, 10, 100, 1000),
    dose_range = c(2, 3),
    data = mock_data,
    protein = "P1",
    n_proteins = 50
  )

  p <- plot_tpr_power_curve(results)
  expect_true(inherits(p, "plotly"))
})