# test_predictIC50.R
# Unit tests for predictIC50.R functions

library(testthat)
library(MSstatsResponse)

# Helper to create multi-protein test data
create_multi_protein_data <- function() {
  data.frame(
    protein = rep(c("P1", "P2", "P3"), each = 12),
    drug = rep(c(rep("DMSO", 3), rep("Drug1", 9)), 3),
    dose = rep(c(0, 0, 0, rep(c(1e-9, 1e-8, 1e-7), each = 3)), 3),
    response = c(
      # P1: Strong responder
      20, 20.1, 19.9, 18, 17.9, 18.1, 15, 14.9, 15.1, 12, 11.9, 12.1,
      # P2: Moderate responder
      21, 20.9, 21.1, 20, 19.9, 20.1, 19, 18.9, 19.1, 18, 17.9, 18.1,
      # P3: Non-responder
      19, 19.1, 18.9, 19, 18.9, 19.1, 18.8, 19.0, 19.2, 19.1, 18.9, 19.0
    )
  )
}

# Tests for main predictIC50 function
test_that("predictIC50 returns correct structure", {
  test_data <- create_multi_protein_data()

  suppressWarnings({
  results <- predictIC50(
    data = test_data,
    n_samples = 50,  # Small for speed
    bootstrap = FALSE
  )
  })

  expect_s3_class(results, "data.frame")
  expect_named(results, c("protein", "drug", "IC50", "IC50_lower_bound", "IC50_upper_bound"))
  expect_equal(nrow(results), 3)  # 3 proteins
})

test_that("predictIC50 processes all protein-drug combinations", {
  test_data <- create_multi_protein_data()
  # Add another drug
  test_data2 <- test_data
  test_data2$drug[test_data2$drug == "Drug1"] <- "Drug2"
  combined_data <- rbind(test_data, test_data2)

  suppressWarnings({
  results <- predictIC50(
    data = combined_data,
    bootstrap = FALSE
  )
  })

  expect_equal(nrow(results), 6)  # 3 proteins × 2 drugs
  expect_equal(sort(unique(results$drug)), c("Drug1", "Drug2"))
})

test_that("predictIC50 handles bootstrap option", {
  test_data <- create_multi_protein_data()[1:12, ]  # Just P1 for speed

  # Without bootstrap
  results_no_boot <- predictIC50(
    data = test_data,
    bootstrap = FALSE
  )

  expect_true(all(is.na(results_no_boot$IC50_lower_bound)))
  expect_true(all(is.na(results_no_boot$IC50_upper_bound)))

  # With bootstrap
  results_boot <- predictIC50(
    data = test_data,
    n_samples = 50,
    bootstrap = TRUE
  )

  expect_false(all(is.na(results_boot$IC50_lower_bound)))
  expect_false(all(is.na(results_boot$IC50_upper_bound)))
})

test_that("predictIC50 respects transform_dose parameter", {
  test_data <- create_multi_protein_data()[1:12, ]

  results_transform <- predictIC50(
    data = test_data,
    transform_dose = TRUE,
    bootstrap = FALSE
  )

  results_no_transform <- predictIC50(
    data = test_data,
    transform_dose = FALSE,
    bootstrap = FALSE
  )

  # Results should differ
  expect_false(isTRUE(all.equal(results_transform$IC50, results_no_transform$IC50)))
})

test_that("predictIC50 handles ratio vs log scale", {
  test_data <- create_multi_protein_data()[1:12, ]

  results_ratio <- predictIC50(
    data = test_data,
    ratio_response = TRUE,
    bootstrap = FALSE
  )

  results_log <- predictIC50(
    data = test_data,
    ratio_response = FALSE,
    bootstrap = FALSE
  )

  # IC50 values should differ between scales
  expect_false(isTRUE(all.equal(results_ratio$IC50, results_log$IC50)))
})

test_that("predictIC50 handles different target responses", {
  test_data <- create_multi_protein_data()[1:12, ]

  ic25 <- predictIC50(
    data = test_data,
    target_response = 0.25,
    bootstrap = FALSE
  )

  ic50 <- predictIC50(
    data = test_data,
    target_response = 0.50,
    bootstrap = FALSE
  )

  ic75 <- predictIC50(
    data = test_data,
    target_response = 0.75,
    bootstrap = FALSE
  )

  # IC values should be ordered (for decreasing response)
  # Lower IC value means more potent, so IC75 < IC50 < IC25
  expect_true(ic75$IC50[1] > ic50$IC50[1])
  expect_true(ic50$IC50[1] > ic25$IC50[1])
})


test_that("predictIC50 uses parallel processing when requested", {
  test_data <- create_multi_protein_data()

  # Time single core
  time_single <- system.time({
    suppressWarnings({
      results_single <- predictIC50(
        data = test_data,
        n_samples = 50,
        bootstrap = TRUE,
        BPPARAM = BiocParallel::SerialParam()
      )
    })
  })

  # Time multi core (if available)
  if (parallel::detectCores() > 1) {
    time_multi <- system.time({
      suppressWarnings({
        results_multi <- predictIC50(
          data = test_data,
          n_samples = 50,
          bootstrap = TRUE,
          BPPARAM = BiocParallel::MulticoreParam(workers = 2)
        )
      })
    })

    # Results should be similar
    expect_equal(nrow(results_single), nrow(results_multi))
  }
})

# In test-predictIC50.R, update the test to suppress the expected warning:

test_that("predictIC50 handles proteins with no dose response", {
  # Flat response data
  flat_data <- data.frame(
    protein = rep("P_flat", 6),
    drug = c(rep("DMSO", 3), rep("Drug1", 3)),
    dose = c(0, 0, 0, 1e-7, 1e-6, 1e-5),
    response = rep(20, 6) + rnorm(6, 0, 0.1)
  )

  # Suppress the expected warning about target not reached
  suppressWarnings({
    results <- predictIC50(
      data = flat_data,
      bootstrap = FALSE
    )
  })

  # Should return NA or high IC50
  expect_true(is.na(results$IC50) || results$IC50 > 5)
})

test_that("predictIC50 filters out DMSO-only proteins", {
  dmso_only <- data.frame(
    protein = "P1",
    drug = "DMSO",
    dose = c(0, 0, 0),
    response = c(20, 20.1, 19.9)
  )

  results <- predictIC50(
    data = dmso_only,
    bootstrap = FALSE
  )

  expect_equal(nrow(results), 0)
})

test_that("predictIC50 handles missing data appropriately", {
  incomplete_data <- data.frame(
    protein = c("P1", "P1", "P1", "P2"),
    drug = c("DMSO", "Drug1", "Drug1", "Drug1"),
    dose = c(0, 1e-7, 1e-6, 1e-6),
    response = c(20, 15, 10, 18)
  )


    results <- predictIC50(
      data = incomplete_data,
      bootstrap = FALSE
    )

  # Should process P1 successfully
  expect_equal(nrow(results), 2)
  expect_false(is.na(results$IC50[results$protein == "P1"]))
  expect_true(is.na(results$IC50[results$protein == "P2"]))
})


test_that("predictIC50 returns pIC50 values", {
  test_data <- create_multi_protein_data()[1:12, ]

  results <- predictIC50(
    data = test_data,
    bootstrap = FALSE
  )

  # IC50 should be in pIC50 format (negative log10)
  # For reasonable IC50 in M (e.g., 1e-7 M), pIC50 should be around 7
  expect_true(results$IC50[1] > 0)  # pIC50 should be positive
  expect_true(results$IC50[1] < 20)  # But reasonable
})

test_that("predictIC50 handles edge case doses", {
  edge_data <- data.frame(
    protein = rep("P1", 6),
    drug = c(rep("DMSO", 3), rep("Drug1", 3)),
    dose = c(0, 0, 0, 1e-12, 1e-3, 1),  # Very wide range
    response = c(20, 20, 20, 18, 10, 5)
  )

  results <- predictIC50(
    data = edge_data,
    bootstrap = FALSE,
    transform_dose = TRUE
  )

  expect_equal(nrow(results), 1)
})


