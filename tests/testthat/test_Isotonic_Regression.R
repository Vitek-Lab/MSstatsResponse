# test_Fit_Isotonic_Regression.R
# Unit tests for Fit_Isotonic_Regression.R functions

library(testthat)
library(MSstatsResponse)

# Helper function to create test data
create_test_protein_data <- function() {
  # Create data with clear dose-response for some proteins
  data.frame(
    protein = rep(c("P1", "P2", "P3"), each = 15),
    drug = rep(c(rep("DMSO", 3), rep("Drug1", 12)), 3),
    dose = rep(c(0, 0, 0, rep(c(1e-9, 1e-8, 1e-7, 1e-6), each = 3)), 3),
    response = c(
      # P1: Strong responder
      20, 20.1, 19.9, 19, 18.9, 19.1, 17, 17.2, 16.8, 14, 14.2, 13.8, 10, 10.1, 9.9,
      # P2: Weak responder
      21, 20.9, 21.1, 20.5, 20.4, 20.6, 19.8, 19.9, 19.7, 19, 19.1, 18.9, 18.5, 18.4, 18.6,
      # P3: Non-responder
      19, 19.2, 18.8, 19.1, 18.9, 19.0, 19.2, 18.8, 19.1, 18.9, 19.0, 19.1, 19.0, 18.9, 19.2
    )
  )
}

# Tests for doseResponseFit
test_that("doseResponseFit returns correct structure", {
  test_data <- create_test_protein_data()

  results <- doseResponseFit(
    data = test_data,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = FALSE
  )

  expect_s3_class(results, "data.frame")
  expect_true("Protein" %in% names(results))
  expect_true("drug" %in% names(results))
  expect_true("F_statistic" %in% names(results))
  expect_true("P_value" %in% names(results))
  expect_true("adj.pvalue" %in% names(results))
})

test_that("doseResponseFit processes all proteins", {
  test_data <- create_test_protein_data()

  results <- doseResponseFit(
    data = test_data,
    increasing = FALSE
  )

  expect_equal(nrow(results), 3)  # 3 proteins
  expect_equal(unique(results$Protein), c("P1", "P2", "P3"))
})

test_that("doseResponseFit handles multiple drugs", {
  test_data <- create_test_protein_data()
  # Add another drug
  test_data2 <- test_data
  test_data2$drug[test_data2$drug == "Drug1"] <- "Drug2"
  combined_data <- rbind(test_data, test_data2)

  results <- doseResponseFit(
    data = combined_data,
    increasing = FALSE
  )

  expect_equal(nrow(results), 6)  # 3 proteins × 2 drugs
  expect_equal(sort(unique(results$drug)), c("Drug1", "Drug2"))
})


test_that("doseResponseFit handles custom weights", {
  test_data <- create_test_protein_data()

  # Create weights for the full dataset
  weights <- rep(1, nrow(test_data))
  weights[1:5] <- 2  # Give more weight to first observations

  results <- doseResponseFit(
    data = test_data,
    weights = weights,
    increasing = FALSE
  )

  # Should return results for all proteins
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 3)
  expect_equal(length(unique(results$Protein)), 3)
})

test_that("doseResponseFit filters out DMSO from drug list", {
  test_data <- data.frame(
    protein = rep("P1", 3),
    drug = rep("DMSO", 3),
    dose = c(0, 0, 0),
    response = c(20, 19.9, 20.1)
  )

  results <- doseResponseFit(data = test_data)

  # Should return empty data frame with correct structure
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 0)
  expect_true("Protein" %in% names(results))
  expect_true("drug" %in% names(results))
  expect_true("P_value" %in% names(results))
})

test_that("doseResponseFit handles ratio_response parameter", {
  test_data <- create_test_protein_data()

  results_log <- doseResponseFit(
    data = test_data,
    ratio_response = FALSE
  )

  results_ratio <- doseResponseFit(
    data = test_data,
    ratio_response = TRUE
  )

  # F-statistics should be different
  expect_false(all(results_log$F_statistic == results_ratio$F_statistic))
})



# Tests for fitIsotonicRegression
test_that("fitIsotonicRegression returns isotonic_model object", {
  x <- c(0, 1e-9, 1e-8, 1e-7, 1e-6)
  y <- c(20, 19, 17, 15, 12)

  fit <- fitIsotonicRegression(
    x = x,
    y = y,
    increasing = FALSE,
    transform_x = TRUE,
    ratio_y = FALSE,
    test_significance = FALSE
  )

  expect_s3_class(fit, "isotonic_model")
  expect_equal(length(fit$y_pred), length(y))
})

test_that("fitIsotonicRegression enforces monotonicity", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(10, 8, 9, 5, 3)  # Non-monotonic

  fit <- fitIsotonicRegression(
    x = x,
    y = y,
    increasing = FALSE,
    transform_x = FALSE,
    test_significance = FALSE
  )

  # Check fitted values are non-increasing
  diffs <- diff(fit$y_pred[order(fit$x)])
  expect_true(all(diffs <= 0 | abs(diffs) < 1e-10))
})

test_that("fitIsotonicRegression handles increasing constraint", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(3, 5, 4, 8, 10)

  fit <- fitIsotonicRegression(
    x = x,
    y = y,
    increasing = TRUE,
    transform_x = FALSE,
    test_significance = FALSE
  )

  # Check fitted values are non-decreasing
  diffs <- diff(fit$y_pred[order(fit$x)])
  expect_true(all(diffs >= 0 | abs(diffs) < 1e-10))
})

test_that("fitIsotonicRegression transforms dose correctly", {
  x <- c(0, 1, 10, 100, 1000)
  y <- c(20, 18, 16, 14, 12)

  fit_transform <- fitIsotonicRegression(
    x = x,
    y = y,
    transform_x = TRUE,
    test_significance = FALSE
  )

  fit_no_transform <- fitIsotonicRegression(
    x = x,
    y = y,
    transform_x = FALSE,
    test_significance = FALSE
  )

  # Internal x values should be different
  expect_false(all(fit_transform$x == fit_no_transform$x))
})

test_that("fitIsotonicRegression handles ratio transformation", {
  x <- c(0, 0, 0, 1e-6, 1e-6, 1e-6)
  y <- c(20, 20.1, 19.9, 18, 18.1, 17.9)

  fit <- fitIsotonicRegression(
    x = x,
    y = y,
    ratio_y = TRUE,
    test_significance = FALSE
  )

  # y_pred should be on ratio scale (around 0-1)
  expect_true(all(fit$y_pred >= 0))
  expect_true(max(fit$y_pred) <= 2)  # Allow some variation
})


test_that("fitIsotonicRegression handles edge cases", {
  # Single point
  expect_error(
    fitIsotonicRegression(x = 1, y = 1),
    class = "insufficient_data_error"
  )

  # Two points
  fit <- fitIsotonicRegression(
    x = c(0, 1),
    y = c(10, 5),
    test_significance = FALSE
  )
  expect_equal(length(fit$y_pred), 2)

  # All same y values
  fit <- fitIsotonicRegression(
    x = c(0, 1, 2, 3),
    y = rep(10, 4),
    test_significance = TRUE
  )
  expect_true(all(fit$y_pred == 10))
})

test_that("fitIsotonicRegression preserves original data", {
  x <- c(3, 1, 4, 2, 0)  # Unsorted
  y <- c(15, 18, 12, 16, 20)

  fit <- fitIsotonicRegression(
    x = x,
    y = y,
    transform_x = FALSE,
    ratio_y = FALSE,
    test_significance = FALSE
  )

  expect_equal(length(fit$original_x), length(x))
  expect_equal(length(fit$original_y), length(y))
})

test_that("fitIsotonicRegression handles NA values", {
  x <- c(0, 1, NA, 3, 4)
  y <- c(10, 8, 6, 4, 2)

  # Should handle or error appropriately
  expect_error(
    fitIsotonicRegression(x, y),
    NA  # Expect it to handle without error or with specific error
  )
})

test_that("fitIsotonicRegression returns NA (not NaN) when the F-test is undefined", {
  # One observation per dose -> no residual df (df_full = 0). With a perfect
  # monotone fit, SSE_Full = 0 and the old formula gave 0/0 = NaN.
  fit <- fitIsotonicRegression(x = c(0, 1e-7, 2.8e-5), y = c(4.5, 2.2, 1.8),
                               increasing = FALSE, test_significance = TRUE)
  expect_true(is.na(fit$f_test$F_statistic))
  expect_false(is.nan(fit$f_test$F_statistic))
  expect_true(is.na(fit$f_test$P_value))
})

test_that("doseResponseFit increasing = 'both' does not fail on an undefined F-statistic", {
  # P_sparse is measured in one well per dose (like a protein seen in only a
  # few wells), so its F is undefined in both directions. This used to stop the
  # whole run with "missing value where TRUE/FALSE needed".
  test_data <- rbind(
    create_test_protein_data(),
    data.frame(protein = "P_sparse", drug = c("DMSO", "Drug1", "Drug1"),
               dose = c(0, 1e-9, 1e-7), response = c(4.5, 2.2, 1.8))
  )

  results <- doseResponseFit(data = test_data, increasing = "both")

  expect_setequal(results$Protein, c("P1", "P2", "P3", "P_sparse"))
  expect_equal(sum(results$Protein == "P_sparse"), 1)
  sparse <- results[results$Protein == "P_sparse", ]
  expect_true(is.na(sparse$F_statistic))
  expect_true(is.na(sparse$P_value))
  expect_false("F_statistic_alt" %in% names(results))

  # Proteins with a defined F still get a direction and a p-value
  expect_true(all(!is.na(results$P_value[results$Protein != "P_sparse"])))
  expect_equal(results$direction[results$Protein == "P1"], "decreasing")
})

test_that("doseResponseFit empty results carry the same columns as populated ones", {
  # doseResponseFit has two separate early returns for "nothing to report".
  # A column rename once updated the populated path and only one of the two,
  # so an empty result disagreed with a populated one. Pin all three together.
  populated <- doseResponseFit(data = create_test_protein_data())
  expect_gt(nrow(populated), 0)

  # Early return 1: every drug is DMSO, so the drug list comes out empty
  no_drugs <- doseResponseFit(data = data.frame(
    protein = rep("P1", 3),
    drug = rep("DMSO", 3),
    dose = c(0, 0, 0),
    response = c(20, 19.9, 20.1)
  ))

  # Early return 2: a real drug is present, but no protein yields a fit
  no_fits <- suppressWarnings(doseResponseFit(data = data.frame(
    protein = "P1",
    drug = "D1",
    dose = 1e-9,
    response = 20
  )))

  expect_equal(nrow(no_drugs), 0)
  expect_equal(nrow(no_fits), 0)

  # Column order differs (the populated path moves `direction` after `drug`),
  # so compare as sets
  expect_setequal(names(no_drugs), names(populated))
  expect_setequal(names(no_fits), names(populated))

  # Guard the specific fields that drifted
  expect_true("Protein" %in% names(no_fits))
  expect_false("protein" %in% names(no_fits))
  expect_true("log2FC" %in% names(no_fits))
})

test_that("fitIsotonicRegression pools across a block of zero-weight points", {
  x <- c(1, 2, 3, 4, 5, 6)
  y <- c(0.1, 0.9, 0.5, 0.6, 0.7, 0.8)
  w <- c(1, 0, 0, 1, 1, 1)
  fit <- fitIsotonicRegression(x, y, w = w, increasing = TRUE,
                               transform_x = FALSE, test_significance = TRUE)
  expect_true(all(is.finite(fit$y_pred)))

  # Same curve as fitting the weighted points alone, interpolated in between
  keep <- w > 0
  fit_kept <- fitIsotonicRegression(x[keep], y[keep], w = w[keep],
                                    increasing = TRUE, transform_x = FALSE)
  expect_equal(fit$y_pred[keep], fit_kept$y_pred, tolerance = 1e-12)
  expect_equal(fit$y_pred[2:3], approx(x[keep], fit_kept$y_pred, xout = 2:3)$y,
               tolerance = 1e-12)
})

test_that("a zero-weight point does not set the fitted maximum", {
  fit <- fitIsotonicRegression(x = c(1, 2, 3, 4), y = c(0.2, 0.4, 0.6, 1.5),
                               w = c(1, 1, 1, 0), increasing = TRUE,
                               transform_x = FALSE)
  expect_equal(max(fit$y_pred), 0.6, tolerance = 1e-12)
})

test_that("fitIsotonicRegression needs 2 doses with weight > 0", {
  expect_error(
    fitIsotonicRegression(x = c(1, 1, 2, 3), y = c(0.1, 0.2, 0.3, 0.4),
                          w = c(1, 1, 0, 0), increasing = TRUE,
                          transform_x = FALSE),
    class = "insufficient_data_error"
  )
  expect_error(
    fitIsotonicRegression(x = c(1, 2, 3), y = c(0.1, 0.2, 0.3),
                          w = c(0, 0, 0), increasing = TRUE,
                          transform_x = FALSE),
    class = "insufficient_data_error"
  )
})

test_that("doseResponseFit skips proteins without enough weighted doses", {
  data <- data.frame(
    protein = rep(c("P_ok", "P_zero"), each = 4),
    drug = "Synthesis",
    dose = rep(c(1, 2, 3, 4), 2),
    response = c(0.1, 0.3, 0.5, 0.7, 0.1, 0.3, 0.5, 0.7)
  )
  weights <- c(1, 1, 1, 1, 0, 0, 0, 0)
  expect_no_warning(
    expect_message(
      results <- doseResponseFit(data, weights = weights, increasing = TRUE,
                                 transform_dose = FALSE,
                                 precalculated_ratios = TRUE),
      "1 protein\\(s\\) not fit"
    )
  )
  expect_equal(results$Protein, "P_ok")
})
