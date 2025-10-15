# test_IC50_prediction.R
# Unit tests for IC50_prediction.R functions

library(testthat)
library(MSstatsResponse)

# Helper function to create test data
create_dose_response_data <- function(seed = 123) {
  set.seed(seed)
  doses <- c(0, 0, 0, 1e-9, 1e-9, 1e-9, 1e-8, 1e-8, 1e-8,
             1e-7, 1e-7, 1e-7, 1e-6, 1e-6, 1e-6)
  responses <- c(25, 24, 26, 22, 23, 22, 18, 19, 17,
                 15, 14, 16, 12, 13, 11)
  return(list(dose = doses, response = responses))
}

# Tests for PredictIC50
test_that("PredictIC50 finds IC50 when response crosses target", {
  # Create a simple monotonic fit
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(1.0, 0.8, 0.6, 0.4, 0.2)
  )

  ic50 <- PredictIC50(mock_fit, target_response = 0.5)

  # Should be between 1e-8 and 1e-7
  expect_true(ic50 > 1e-8 && ic50 < 1e-7)
})

test_that("PredictIC50 returns exact dose when response matches target", {
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(1.0, 0.8, 0.5, 0.3, 0.1)  # Exactly 0.5 at 1e-8
  )

  ic50 <- PredictIC50(mock_fit, target_response = 0.5)
  expect_equal(ic50, 1e-8)
})

test_that("PredictIC50 handles different target responses", {
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(1.0, 0.8, 0.6, 0.4, 0.2)
  )

  ic25 <- PredictIC50(mock_fit, target_response = 0.25)
  ic50 <- PredictIC50(mock_fit, target_response = 0.50)
  ic75 <- PredictIC50(mock_fit, target_response = 0.75)

  # IC75 should be lower dose than IC50 than IC25 (for decreasing response)
  expect_true(ic75 < ic50)
  expect_true(ic50 < ic25)
})

test_that("PredictIC50 returns NA when target not reached", {
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(0.9, 0.8, 0.7, 0.6, 0.55)  # Never reaches 0.5
  )

  expect_warning(
    ic50 <- PredictIC50(mock_fit, target_response = 0.5),
    "Target response value not reached"
  )
  expect_true(is.na(ic50))
})

test_that("PredictIC50 handles increasing response", {
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(0.2, 0.4, 0.6, 0.8, 1.0)  # Increasing
  )

  ic50 <- PredictIC50(mock_fit, target_response = 0.5)
  expect_true(ic50 > 1e-9 && ic50 < 1e-8)
})

# Tests for bootstrapIC50
test_that("bootstrapIC50 returns correct structure", {
  data <- create_dose_response_data()

  result <- bootstrapIC50(
    dose = data$dose,
    response = data$response,
    n_samples = 100,  # Fewer samples for speed
    alpha = 0.10,
    increasing = FALSE,
    target_response = 0.5
  )

  expect_type(result, "list")
  expect_named(result, c("ic50_values", "mean_ic50", "ci_lower", "ci_upper",
                         "mean_ic50_transform", "ci_lower_transform", "ci_upper_transform"))
  expect_length(result$ic50_values, 100)
})

test_that("bootstrapIC50 confidence intervals contain mean", {
  data <- create_dose_response_data()

  result <- bootstrapIC50(
    dose = data$dose,
    response = data$response,
    n_samples = 100,
    alpha = 0.10,
    increasing = FALSE
  )

  # Mean should be within confidence intervals
  expect_true(result$mean_ic50 >= result$ci_lower)
  expect_true(result$mean_ic50 <= result$ci_upper)
})



test_that("bootstrapIC50 handles missing values", {
  doses <- c(0, 0, 1e-8, 1e-7, 1e-6, NA)
  responses <- c(25, 24, 18, 15, 12, 20)

  # Should handle NA in dose
  result <- bootstrapIC50(
    dose = doses,
    response = responses,
    n_samples = 50
  )

  expect_type(result$mean_ic50, "double")
})

# Tests for bootstrapIC50LogScale
test_that("bootstrapIC50LogScale returns correct structure", {
  data <- create_dose_response_data()

  result <- bootstrapIC50LogScale(
    x = data$dose,
    y = data$response,
    n_samples = 100,
    alpha = 0.05,
    increasing = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("ic50_values", "mean_ic50", "ci_lower", "ci_upper",
                         "mean_ic50_transform", "ci_lower_transform", "ci_upper_transform"))
})

test_that("bootstrapIC50LogScale handles no DMSO samples", {
  # Data without dose = 0
  doses <- c(1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7)
  responses <- c(22, 23, 18, 19, 15, 14)

  result <- bootstrapIC50LogScale(
    x = doses,
    y = responses,
    n_samples = 50,
    increasing = FALSE
  )

  # Should have NAs since no DMSO
  expect_true(sum(!is.na(result$ic50_values)) < 50)
})

test_that("bootstrapIC50LogScale adjusts target response correctly", {
  data <- create_dose_response_data()

  result <- bootstrapIC50LogScale(
    x = data$dose,
    y = data$response,
    n_samples = 100,
    target_response = 0.25  # IC25 instead of IC50
  )

  expect_type(result$mean_ic50, "double")
  expect_false(is.na(result$mean_ic50))
})

# Tests for helper functions (.calcSingleIC50)
test_that(".calcSingleIC50 handles single protein-drug pair", {
  test_df <- data.frame(
    dose = c(0, 0, 0, 1e-9, 1e-9, 1e-8, 1e-8, 1e-7, 1e-7),
    response = c(20, 20.1, 19.9, 19, 18.9, 17, 17.1, 15, 14.9)
  )

  result <- .calcSingleIC50(
    df = test_df,
    n_samples = 50,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = TRUE,
    bootstrap = FALSE,
    prot = "TestProtein",
    drug_type = "TestDrug",
    target_response = 0.5
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("protein", "drug", "IC50", "IC50_lower_bound", "IC50_upper_bound"))
  expect_equal(result$protein, "TestProtein")
  expect_equal(result$drug, "TestDrug")
})

test_that(".calcSingleIC50 handles bootstrap option", {
  test_df <- data.frame(
    dose = c(0, 0, 1e-9, 1e-8, 1e-7),
    response = c(20, 19.9, 18, 16, 14)
  )

  # Without bootstrap
  result_no_boot <- .calcSingleIC50(
    df = test_df,
    n_samples = 50,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = TRUE,
    bootstrap = FALSE,
    prot = "P1",
    drug_type = "D1",
    target_response = 0.5
  )

  expect_true(is.na(result_no_boot$IC50_lower_bound))
  expect_true(is.na(result_no_boot$IC50_upper_bound))

  # With bootstrap
  result_boot <- .calcSingleIC50(
    df = test_df,
    n_samples = 50,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = TRUE,
    bootstrap = TRUE,
    prot = "P1",
    drug_type = "D1",
    target_response = 0.5
  )

  expect_false(is.na(result_boot$IC50_lower_bound))
  expect_false(is.na(result_boot$IC50_upper_bound))
})

test_that(".calcSingleIC50 handles ratio vs log scale", {
  test_df <- data.frame(
    dose = c(0, 0, 1e-9, 1e-8, 1e-7),
    response = c(20, 19.9, 18, 16, 14)
  )

  result_ratio <- .calcSingleIC50(
    df = test_df,
    n_samples = 30,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = TRUE,
    bootstrap = FALSE,
    prot = "P1",
    drug_type = "D1",
    target_response = 0.5
  )

  result_log <- .calcSingleIC50(
    df = test_df,
    n_samples = 30,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = FALSE,
    bootstrap = FALSE,
    prot = "P1",
    drug_type = "D1",
    target_response = 0.5
  )

  # Results should be different
  expect_false(isTRUE(all.equal(result_ratio$IC50, result_log$IC50)))
})

test_that(".calcSingleIC50 handles failed fits", {
  # Data that might cause fitting issues
  test_df <- data.frame(
    dose = c(0, 1e-9),
    response = c(20, 20)  # Flat response
  )

  suppressWarnings({
  result <- .calcSingleIC50(
    df = test_df,
    n_samples = 10,
    alpha = 0.10,
    increasing = FALSE,
    transform_dose = TRUE,
    ratio_response = TRUE,
    bootstrap = FALSE,
    prot = "P1",
    drug_type = "D1",
    target_response = 0.5
  )
  })

  expect_true(is.na(result$IC50))
})
