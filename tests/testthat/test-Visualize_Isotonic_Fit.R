# test_Visualize_Isotonic_Fit.R
# Unit tests for Visualize_Isotonic_Fit.R functions

library(testthat)
library(MSstatsResponse)
library(ggplot2)

# Create test data
create_test_data <- function() {
  data.frame(
    protein = rep(c("P1", "P2"), each = 8),
    drug = rep(c("DMSO", "DMSO", "DMSO", "Drug1", "Drug1", "Drug1", "Drug1", "Drug1"), 2),
    dose = rep(c(0, 0, 0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5), 2),
    response = c(20, 20.1, 19.9, 19, 18, 16, 14, 12,
                 21, 20.9, 21.1, 20.5, 20, 19.5, 19, 18.5)
  )
}

test_that("VisualizeResponseProtein returns ggplot object", {
  test_data <- create_test_data()

  plot <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  expect_s3_class(plot, "ggplot")
})

test_that("VisualizeResponseProtein handles ratio_response parameter", {
  test_data <- create_test_data()

  plot_ratio <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    ratio_response = TRUE,
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  plot_log <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    ratio_response = FALSE,
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  expect_s3_class(plot_ratio, "ggplot")
  expect_s3_class(plot_log, "ggplot")
})

test_that("VisualizeResponseProtein filters data correctly", {
  test_data <- create_test_data()
  test_data$protein[1] <- "P3"  # Add another protein

  # This should only use P1 and Drug1/DMSO data
  plot <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  expect_s3_class(plot, "ggplot")
})

test_that("VisualizeResponseProtein handles IC50 display", {
  test_data <- create_test_data()

  plot_with_ic50 <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = TRUE,
    add_ci = FALSE
  )

  expect_s3_class(plot_with_ic50, "ggplot")
})

test_that("VisualizeResponseProtein handles custom y_lab", {
  test_data <- create_test_data()

  custom_label <- "Custom Y Label"
  plot <- VisualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE,
    y_lab = custom_label
  )

  expect_s3_class(plot, "ggplot")
})

test_that("VisualizeResponseProtein handles missing protein", {
  test_data <- create_test_data()

  # Test with non-existent protein (will result in empty data)
  expect_error(
    plot <- VisualizeResponseProtein(
      data = test_data,
      protein_name = "NonExistent",
      drug_name = "Drug1",
      show_ic50 = FALSE,
      add_ci = FALSE
    ),
    "need at least two"  # Expect the specific error
  )
})

# Tests for plotIsotonic (if exported)
test_that("plotIsotonic creates valid plot from fit object", {
  # Create mock fit object
  mock_fit <- list(
    x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    y_pred = c(1.0, 0.9, 0.7, 0.4, 0.2),
    original_x = c(0, 1e-9, 1e-8, 1e-7, 1e-6),
    original_y = c(1.0, 0.88, 0.72, 0.38, 0.18),
    increasing = FALSE,
    transform_x = TRUE,
    ratio_y = TRUE
  )
  class(mock_fit) <- "isotonic_model"

  plot <- plotIsotonic(
    fit = mock_fit,
    ratio = TRUE,
    show_ic50 = FALSE,
    drug_name = "TestDrug",
    protein_name = "TestProtein"
  )

  expect_s3_class(plot, "ggplot")
})
