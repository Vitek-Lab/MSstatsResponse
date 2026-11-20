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

test_that("visualizeResponseProtein returns ggplot object", {
  test_data <- create_test_data()

  plot <- visualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  expect_s3_class(plot, "ggplot")
})

test_that("visualizeResponseProtein handles ratio_response parameter", {
  test_data <- create_test_data()

  plot_ratio <- visualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    ratio_response = TRUE,
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  plot_log <- visualizeResponseProtein(
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

test_that("visualizeResponseProtein filters data correctly", {
  test_data <- create_test_data()
  test_data$protein[1] <- "P3"  # Add another protein

  # This should only use P1 and Drug1/DMSO data
  plot <- visualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE
  )

  expect_s3_class(plot, "ggplot")
})

test_that("visualizeResponseProtein handles IC50 display", {
  test_data <- create_test_data()

  plot_with_ic50 <- visualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = TRUE,
    add_ci = FALSE
  )

  expect_s3_class(plot_with_ic50, "ggplot")
})

test_that("visualizeResponseProtein handles custom y_lab", {
  test_data <- create_test_data()

  custom_label <- "Custom Y Label"
  plot <- visualizeResponseProtein(
    data = test_data,
    protein_name = "P1",
    drug_name = "Drug1",
    show_ic50 = FALSE,
    add_ci = FALSE,
    y_lab = custom_label
  )

  expect_s3_class(plot, "ggplot")
})

test_that("visualizeResponseProtein handles missing protein", {
  test_data <- create_test_data()

  # Test with non-existent protein (will result in empty data)
  expect_error(
    plot <- visualizeResponseProtein(
      data = test_data,
      protein_name = "NonExistent",
      drug_name = "Drug1",
      show_ic50 = FALSE,
      add_ci = FALSE
    ),
    class = "insufficient_data_error"
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

create_turnover_plot_data <- function() {
  dose <- rep(c(1, 2, 6, 12, 24, 48), times = 2)
  data.frame(
    protein = "P1",
    drug = "Synthesis",
    dose = dose,
    response = c(0.05, 0.15, 0.35, 0.55, 0.75, 0.9,
                 0.08, 0.2, 0.3, 0.6, 0.7, 1.6),
    BaseSequence = rep(c("PEPA", "PEPB"), each = 6),
    weight = c(rep(1, 6), rep(0.5, 5), 0)
  )
}

test_that("visualizeResponseProtein supports the turnover plotting options", {
  data <- create_turnover_plot_data()
  plot <- visualizeResponseProtein(
    data = data, protein_name = "P1", drug_name = "Synthesis",
    precalculated_ratios = TRUE, ratio_response = FALSE,
    transform_dose = FALSE, weights = data$weight, increasing = TRUE,
    show_ic50 = TRUE, ic50_label = "t50", add_ci = FALSE,
    target_response = 0.5, show_weights = TRUE, color_by = "BaseSequence"
  )
  expect_s3_class(plot, "ggplot")

  built <- ggplot2::ggplot_build(plot)
  labels <- unlist(lapply(built$data, function(layer) layer$label))
  expect_true(any(grepl("t50", labels)))
  point_layer <- built$data[[which(vapply(built$plot$layers, function(l)
    inherits(l$geom, "GeomPoint"), logical(1)))[1]]]
  expect_equal(length(unique(point_layer$colour)), 2)
})

test_that("visualizeResponseProtein errors on an unknown color_by column", {
  data <- create_turnover_plot_data()
  expect_error(
    visualizeResponseProtein(
      data = data, protein_name = "P1", drug_name = "Synthesis",
      precalculated_ratios = TRUE, ratio_response = FALSE,
      transform_dose = FALSE, increasing = TRUE, show_ic50 = FALSE,
      add_ci = FALSE, color_by = "NotAColumn"
    ),
    "NotAColumn"
  )
})
