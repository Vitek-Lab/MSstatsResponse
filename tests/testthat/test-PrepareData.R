# test_PrepareData.R
# Unit tests for PrepareData.R functions

library(testthat)
library(MSstatsResponse)

# Tests for convertGroupToNumericDose
test_that("convertGroupToNumericDose handles DMSO correctly", {
  result <- convertGroupToNumericDose("DMSO")
  expect_equal(result$dose_nM, 0)
  expect_equal(result$drug, "DMSO")
})

test_that("convertGroupToNumericDose converts uM to nM correctly", {
  groups <- c("Dasatinib_001uM", "Dasatinib_010uM", "Dasatinib_100uM")
  result <- convertGroupToNumericDose(groups)
  expect_equal(result$dose_nM, c(1000, 10000, 100000))
  expect_equal(unique(result$drug), "Dasatinib")
})

test_that("convertGroupToNumericDose handles nM units correctly", {
  groups <- c("Imatinib_100nM", "Imatinib_1000nM")
  result <- convertGroupToNumericDose(groups)
  expect_equal(result$dose_nM, c(100, 1000))
  expect_equal(unique(result$drug), "Imatinib")
})

test_that("convertGroupToNumericDose handles mixed drugs", {
  groups <- c("DMSO", "Dasatinib_001uM", "Imatinib_100nM")
  result <- convertGroupToNumericDose(groups)
  expect_equal(result$dose_nM, c(0, 1000, 100))
  expect_equal(result$drug, c("DMSO", "Dasatinib", "Imatinib"))
})

test_that("convertGroupToNumericDose handles decimal values", {
  groups <- c("Drug_0.5uM", "Drug_1.5uM")
  result <- convertGroupToNumericDose(groups)
  expect_equal(result$dose_nM, c(500, 1500))
})

test_that("convertGroupToNumericDose returns data.frame with correct structure", {
  result <- convertGroupToNumericDose(c("DMSO", "Drug_001uM"))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("drug", "dose_nM"))
  expect_type(result$dose_nM, "double")
  expect_type(result$drug, "character")
})

# Tests for MSstatsPrepareDoseResponseFit
test_that("MSstatsPrepareDoseResponseFit handles basic input correctly", {
  test_data <- data.frame(
    Protein = c("P1", "P1", "P2", "P2"),
    drug = c("DMSO", "Drug1", "DMSO", "Drug1"),
    dose = c(0, 1e-6, 0, 1e-6),
    LogIntensities = c(20, 18, 21, 19)
  )

  result <- MSstatsPrepareDoseResponseFit(test_data)

  expect_equal(nrow(result), 4)
  expect_named(result, c("protein", "drug", "dose", "response"))
  expect_equal(result$protein, c("P1", "P1", "P2", "P2"))
  expect_equal(result$response, c(20, 18, 21, 19))
})

test_that("MSstatsPrepareDoseResponseFit handles custom column names", {
  test_data <- data.frame(
    ProteinID = c("P1", "P2"),
    Treatment = c("DMSO", "Drug1"),
    Concentration = c(0, 100),
    Log2Abundance = c(20, 18)
  )

  result <- MSstatsPrepareDoseResponseFit(
    test_data,
    dose_column = "Concentration",
    drug_column = "Treatment",
    protein_column = "ProteinID",
    log_abundance_column = "Log2Abundance"
  )

  expect_named(result, c("protein", "drug", "dose", "response"))
  expect_equal(result$protein, c("P1", "P2"))
  expect_equal(result$dose, c(0, 100))
})

test_that("MSstatsPrepareDoseResponseFit handles nM to M conversion", {
  test_data <- data.frame(
    Protein = c("P1", "P1"),
    drug = c("DMSO", "Drug1"),
    dose_nM = c(0, 1000),
    LogIntensities = c(20, 18)
  )

  result <- MSstatsPrepareDoseResponseFit(
    test_data,
    dose_column = "dose_nM",
    transform_nM_to_M = TRUE
  )

  expect_equal(result$dose, c(0, 1e-6))
  expect_equal(result$dose_nM, c(0, 1000))
})

test_that("MSstatsPrepareDoseResponseFit preserves data when transform_nM_to_M is FALSE", {
  test_data <- data.frame(
    Protein = "P1",
    drug = "Drug1",
    dose = 1e-6,
    LogIntensities = 20
  )

  result <- MSstatsPrepareDoseResponseFit(
    test_data,
    transform_nM_to_M = FALSE
  )

  expect_equal(result$dose, 1e-6)
  expect_false("dose_nM" %in% names(result))
})

test_that("MSstatsPrepareDoseResponseFit throws error for missing columns", {
  test_data <- data.frame(
    Protein = c("P1", "P2"),
    dose = c(0, 1)
  )

  expect_error(
    MSstatsPrepareDoseResponseFit(test_data),
    "Missing required column"
  )
})

test_that("MSstatsPrepareDoseResponseFit validates input is data.frame", {
  expect_error(
    MSstatsPrepareDoseResponseFit(list(a = 1, b = 2)),
    class = "simpleError"
  )
})
