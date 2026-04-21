# test_ExperimentalDesignSimulation.R
# Unit tests for ExperimentalDesignSimulation.R functions

library(testthat)
library(MSstatsResponse)

# Helper to check if template structure is valid
validate_template_structure <- function(template) {
  required_cols <- c("dose", "Intensity", "LogIntensities", "ratio")
  all(sapply(template, function(x) all(required_cols %in% names(x))))
}

# Tests for futureExperimentSimulation
test_that("futureExperimentSimulation returns correct structure", {
  sim_results <- futureExperimentSimulation(
    N_proteins = 30,  # Small for speed
    N_rep = 2,
    N_Control_Rep = 3,
    Concentrations = c(0, 10, 100, 1000),
    IC50_Prediction = FALSE
  )

  expect_type(sim_results, "list")
  expect_named(sim_results, c("Simulated_Data", "MSstats_Simulated_Data",
                              "DoseResponseFit_Results", "Hit_Rates_Plot",
                              "Hit_Rates_Data", "PredictIC50_Results",
                              "Template_Used"))

  expect_s3_class(sim_results$Simulated_Data, "data.frame")
  expect_s3_class(sim_results$MSstats_Simulated_Data, "data.frame")
  expect_s3_class(sim_results$DoseResponseFit_Results, "data.frame")
  expect_s3_class(sim_results$Hit_Rates_Plot, "ggplot")
})


test_that("futureExperimentSimulation handles custom control replicates", {
  sim_results <- futureExperimentSimulation(
    N_proteins = 30,
    N_rep = 2,
    N_Control_Rep = 5,  # Different from N_rep
    Concentrations = c(0, 100, 1000),
    IC50_Prediction = FALSE
  )

  control_data <- sim_results$Simulated_Data[sim_results$Simulated_Data$dose == 0, ]
  n_control_per_protein <- table(control_data$protein)

  expect_true(all(n_control_per_protein == 5))
})



test_that("futureExperimentSimulation uses custom templates from user data", {
  # Create mock user data
  mock_data <- data.frame(
    protein = rep(c("P1", "P2", "P3"), each = 4),
    drug = rep(c("DMSO", "DMSO", "Drug1", "Drug1"), 3),
    dose = rep(c(0, 0, 1, 1000), 3),  # nM scale to match Concentrations
    response = c(20, 20, 18, 15,  # P1
                 21, 21, 20, 19,  # P2
                 19, 19, 19, 19)  # P3
  )

  sim_results <- futureExperimentSimulation(
    N_proteins = 30,
    N_rep = 2,
    Concentrations = c(0, 1, 1000),  # nM scale
    data = mock_data,
    strong_proteins = "P1",
    weak_proteins = "P2",
    no_interaction_proteins = "P3",
    drug_name = "Drug1",
    IC50_Prediction = FALSE
  )

  expect_type(sim_results$Template_Used, "list")
  expect_true(validate_template_structure(sim_results$Template_Used))
})

test_that("futureExperimentSimulation requires protein specification with custom data", {
  mock_data <- data.frame(
    protein = "P1",
    drug = c("DMSO", "Drug1"),
    dose = c(0, 1e-6),
    response = c(20, 18)
  )

  expect_error(
    futureExperimentSimulation(
      N_proteins = 30,
      data = mock_data
      # No protein specifications
    ),
    "must specify at least one"
  )
})

# Tests for .extractTemplatesFromData
test_that(".extractTemplatesFromData creates valid template structure", {
  mock_data <- data.frame(
    protein = rep(c("P1", "P2"), each = 4),
    drug = rep(c("DMSO", "DMSO", "Drug1", "Drug1"), 2),
    dose = rep(c(0, 0, 1e-9, 1e-6), 2),
    response = c(20, 20.1, 18, 15,  # P1
                 21, 20.9, 20, 19)   # P2
  )

  template <- .extractTemplatesFromData(
    data = mock_data,
    strong_proteins = "P1",
    weak_proteins = "P2",
    no_interaction_proteins = NULL,
    drug_name = "Drug1",
    concentrations = c(0, 1, 1000)  # nM
  )

  expect_type(template, "list")
  expect_named(template, c("strong_interaction", "weak_interaction", "no_interaction"))
  expect_true(validate_template_structure(template))
})

test_that(".extractTemplatesFromData handles missing doses via interpolation", {
  mock_data <- data.frame(
    protein = rep("P1", 4),
    drug = c("DMSO", "DMSO", "Drug1", "Drug1"),
    dose = c(0, 0, 1e-7, 1e-5),  # Missing middle doses
    response = c(20, 20, 15, 10)
  )

  template <- .extractTemplatesFromData(
    data = mock_data,
    strong_proteins = "P1",
    weak_proteins = NULL,
    no_interaction_proteins = NULL,
    drug_name = "Drug1",
    concentrations = c(0, 10, 100, 1000, 3000)  # All within default template range
  )

  # Should have all requested concentrations
  expect_equal(nrow(template$strong_interaction), 5)
  expect_equal(template$strong_interaction$dose, c(0, 10, 100, 1000, 3000))
})

test_that(".extractTemplatesFromData handles empty protein lists", {
  mock_data <- data.frame(
    protein = "P1",
    drug = c("DMSO", "Drug1"),
    dose = c(0, 1e-6),
    response = c(20, 18)
  )

  template <- .extractTemplatesFromData(
    data = mock_data,
    strong_proteins = "P1",
    weak_proteins = NULL,
    no_interaction_proteins = NULL,
    drug_name = "Drug1",
    concentrations = c(0, 1000)
  )

  # Null protein lists get flat baseline profiles
  expect_equal(nrow(template$weak_interaction), 2)
  expect_equal(nrow(template$no_interaction), 2)
  expect_true(all(c("dose", "Intensity", "LogIntensities", "ratio") %in% names(template$weak_interaction)))
  expect_true(all(c("dose", "Intensity", "LogIntensities", "ratio") %in% names(template$no_interaction)))
  # Flat profiles should have ratio = 1 for all doses
  expect_equal(template$weak_interaction$ratio, c(1, 1))
  expect_equal(template$no_interaction$ratio, c(1, 1))
})


test_that("simulateChemoProteinLevelNonParametric subsets to matching concentrations", {
  # Template only has dose = 0; concentrations c(0, 5, 15) includes non-matching
  # values. The function silently subsets rather than erroring — verify it runs
  # and returns only the matching dose.
  template <- list(
    strong_interaction = data.frame(dose = 0, LogIntensities = 20, Intensity = 2^20, ratio = 1),
    weak_interaction   = data.frame(dose = 0, LogIntensities = 20, Intensity = 2^20, ratio = 1),
    no_interaction     = data.frame(dose = 0, LogIntensities = 20, Intensity = 2^20, ratio = 1)
  )

  result <- simulateChemoProteinLevelNonParametric(
    N_proteins     = 30,
    concentrations = c(0, 5, 15),
    template       = template,
    outlier_prob   = 0
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(result$dose == 0))  # Only dose = 0 was in the template
})


test_that("simulateChemoProteinLevelNonParametric adds outliers", {
  template <- list(
    strong_interaction = data.frame(
      dose = 0,
      LogIntensities = 20,
      Intensity = 2^20,
      ratio = 1
    ),
    weak_interaction = data.frame(
      dose = 0,
      LogIntensities = 20,
      Intensity = 2^20,
      ratio = 1
    ),
    no_interaction = data.frame(
      dose = 0,
      LogIntensities = 20,
      Intensity = 2^20,
      ratio = 1
    )
  )

  sim_data <- simulateChemoProteinLevelNonParametric(
    N_proteins = 100,
    concentrations = c(0),
    rep = 3,
    template = template,
    outlier_prob = 0.1
  )

  expect_true("outlier" %in% names(sim_data))
  expect_true(sum(sim_data$outlier) > 0)
})

# Tests for plotHitRateMSstatsResponse
test_that("plotHitRateMSstatsResponse calculates hit rates correctly", {
  # Create mock results with known outcomes
  mock_results <- data.frame(
    Protein = c(paste0("p_strong_interaction_", 1:10),
                paste0("p_weak_interaction_", 1:10),
                paste0("p_no_interaction_", 1:10)),
    drug = "Drug1",
    adj.pvalue = c(rep(0.01, 8), rep(0.1, 2),  # 80% strong detected
                    rep(0.01, 5), rep(0.1, 5),   # 50% weak detected
                    rep(0.01, 1), rep(0.1, 9))   # 10% false positive
  )

  plot_result <- plotHitRateMSstatsResponse(
    mock_results,
    rep_count = 3,
    concentration_count = 5
  )

  expect_s3_class(plot_result$plot, "ggplot")
  expect_equal(plot_result$plot_data$Percent[1], 80)  # TPR Strong
  expect_equal(plot_result$plot_data$Percent[2], 50)  # TPR Weak
  expect_equal(plot_result$plot_data$Percent[3], 10)  # FPR
})

test_that("plotHitRateMSstatsResponse handles empty groups", {
  # Results with no weak proteins
  mock_results <- data.frame(
    Protein = c(paste0("p_strong_interaction_", 1:5),
                paste0("p_no_interaction_", 1:5)),
    drug = "Drug1",
    adj.pvalue = c(rep(0.01, 3), rep(0.1, 2),
                    rep(0.1, 5))
  )

  plot_result <- plotHitRateMSstatsResponse(
    mock_results,
    rep_count = 3,
    concentration_count = 5
  )

  expect_true(is.na(plot_result$plot_data$Percent[2]))  # No weak proteins
})

test_that("futureExperimentSimulation achieves expected TPR/FPR with built-in data", {
  # Prepare built-in DIA data the same way as the vignette
  data("DIA_MSstats_Normalized", package = "MSstatsResponse")
  dia_normalized <- DIA_MSstats_Normalized
  dose_info <- convertGroupToNumericDose(dia_normalized$ProteinLevelData$GROUP)
  dia_normalized$ProteinLevelData$dose_nM <- dose_info$dose_nM
  dia_normalized$ProteinLevelData$drug    <- dose_info$drug

  dia_prepared <- MSstatsPrepareDoseResponseFit(
    data                 = dia_normalized$ProteinLevelData,
    dose_column          = "dose_nM",
    drug_column          = "drug",
    protein_column       = "Protein",
    log_abundance_column = "LogIntensities",
    transform_nM_to_M    = TRUE
  )

  dose_levels <- as.numeric(unique(dia_prepared$dose))

  set.seed(42)
  custom_simulation <- futureExperimentSimulation(
    N_proteins              = 300,
    N_rep                   = 3,
    N_Control_Rep           = 3,
    Concentrations          = dose_levels,
    data                    = dia_prepared,
    strong_proteins         = "PROTEIN_B",
    weak_proteins           = "PROTEIN_A",
    no_interaction_proteins = "PROTEIN_C",
    drug_name               = "Drug1",
    IC50_Prediction         = FALSE
  )

  hit_rates <- custom_simulation$Hit_Rates_Data

  tpr_strong <- hit_rates$Percent[hit_rates$Category == "TPR (Strong)"]
  tpr_weak   <- hit_rates$Percent[hit_rates$Category == "TPR (Weak)"]
  fpr        <- hit_rates$Percent[hit_rates$Category == "FPR"]

  expect_gte(tpr_strong, 80)
  expect_gte(tpr_weak,   80)
  expect_lt(fpr,         10)
})
