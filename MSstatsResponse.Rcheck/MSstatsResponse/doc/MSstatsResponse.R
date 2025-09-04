## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----install, eval=FALSE------------------------------------------------------
# # Install from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MSstatsResponse")

## ----load_libraries, message=FALSE--------------------------------------------
library(MSstatsResponse)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(data.table)

# Optional: for upstream data processing
# library(MSstats)
# library(MSstatsTMT)

## ----msstats_preprocessing, eval=FALSE----------------------------------------
# # Read raw data (example with Spectronaut output)
# raw_data <- read_tsv("path/to/spectronaut_report.tsv")
# 
# # Convert to MSstats format
# msstats_data <- SpectronauttoMSstatsFormat(raw_data)
# 
# # Process data: normalization and protein summarization
# processed_data <- dataProcess(
#   msstats_data,
#   normalization = "equalizeMedians",  # or FALSE for no normalization
#   summaryMethod = "TMP",              # Tukey's median polish
#   MBimpute = TRUE,                     # Impute missing values
#   maxQuantileforCensored = 0.999
# )
# 
# # Extract protein-level data for dose-response analysis
# protein_level_data <- processed_data$ProteinLevelData

## ----load_example_data--------------------------------------------------------
# Load pre-processed DIA-MS data example
data_path <- system.file("extdata", "DIA_MSstats_Normalized.RDS", 
                         package = "MSstatsResponse")
dia_normalized <- readRDS(data_path)

# Examine data structure
str(dia_normalized$ProteinLevelData[1:5, ])

## ----convert_doses------------------------------------------------------------
# Convert GROUP labels to numeric doses
dose_info <- ConvertGroupToNumericDose(dia_normalized$ProteinLevelData$GROUP)

# Add dose and drug information to the dataset
dia_normalized$ProteinLevelData <- dia_normalized$ProteinLevelData %>%
  mutate(
    dose_nM = dose_info$dose_nM,  # Dose in nanomolar
    drug = dose_info$drug          # Drug name
  )


# View converted data
dia_normalized$ProteinLevelData %>%
  select(Protein, GROUP, drug, dose_nM) %>%
  head(10)

## ----prepare_data-------------------------------------------------------------
# Prepare data for dose-response fitting
dia_prepared <- MSstatsPrepareDoseResponseFit(
  data = dia_normalized$ProteinLevelData,
  dose_column = "dose_nM",
  drug_column = "drug",
  protein_column = "Protein",
  log_abundance_column = "LogIntensities",
  transform_nM_to_M = TRUE  
)

# Examine prepared data structure
str(dia_prepared)

# View sample of prepared data
dia_prepared %>%
  filter(protein %in% c("P00519", "P12931")) %>%
  arrange(protein, drug, dose) %>%
  head(20)

## ----fit_dose_response, message=FALSE, warning=FALSE--------------------------
# Detect drug-protein interactions
response_results <- DoseResponseFit(
  data = dia_prepared,
  increasing = FALSE,        # Expect decreasing response
  transform_dose = TRUE,     # Apply log10(dose + 1) transformation
  ratio_response = FALSE     # Stay on log2 scale for testing
)

# Examine results
response_results %>%
  select(protein, drug, F_statistic, P_value, adjust_pval) %>%
  arrange(adjust_pval) %>%
  head(10)

## ----summarize_results--------------------------------------------------------
# Count significant interactions
n_total <- nrow(response_results)
n_significant <- sum(response_results$adjust_pval < 0.05)

cat("Total protein-drug pairs tested:", n_total, "\n")
cat("Significant interactions (FDR < 0.05):", n_significant, "\n")
cat("Percentage significant:", round(100 * n_significant/n_total, 1), "%\n")

# Top hits
top_hits <- response_results %>%
  filter(adjust_pval < 0.05) %>%
  arrange(adjust_pval) %>%
  head(5)

print(top_hits)

## ----predict_ic50, message=FALSE, warning=FALSE-------------------------------
# Estimate IC50 with bootstrap confidence intervals
ic50_predictions <- PredictIC50(
  data = dia_prepared,
  ratio_response = TRUE,     # Use ratio scale for IC50
  transform_dose = TRUE,     # Log-transform doses
  increasing = FALSE,        # Decreasing response expected
  bootstrap = TRUE,          # Compute confidence intervals
  n_samples = 1000,         # Number of bootstrap samples
  alpha = 0.10              # 90% confidence intervals
)

# View IC50 estimates
ic50_predictions %>%
  arrange(IC50) %>%
  head(10)

## ----parallel_ic50, eval=FALSE------------------------------------------------
# # Parallel IC50 estimation (faster for large datasets)
# ic50_parallel <- PredictIC50Parallel(
#   data = dia_prepared,
#   ratio_response = TRUE,
#   transform_dose = TRUE,
#   bootstrap = TRUE,
#   n_samples = 1000,
#   numberOfCores = 4  # Use 4 CPU cores
# )

## ----visualize_single, fig.height=5, fig.width=7------------------------------
# Visualize strong responder
VisualizeResponseProtein(
  data = dia_prepared,
  protein_name = "P00519",
  drug_name = "Dasatinib",
  ratio_response = TRUE,
  show_ic50 = TRUE,
  add_ci = TRUE,
  n_samples = 1000
)

## ----visualize_another, fig.height=5, fig.width=7-----------------------------
# Visualize another target
VisualizeResponseProtein(
  data = dia_prepared,
  protein_name = "P12931",
  drug_name = "Dasatinib",
  ratio_response = TRUE,
  show_ic50 = TRUE,
  add_ci = TRUE
)

## ----compare_scales, fig.height=4, fig.width=12-------------------------------
# Log2 scale (left panel)
p1 <- VisualizeResponseProtein(
  data = dia_prepared,
  protein_name = "P00519",
  drug_name = "Dasatinib",
  ratio_response = FALSE,  # Log2 scale
  show_ic50 = TRUE,
  add_ci = FALSE
)

# Ratio scale (right panel)  
p2 <- VisualizeResponseProtein(
  data = dia_prepared,
  protein_name = "P00519",
  drug_name = "Dasatinib",
  ratio_response = TRUE,   # Ratio scale
  show_ic50 = TRUE,
  add_ci = FALSE
)

# Combine plots (requires gridExtra)
 gridExtra::grid.arrange(p1, p2, ncol = 2)

## ----simulation_example, eval=TRUE--------------------------------------------
# Simulate experiment with specific parameters
simulation_results <- FutureExperimentSimulation(
  N_proteins = 3000,           # Total number of proteins
  N_rep = 3,                   # Replicates per dose
  N_Control_Rep = 3,           # Control (DMSO) replicates
  Concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),  # nM
  IC50_Prediction = FALSE       # Also perform IC50 estimation
)

# View performance metrics
print(simulation_results$Hit_Rates_Plot)

## -----------------------------------------------------------------------------
# First, identify proteins from your results to use as templates
# Run simulation using your data as templates

custom_simulation <- FutureExperimentSimulation(
  N_proteins = 3000,
  N_rep = 2,
  N_Control_Rep = 3,
  Concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
  data = dia_prepared,                          # Your data
  strong_proteins = 'P12931',                # User defined strong responder
  weak_proteins = 'P00519',                    # User defined weak responder 
  no_interaction_proteins = 'Q9P2K8',          # User defined negative control 
  drug_name = "Dasatinib",                     # Specify which drug to model
  IC50_Prediction = FALSE
)

# Compare performance metrics
print(custom_simulation$Hit_Rates_Plot)

# Examine the templates that were extracted from your data
print(custom_simulation$Template_Used)

## ----power_curves, fig.height=5, fig.width=8, eval=FALSE----------------------
# # Plot power curves for strong interactions only
# power_plot_strong <- PlotExperimentPowerCurve(
#   rep_grid = 1:5,
#   N_proteins = 300,
#   interaction_type = "Strong"
# )
# print(power_plot_strong)
# 
# # Plot power curves for weak interactions only
# power_plot_weak <- PlotExperimentPowerCurve(
#   rep_grid = 1:5,
#   N_proteins = 300,
#   interaction_type = "Weak"
# )
# print(power_plot_weak)

## ----session_info-------------------------------------------------------------
sessionInfo()

