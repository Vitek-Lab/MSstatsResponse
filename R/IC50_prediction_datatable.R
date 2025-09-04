#' Parallel version of PredictIC50 function
#'
#' Runs PredictIC50 on the entire dataset in parallel across proteins.
#'
#' @param data A data frame with columns: protein, drug, dose, response.
#' @param n_samples Number of bootstrap samples. Default = 1000.
#' @param alpha Confidence level. Default = 0.10.
#' @param increasing Logical. If TRUE, fit non-decreasing trend. Default = FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default = TRUE.
#' @param ratio_response Logical. If TRUE, use ratio response; else use log2 scale. Default = TRUE.
#' @param bootstrap Logical. If TRUE, compute bootstrap CIs. Default = TRUE.
#' @param numberOfCores Number of cores for parallel processing. Default = 2.
#'
#' @return A data frame with columns: protein, drug, IC50, lower CI, upper CI.
#' @examples
#' # Load example data
#' data_path <- system.file("extdata", "DIA_MSstats_Normalized.RDS",
#'                          package = "MSstatsResponse")
#' dia_data <- readRDS(data_path)
#'
#' # Convert GROUP to dose
#' dose_info <- ConvertGroupToNumericDose(dia_data$ProteinLevelData$GROUP)
#' dia_data$ProteinLevelData$dose <- dose_info$dose_nM * 1e-9
#' dia_data$ProteinLevelData$drug <- dose_info$drug
#'
#' # Prepare data for analysis
#' prepared_data <- MSstatsPrepareDoseResponseFit(
#'   dia_data$ProteinLevelData,
#'   dose_column = "dose",
#'   drug_column = "drug",
#'   protein_column = "Protein",
#'   log_abundance_column = "LogIntensities"
#' )
#'
#' # Subset for quick example
#' example_data <- prepared_data[prepared_data$protein %in%
#'                               unique(prepared_data$protein)[1:5], ]
#'
#' # Example 1: Quick parallel IC50 without bootstrap (2 cores)
#' ic50_quick_parallel <- PredictIC50Parallel(
#'   data = example_data,
#'   bootstrap = FALSE,
#'   numberOfCores = 2
#' )
#' print(ic50_quick_parallel)
#'
#' @export
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#' @import dplyr
PredictIC50Parallel = function(data,
                               n_samples = 1000,
                               alpha = 0.10,
                               increasing = FALSE,
                               transform_dose = TRUE,
                               ratio_response = TRUE,
                               bootstrap = TRUE,
                               numberOfCores = 2) {

  protein_list = unique(data$protein)

  # Create cluster for parallel processing
  cl = parallel::makeCluster(numberOfCores)

  # Set up environment for cluster
  function_environment = environment()

  # Export all required objects and functions to cluster workers
  parallel::clusterExport(cl,
                          varlist = c("data", "PredictIC50", "n_samples", "alpha",
                                      "increasing", "transform_dose", "ratio_response", "bootstrap"),
                          envir = function_environment)

  # Log progress
  message(paste0("Number of proteins to process: ", length(protein_list)))

  # Create list of protein-drug combinations to process
  loop_list = data %>%
    dplyr::distinct(drug, protein) %>%
    dplyr::filter(drug != "DMSO")

  # Run parallel processing
  results_list = parallel::parLapply(cl, seq_len(nrow(loop_list)), function(i) {
    temp = loop_list[i, ]

    # Log progress every 50 proteins
    if (i %% 50 == 0) {
      message(paste("Finished processing", i, "protein-drug combinations"))
    }

    # Subset data for this protein-drug combination
    data_subset = data[data$protein == temp[[2]] &
                         data$drug %in% c("DMSO", temp[[1]]), ]

    # Run PredictIC50 on this subset
    PredictIC50(data_subset,
                n_samples = n_samples,
                alpha = alpha,
                increasing = increasing,
                transform_dose = transform_dose,
                ratio_response = ratio_response,
                bootstrap = bootstrap)
  })

  # Clean up cluster
  parallel::stopCluster(cl)

  # Combine results
  final_df = do.call(rbind, results_list)
  return(final_df)
}
