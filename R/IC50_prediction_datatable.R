#' Parallel version of PredictIC50 function
#'
#' Runs PredictIC50 on the entire data set in parallel across proteins.
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
#' @export
PredictIC50Parallel = function(data,
                               n_samples = 1000,
                               alpha = 0.10,
                               increasing = FALSE,
                               transform_dose = TRUE,
                               ratio_response = TRUE,
                               bootstrap = TRUE,
                               numberOfCores = 2) {

  library(parallel)

  protein_list = unique(data$protein)
  cl = makeCluster(numberOfCores)
  function_environment = environment()

  # Export all required objects and functions
  clusterExport(cl,
                varlist = c("data", "PredictIC50", "n_samples", "alpha",
                            "increasing", "transform_dose", "ratio_response", "bootstrap"),
                envir = function_environment)

  cat(paste0("Number of proteins to process: ", length(protein_list)),
      sep = "\n", file = "PredictIC50_parallel_log_progress.log")

  results_list = parLapply(cl, seq_along(protein_list), function(i) {
    prot = protein_list[i]
    if (i %% 50 == 0) {
      cat("Finished processing 50 additional proteins\n",
          file = "PredictIC50_parallel_log_progress.log", append = TRUE)
    }

    # Subset data for this protein
    data_subset = data[data$protein == prot, ]

    # Run PredictIC50 on this protein's data
    PredictIC50(data_subset,
                n_samples = n_samples,
                alpha = alpha,
                increasing = increasing,
                transform_dose = transform_dose,
                ratio_response = ratio_response)
               # bootstrap = bootstrap
  })

  stopCluster(cl)

  final_df = do.call(rbind, results_list)
  return(final_df)
}
