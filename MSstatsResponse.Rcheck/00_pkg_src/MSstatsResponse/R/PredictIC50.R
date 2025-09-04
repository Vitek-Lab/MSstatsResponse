#' Predict IC50 (dose where response = target) for each protein and drug
#'
#' @param data A data frame with columns: protein, drug, dose, response.
#' @param n_samples Number of bootstrap samples. Default = 1000.
#' @param alpha Confidence level. Default = 0.10.
#' @param increasing Logical. If TRUE, fit a non-decreasing trend. Default = FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default = TRUE.
#' @param ratio_response Logical. If TRUE, use ratio response; else use log2 scale. Default = TRUE.
#' @param bootstrap Logical. If FALSE, skip confidence interval bootstrap estimation and only return IC50. Default = TRUE.
#' @param numberOfCores Number of cores for parallel processing. Default = 1.
#' @param target_response Numeric, the response fraction (e.g., 0.5, 0.25, 0.75). Default = 0.5.
#'
#' @return A data frame with columns: protein, drug, IC50, IC50_lower_bound, IC50_upper_bound.
#' @export
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport
#' @importFrom data.table rbindlist
#' @import dplyr
PredictIC50 = function(data,
                       n_samples = 1000,
                       alpha = 0.10,
                       increasing = FALSE,
                       transform_dose = TRUE,
                       ratio_response = TRUE,
                       bootstrap = TRUE,
                       numberOfCores = 1,
                       target_response = 0.5) {

  if (numberOfCores == 1){
    results = .singleCoreIC50(data,
                              n_samples,
                              alpha,
                              increasing,
                              transform_dose,
                              ratio_response,
                              bootstrap,
                              target_response)
  } else{
    results = .multiCoreIC50(data,
                             n_samples,
                             alpha,
                             increasing,
                             transform_dose,
                             ratio_response,
                             bootstrap,
                             numberOfCores,
                             target_response)
  }
  return(results)
}
