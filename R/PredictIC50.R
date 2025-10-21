#' Predict IC50 (dose where response = target) for each protein and drug
#'
#' @param data A data frame with columns: protein, drug, dose, response.
#' @param n_samples Number of bootstrap samples. Default = 1000.
#' @param alpha Confidence level. Default = 0.10.
#' @param increasing Logical. If TRUE, fit a non-decreasing trend. Default = FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default = TRUE.
#' @param ratio_response Logical. If TRUE, use ratio response; else use log2 scale. Default = TRUE.
#' @param bootstrap Logical. If FALSE, skip confidence interval bootstrap estimation and only return IC50. Default = TRUE.
#'@param BPPARAM A \code{BiocParallelParam} for parallel processing. The recommended usage
#'   is to \emph{register} a backend once (e.g., \code{register(MulticoreParam(workers=4))} on
#'   Linux/macOS or \code{register(SnowParam(workers=4, type="SOCK"))} on Windows) and pass
#'   \code{BPPARAM = bpparam()}. Default \code{bpparam()}.
#' @param target_response Numeric, the response fraction (e.g., 0.5, 0.25, 0.75). Default = 0.5.
#'
#' @return A data frame with columns: protein, drug, IC50, IC50_lower_bound, IC50_upper_bound.
#'
#' @examples
#' # Load example data
#' data_path <- system.file("extdata", "DIA_MSstats_Normalized.RDS",
#'                          package = "MSstatsResponse")
#' dia_data <- readRDS(data_path)
#'
#' # Convert GROUP to dose
#' dose_info <- convertGroupToNumericDose(dia_data$ProteinLevelData$GROUP)
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
#' # Subset to fewer proteins for example
#' example_data <- prepared_data[prepared_data$protein %in%
#'                               unique(prepared_data$protein)[1:3], ]
#'
#'
#' # Example 1: Quick IC50 estimation without bootstrap (fast)
#' ic50_quick <- predictIC50(
#'   data = example_data,
#'   bootstrap = FALSE
#' )
#' print(ic50_quick)
#'
#' \dontrun{
#' # Example 2: Full IC50 estimation with bootstrap confidence intervals
#' ic50_results <- predictIC50(
#'   data = prepared_data,
#'   n_samples = 1000,
#'   alpha = 0.10,
#'   ratio_response = TRUE,
#'   bootstrap = TRUE
#' )
#'
#' # Example 3: Parallel processing for large datasets
#' library(BiocParallel)
#' ic50_parallel <- predictIC50(
#'   data = prepared_data,
#'   n_samples = 1000,
#'   BPPARAM = bpparam(),
#'   bootstrap = TRUE
#' )
#'
#' # Example 4: IC50 at different response levels (IC25, IC75)
#' ic25_results <- predictIC50(
#'   data = prepared_data,
#'   target_response = 0.25,
#'   bootstrap = TRUE
#' )
#' }
#'
#' @export
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom data.table rbindlist
#' @importFrom dplyr filter distinct
predictIC50 = function(data,
                       n_samples = 1000,
                       alpha = 0.10,
                       increasing = FALSE,
                       transform_dose = TRUE,
                       ratio_response = TRUE,
                       bootstrap = TRUE,
                       BPPARAM = bpparam(),
                       target_response = 0.5) {

  # Create list of protein-drug combinations to process
  loop_list = data %>%
    distinct(drug, protein) %>%
    filter(drug != "DMSO")

  # Process each combination
  test_results = BiocParallel::bplapply(seq_len(nrow(loop_list)), function(i) {
    temp = loop_list[i, ]
    data_subset = data %>%
      dplyr::filter(drug %in% c("DMSO", temp[[1]]) & protein == temp[[2]])

    .calcSingleIC50(
      df = data_subset,
      n_samples = n_samples,
      alpha = alpha,
      increasing = increasing,
      transform_dose = transform_dose,
      ratio_response = ratio_response,
      bootstrap = bootstrap,
      prot = temp[[2]],
      drug_type = temp[[1]],
      target_response = target_response
    )
  }, BPPARAM = BPPARAM)

  # Combine results
  return(rbindlist(test_results))
}
