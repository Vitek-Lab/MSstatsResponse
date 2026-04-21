#' Predict IC50 (dose where response = target) for each protein and drug
#'
#' @param data A data frame with columns: protein, drug, dose, response.
#' @param n_samples Number of bootstrap samples. Default = 1000.
#' @param alpha Confidence level. Default = 0.10.
#' @param increasing Logical or character. If TRUE, fit a non-decreasing trend. If FALSE, fit non-increasing.
#'   If "both", fits both directions and selects the better fit. Default = FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default = TRUE.
#' @param ratio_response Logical. If TRUE, use ratio response; else use log2 scale. Default = TRUE.
#' @param precalculated_ratios Logical. If TRUE, assumes response column contains pre-calculated ratios.
#'   Default is FALSE.
#' @param bootstrap Logical. If FALSE, skip confidence interval bootstrap estimation and only return IC50. Default = TRUE.
#' @param BPPARAM A \code{BiocParallelParam} for parallel processing. Default \code{bpparam()}.
#' @param target_response Numeric, the response fraction (e.g., 0.5, 0.25, 0.75, 1.5).
#'   For decreasing curves: use values < 1 (e.g., 0.5 for IC50).
#'   For increasing curves: use values > 1 (e.g., 1.5 for 50% increase).
#'   Default = 0.5 (auto-adjusts to 1.5 for increasing curves).
#'
#' @return A data frame with columns: protein, drug, IC50, IC50_lower_bound, IC50_upper_bound, direction.
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
#' # Example 2: Fit both increasing and decreasing curves
#' ic50_both <- predictIC50(
#'   data = example_data,
#'   increasing = "both",
#'   bootstrap = FALSE
#' )
#' print(ic50_both)
#'
#' \dontrun{
#' # Example 3: Full IC50 estimation with bootstrap confidence intervals
#' ic50_results <- predictIC50(
#'   data = prepared_data,
#'   n_samples = 1000,
#'   alpha = 0.10,
#'   ratio_response = TRUE,
#'   bootstrap = TRUE
#' )
#'
#' # Example 4: Parallel processing for large datasets
#' library(BiocParallel)
#' ic50_parallel <- predictIC50(
#'   data = prepared_data,
#'   n_samples = 1000,
#'   BPPARAM = bpparam(),
#'   bootstrap = TRUE
#' )
#'
#' # Example 5: IC50 at different response levels (IC25, IC75)
#' ic25_results <- predictIC50(
#'   data = prepared_data,
#'   target_response = 0.25,
#'   bootstrap = TRUE
#' )
#'
#' # Example 6: Increasing curves with EC50 at 1.5x baseline
#' ec50_results <- predictIC50(
#'   data = prepared_data,
#'   increasing = TRUE,
#'   target_response = 1.5,
#'   bootstrap = TRUE
#' )
#'
#' # Example 7: Using pre-calculated ratios
#' turnover_data <- prepared_data
#' turnover_data$response <- 2^turnover_data$response /
#'                           mean(2^turnover_data$response[turnover_data$dose == 0])
#'
#' ic50_precalc <- predictIC50(
#'   data = turnover_data,
#'   precalculated_ratios = TRUE,
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
                       precalculated_ratios = FALSE,
                       bootstrap = TRUE,
                       BPPARAM = bpparam(),
                       target_response = NULL,
                       weights = NULL) {

  # Handle target_response
  if (is.null(target_response)) {
    if (precalculated_ratios) {
      stop("When using precalculated_ratios = TRUE, you must explicitly specify target_response.\n",
           "For decreasing responses (e.g., inhibition): use values < 1 (e.g., 0.5 for IC50)\n",
           "For increasing responses (e.g., activation): use values > 1 (e.g., 1.5 for EC50)")
    } else {
      # Set defaults for package-calculated ratios
      if (is.logical(increasing) && increasing == TRUE) {
        target_response = 1.5
        message("Using target_response = 1.5 for increasing curves (EC50 at 50% increase)")
      } else {
        target_response = 0.5
        message("Using target_response = 0.5 for decreasing curves (IC50 at 50% reduction)")
      }
    }
  }

  # Create list of protein-drug combinations to process
  loop_list = data %>%
    distinct(drug, protein) %>%
    filter(drug != "DMSO")

  # Process each combination
  test_results = BiocParallel::bplapply(seq_len(nrow(loop_list)), function(i) {
    temp = loop_list[i, ]
    row_idx = which(data$drug %in% c("DMSO", temp[[1]]) & data$protein == temp[[2]])
    data_subset = data[row_idx, ]
    w_subset = if (!is.null(weights)) weights[row_idx] else NULL

    .calcSingleIC50(
      df = data_subset,
      n_samples = n_samples,
      alpha = alpha,
      increasing = increasing,
      transform_dose = transform_dose,
      ratio_response = ratio_response,
      precalculated_ratios = precalculated_ratios,
      bootstrap = bootstrap,
      prot = temp[[2]],
      drug_type = temp[[1]],
      target_response = target_response,
      weights = w_subset
    )
  }, BPPARAM = BPPARAM)

  # Combine results
  return(rbindlist(test_results))
}
