#' Example pre-processed DIA-MS dataset
#'
#' This dataset contains normalized protein-level data from a DIA-MS
#' chemoproteomics experiment, pre-processed using MSstats.
#'
#' It is used in the MSstatsResponse vignette to demonstrate data
#' formatting and downstream dose–response analysis.
#'
#' @docType data
#' @name DIA_MSstats_Normalized
#' @keywords datasets
#'
#' @format A data frame with protein-level abundance values and
#' associated MSstats metadata column names.
#'
#' @details
#' The dataset is formatted using the standard MSstats preprocessing
#' workflow. For more information on preprocessing mass spectrometry–based
#' proteomics experiments, see the vignettes for MSstats and/or MSstatsTMT.
#'
#' Below is an example of how such data can be prepared:
#'
#' \preformatted{
#' # Read raw data (example with Spectronaut output)
#' raw_data <- readr::read_tsv("path/to/spectronaut_report.tsv")
#'
#' # Convert to MSstats format
#' msstats_data <- MSstats::SpectronauttoMSstatsFormat(raw_data)
#'
#' # Process data: normalization and protein summarization
#' processed_data <- MSstats::dataProcess(
#'   msstats_data,
#'   normalization = "equalizeMedians",  # or FALSE for no normalization
#'   summaryMethod = "TMP",              # Tukey's median polish
#'   MBimpute = TRUE,                    # Impute missing values
#'   maxQuantileforCensored = 0.999
#' )
#' }
#'
#' @examples
#' data("DIA_MSstats_Normalized")
#' head(DIA_MSstats_Normalized)
"DIA_MSstats_Normalized"
