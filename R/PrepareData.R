#' Convert MSstats GROUP labels to numeric dose in nM and extract drug name
#'
#' @param group_vector A character or factor vector with GROUP labels (e.g., "Dasatinib_003uM")
#'
#' @return A data frame with two columns: dose_nM (numeric), and drug (character).
#'
#' @examples
#' # Example 1: Basic conversion with mixed units
#' groups <- c("DMSO", "Dasatinib_001uM", "Dasatinib_010uM",
#'             "Dasatinib_100nM", "Dasatinib_1000nM")
#' dose_info <- convertGroupToNumericDose(groups)
#' print(dose_info)
#'
#' # Example 2: Handle multiple drugs
#' multi_drug_groups <- c("DMSO",
#'                       "Dasatinib_001uM", "Dasatinib_010uM",
#'                       "Imatinib_001uM", "Imatinib_010uM")
#' multi_dose_info <- convertGroupToNumericDose(multi_drug_groups)
#' print(multi_dose_info)
#'
#' # Show unique drugs found
#' print(unique(multi_dose_info$drug))
#'
#' @export
convertGroupToNumericDose = function(group_vector) {
  group_vector = as.character(group_vector)

  dose_numeric = numeric(length(group_vector))
  drug_vector = character(length(group_vector))

  for (i in seq_along(group_vector)) {
    label = group_vector[i]

    if (label == "DMSO") {
      dose_numeric[i] = 0
      drug_vector[i] = "DMSO"
    } else {
      parts = strsplit(label, "_")[[1]]
      drug_vector[i] = parts[1]

      value = as.numeric(gsub("[^0-9.]", "", parts[2]))
      unit = gsub("[0-9.]", "", parts[2])

      dose_numeric[i] = if (unit == "uM") value * 1000 else value
    }
  }

  return(data.frame(
    drug = drug_vector,
    dose_nM = dose_numeric,
    stringsAsFactors = FALSE
  ))
}



#' Prepare data for dose-response fitting with isotonic regression
#'
#' @param data A data.frame (e.g. data$ProteinLevelData from MSstats)
#' @param dose_column Name of the column containing dose values (e.g., "dose")
#' @param drug_column Name of column containing treatment name (e.g. drug name )
#' @param protein_column Name of the column containing protein identifiers (e.g., "Protein")
#' @param log_abundance_column Name of the column with log-transformed abundance values (e.g., "LogIntensities")
#' @param transform_nM_to_M Logical. If TRUE, converts dose values from nanomolar (nM) to molar (M)
#'   by multiplying by 10^-9. Use when dose_column contains nM values but analysis requires M units.
#'   Default is NULL (no transformation applied).
#'
#' @return A standardized data.frame with columns: dose, response, protein
#'
#' @examples
#' # Load example data
#' data_path <- system.file("extdata", "DIA_MSstats_Normalized.RDS",
#'                          package = "MSstatsResponse")
#' dia_data <- readRDS(data_path)
#'
#' # Example 1: Basic data preparation with dose already in M
#' # First add dose column if using GROUP labels
#' dose_info <- convertGroupToNumericDose(dia_data$ProteinLevelData$GROUP)
#' dia_data$ProteinLevelData$dose <- dose_info$dose_nM * 1e-9  # Convert to M
#' dia_data$ProteinLevelData$drug <- dose_info$drug
#'
#' prepared_data <- MSstatsPrepareDoseResponseFit(
#'   data = dia_data$ProteinLevelData,
#'   dose_column = "dose",
#'   drug_column = "drug",
#'   protein_column = "Protein",
#'   log_abundance_column = "LogIntensities"
#' )
#'
#' # Check structure
#' str(prepared_data)
#' head(prepared_data)
#'
#' # Example 2: Convert dose from nM to M during preparation
#' dia_data$ProteinLevelData$dose_nM <- dose_info$dose_nM  # Keep in nM
#'
#' prepared_data_converted <- MSstatsPrepareDoseResponseFit(
#'   data = dia_data$ProteinLevelData,
#'   dose_column = "dose_nM",
#'   drug_column = "drug",
#'   protein_column = "Protein",
#'   log_abundance_column = "LogIntensities",
#'   transform_nM_to_M = TRUE  # Convert nM to M
#' )
#'
#' # Verify conversion
#' print(unique(prepared_data_converted$dose))
#'
#' \dontrun{
#' # Example 3: Working with custom column names
#' custom_data <- data.frame(
#'   ProteinID = rep(c("P1", "P2"), each = 10),
#'   Treatment = rep(c("DMSO", "Drug1"), 10),
#'   Concentration = rep(c(0, 1, 10, 100, 1000), 4),
#'   Log2Abundance = rnorm(20, mean = 20, sd = 1)
#' )
#'
#' prepared_custom <- MSstatsPrepareDoseResponseFit(
#'   data = custom_data,
#'   dose_column = "Concentration",
#'   drug_column = "Treatment",
#'   protein_column = "ProteinID",
#'   log_abundance_column = "Log2Abundance",
#'   transform_nM_to_M = TRUE
#' )
#'}
#' @export
MSstatsPrepareDoseResponseFit = function(data,
                                         dose_column = "dose",
                                         drug_column = "drug",
                                         protein_column = "Protein",
                                         log_abundance_column = "LogIntensities",
                                         transform_nM_to_M = NULL) {
  # Check input
  stopifnot(is.data.frame(data))
  required_cols = c(dose_column,drug_column, protein_column, log_abundance_column)

  # Ensure all required columns are present
  missing_cols = setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # Select and rename relevant columns
  subset_df = data[, c(protein_column, drug_column, dose_column, log_abundance_column)]
  colnames(subset_df) = c("protein", "drug", "dose", "response")


  if (!is.null(transform_nM_to_M) && transform_nM_to_M) {
    subset_df$dose_nM = subset_df$dose
    subset_df$dose = subset_df$dose * 10^-9
  }

  return(subset_df)
}
