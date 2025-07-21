#' Convert MSstats GROUP labels to numeric dose in nM and extract drug name
#'
#' @param group_vector A character or factor vector with GROUP labels (e.g., "Dasatinib_003uM")
#'
#' @return A data frame with two columns: dose_nM (numeric), and drug (character).
#' @export
ConvertGroupToNumericDose = function(group_vector) {
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
#'
#' @return A standardized data.frame with columns: dose, response, protein
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

  if (transform_nM_to_M){
    subset_df$dose_nM = subset_df$dose
    subset_df$dose = subset_df$dose * 10^-9
  }

  return(subset_df)
}
