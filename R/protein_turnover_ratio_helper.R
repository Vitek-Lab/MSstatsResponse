#' Calculate turnover ratios from MSstats FeatureLevelData
#'
#' @param feature_data Data frame from MSstats dataProcess()$FeatureLevelData
#' @param channel_col Character. Name of column containing Heavy/Light labels. Default = "LABEL"
#' @param heavy_label Character. Value in channel_col indicating heavy channel. Default = "H"
#' @param light_label Character. Value in channel_col indicating light channel. Default = "L"
#' @param time_col Character. Column containing timepoint information. Default = "GROUP"
#' @param peptide_col Character. Column containing peptide sequences. Default = "PEPTIDE"
#' @param protein_col Character. Column containing protein identifiers. Default = "PROTEIN"
#' @param intensity_col Character. Column with intensity values. Default = "INTENSITY"
#' @param run_col Character. Column identifying technical replicates. Default = "RUN"
#' @param peptide_selector Optional function to select peptides.
#'   Function should take a data frame (grouped by protein) and return filtered data frame.
#'   Default = NULL (use all peptides).
#' @param agg_function Function to aggregate duplicate peptide measurements (multiple charge states,
#'   transitions, etc.). Default = max (takes highest signal).
#' @param normalize_tracer Logical. If TRUE, normalize by tracer incorporation. Default = FALSE
#' @param tracer_constants Named numeric vector. Tracer constants for each timepoint.
#'   Required if normalize_tracer = TRUE
#'
#' @return Data frame with columns: Protein, BaseSequence, TimeVal, Run, Heavy, Light, Total, H_frac, L_frac
#'
#' @examples
#' \dontrun{
#' # Basic usage - all proteins, all peptides
#' ratios <- calculateTurnoverRatios(
#'   feature_data = quant_data$FeatureLevelData
#' )
#'
#' # With peptide selection (top 10 by signal per protein)
#' top10_selector <- function(df) {
#'   top_peptides <- df %>%
#'     group_by(BaseSequence) %>%
#'     summarise(total_signal = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
#'     arrange(desc(total_signal)) %>%
#'     slice_head(n = 10) %>%
#'     pull(BaseSequence)
#'
#'   df %>% filter(BaseSequence %in% top_peptides)
#' }
#'
#' ratios_top10 <- calculateTurnoverRatios(
#'   feature_data = quant_data$FeatureLevelData,
#'   peptide_selector = top10_selector
#' )
#'
#' # Use different aggregation function (e.g., median for robustness)
#' ratios_median <- calculateTurnoverRatios(
#'   feature_data = quant_data$FeatureLevelData,
#'   agg_function = median
#' )
#'
#' # With tracer normalization
#' tracer_consts <- c("0hr" = 1.0, "1hr" = 0.95, "12hrs" = 0.85, "168hrs" = 0.75)
#' ratios_norm <- calculateTurnoverRatios(
#'   feature_data = quant_data$FeatureLevelData,
#'   normalize_tracer = TRUE,
#'   tracer_constants = tracer_consts
#' )
#'
#' # Filter for specific protein after calculation
#' protein_ratios <- ratios %>% filter(Protein == "A0A2K5TXF6")
#' }
#'
#' @export
#' @importFrom dplyr filter mutate group_by summarise arrange slice_head pull ungroup select rename
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_remove_all str_extract str_detect
calculateTurnoverRatios <- function(
    feature_data,
    channel_col = "LABEL",
    heavy_label = "H",
    light_label = "L",
    time_col = "GROUP",
    peptide_col = "PEPTIDE",
    protein_col = "PROTEIN",
    intensity_col = "INTENSITY",
    run_col = "RUN",
    peptide_selector = NULL,
    agg_function = max,
    normalize_tracer = FALSE,
    tracer_constants = NULL
) {

  # Check required columns
  required_cols <- c(protein_col, channel_col, time_col, peptide_col, intensity_col, run_col)
  missing_cols <- setdiff(required_cols, colnames(feature_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Process all proteins
  df <- feature_data %>%
    mutate(
      Protein = .data[[protein_col]],
      BaseSequence = str_remove_all(.data[[peptide_col]], "\\[.*?\\]"),
      Label = .data[[channel_col]],
      TimeVal = parse_timepoint(.data[[time_col]]),
      Intensity = .data[[intensity_col]],
      Run = .data[[run_col]]
    ) %>%
    filter(!is.na(TimeVal)) %>%
    filter(Label %in% c(heavy_label, light_label))

  if (nrow(df) == 0) {
    warning("No data found matching Heavy/Light labels")
    return(data.frame())
  }

  # Apply optional peptide selector per protein
  if (!is.null(peptide_selector)) {
    df <- df %>%
      group_by(Protein) %>%
      group_modify(~ peptide_selector(.x)) %>%
      ungroup()
  }

  # Aggregate duplicates (multiple features/transitions for same peptide)
  # This can happen when there are multiple charge states or modifications
  df <- df %>%
    group_by(Protein, BaseSequence, TimeVal, Run, Label) %>%
    summarise(Intensity = agg_function(Intensity, na.rm = TRUE), .groups = "drop")

  # Pivot to wide format (keep replicates separate)
  df_wide <- df %>%
    select(Protein, BaseSequence, TimeVal, Run, Label, Intensity) %>%
    pivot_wider(names_from = Label, values_from = Intensity)

  # Rename heavy/light columns if they're not "Heavy" and "Light"
  if (heavy_label != "Heavy") {
    df_wide <- df_wide %>% rename(Heavy = all_of(heavy_label))
  }
  if (light_label != "Light") {
    df_wide <- df_wide %>% rename(Light = all_of(light_label))
  }

  df_wide <- df_wide %>%
    filter(!is.na(Heavy) & !is.na(Light)) %>%
    mutate(
      Total = Heavy + Light,
      H_frac = Heavy / Total,
      L_frac = Light / Total
    )

  # Optional tracer normalization
  if (normalize_tracer) {
    if (is.null(tracer_constants)) {
      stop("normalize_tracer = TRUE but tracer_constants was not provided")
    }
    
    names(tracer_constants) = as.character(parse_timepoint(names(tracer_constants)))

    df_wide <- df_wide %>%
      mutate(
        tracer_factor = tracer_constants[as.character(TimeVal)],
        H_frac = H_frac / tracer_factor
      )
  }

  df_wide
}

#' Parse timepoint strings to numeric hours
#'
#' @param time_strings Character vector of timepoint labels
#'
#' @return Numeric vector of hours
#'
#' @examples
#' MSstatsResponse:::parse_timepoint(c("0hr", "1hr", "12hrs", "168hrs"))
#' # Returns: 0 1 12 168
#'
#' @keywords internal
parse_timepoint <- function(time_strings) {
  # Extract numeric part
  numeric_part <- as.numeric(str_extract(time_strings, "^[0-9]+"))

  # Handle different units
  is_days <- str_detect(time_strings, "d|day")
  is_weeks <- str_detect(time_strings, "w|week")

  hours <- numeric_part
  hours[is_days] <- numeric_part[is_days] * 24
  hours[is_weeks] <- numeric_part[is_weeks] * 24 * 7

  return(hours)
}


#' Calculate quality-based weights for peptide measurements
#'
#' Calculates weights based on coverage, signal intensity, monotonicity, and data validity.
#' Designed for protein turnover data but applicable to any dose/time-response data.
#'
#' @param data Data frame with peptide-level measurements (output from calculateTurnoverRatios)
#' @param protein_col Character. Column containing protein identifiers. Default = "Protein"
#' @param peptide_col Character. Column containing peptide identifiers. Default = "BaseSequence"
#' @param time_col Character. Column containing timepoint values. Default = "TimeVal"
#' @param response_col Character. Column containing response values (e.g., H_frac). Default = "H_frac"
#' @param light_intensity_col Character. Column containing light channel intensity. Default = "Light"
#' @param validity_threshold Numeric. Maximum allowed response value. Default = 1.3
#' @param top_n_peptides Numeric. Top n peptides based on median light channel intensity. If NULL, no filtering (all get weight 1).
#' @return Input data frame with added columns:
#'   - n_obs: Number of observations for this peptide
#'   - coverage_score: Proportion of timepoints observed
#'   - light_intensity_score: Normalized median light intensity (per protein)
#'   - monotonicity_score: Kendall correlation (time vs response), 0 if decreasing
#'   - validity_flag: 0 if any invalid values, 1 otherwise
#'   - weight: Combined quality weight (product of all components)
#'
#' @examples
#' \dontrun{
#' # Calculate ratios first
#' ratios <- calculateTurnoverRatios(feature_data)
#'
#' # Add quality weights
#' ratios_weighted <- calculatePeptideWeights(ratios)
#'
#' # Inspect weights
#' ratios_weighted %>%
#'   group_by(Protein, BaseSequence) %>%
#'   slice(1) %>%
#'   select(Protein, BaseSequence, coverage_score, monotonicity_score, weight)
#'
#' # Use with doseResponseFit
#' result <- doseResponseFit(
#'   data = ratios_weighted,
#'   weights = ratios_weighted$weight,
#'   increasing = TRUE
#' )
#'
#' # Use stricter validity threshold
#' ratios_weighted_strict <- calculatePeptideWeights(ratios, validity_threshold = 1.0)
#' }
#'
#' @export
#' @importFrom dplyr group_by mutate ungroup across all_of
#' @importFrom stats cor
calculatePeptideWeights <- function(
    data,
    protein_col = "Protein",
    peptide_col = "BaseSequence",
    time_col = "TimeVal",
    response_col = "H_frac",
    light_intensity_col = "Light",
    validity_threshold = 1.3,
    top_n_peptides = NULL
) {
  # Get total number of unique timepoints in the dataset
  n_total_timepoints <- length(unique(data[[time_col]]))

  data_weighted <- data %>%
    # Calculate peptide-level metrics
    group_by(across(all_of(c(protein_col, peptide_col)))) %>%
    arrange(across(all_of(time_col))) %>%
    mutate(
      # Coverage: proportion of timepoints observed (peptide-level)
      n_obs = n(),
      coverage_score = n_obs / n_total_timepoints,
      # Light intensity: median for this peptide (peptide-level)
      peptide_median_light = median(.data[[light_intensity_col]], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Light intensity filter within each protein
    group_by(across(all_of(protein_col))) %>%
    mutate(
      # Create peptide-level rank (same rank for all observations of a peptide)
      peptide_rank = dense_rank(-peptide_median_light)
    ) %>%
    ungroup() %>%
    mutate(
      # Binary top-N filter: 1 if in top N, 0 otherwise
      light_intensity_score = if (is.null(top_n_peptides)) {
        1  # No filtering applied
      } else {
        if_else(peptide_rank <= top_n_peptides, 1, 0)
      }
    ) %>%
    # Back to peptide level for monotonicity
    group_by(across(all_of(c(protein_col, peptide_col)))) %>%
    mutate(
      # Monotonicity: Kendall correlation between time and response (peptide-level)
      monotonicity_score = pmax(0,
                                cor(.data[[time_col]], .data[[response_col]],
                                    method = "kendall", use = "complete.obs")
      ),
      # Handle NA monotonicity (peptides with <2 observations)
      monotonicity_score = if_else(is.na(monotonicity_score), 0, monotonicity_score)
    ) %>%
    ungroup() %>%
    # Validity check at observation level (row-by-row)
    mutate(
      validity_flag = if_else(
        .data[[response_col]] > validity_threshold |
          .data[[response_col]] < 0 |
          is.na(.data[[response_col]]),
        0,
        1
      ),
      # Combined weight (multiplicative)
      weight = coverage_score * light_intensity_score * monotonicity_score * validity_flag
    )

  return(data_weighted)
}
