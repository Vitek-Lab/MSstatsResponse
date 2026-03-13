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
#' parse_timepoint(c("0hr", "1hr", "12hrs", "168hrs"))
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

#' Example peptide selector: Top N by total signal
#'
#' @param df Data frame with BaseSequence and Intensity columns (for one protein)
#' @param n Number of top peptides to select
#'
#' @return Filtered data frame
#'
#' @examples
#' # Use as peptide_selector
#' ratios <- calculateTurnoverRatios(
#'   feature_data = quant_data$FeatureLevelData,
#'   peptide_selector = function(df) selectTopNPeptides(df, n = 10)
#' )
#'
#' @export
selectTopNPeptides <- function(df, n = 10) {
  top_peptides <- df %>%
    group_by(BaseSequence) %>%
    summarise(total_signal = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_signal)) %>%
    slice_head(n = n) %>%
    pull(BaseSequence)

  df %>% filter(BaseSequence %in% top_peptides)
}
