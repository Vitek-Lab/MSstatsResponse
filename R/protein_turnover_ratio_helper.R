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
        H_frac = H_frac / tracer_factor,
        L_frac =  pmax(1 - H_frac, 0)
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
#' Coverage is scored via a binomial CDF: for each peptide, P(X <= k | n, p) where k is
#' the number of non-zero timepoints detected, n is the total non-zero timepoints in the
#' experiment, and p is the protein-level mean detection rate across all its peptides.
#' This penalizes peptides with unusually low coverage relative to the protein's norm.
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
#'   - n_obs: Total number of observations for this peptide
#'   - k_obs: Number of non-zero timepoints where this peptide is detected
#'   - coverage_per_peptide: k_obs / n (per-peptide detection proportion)
#'   - p_protein: Protein-level mean detection rate across all its peptides
#'   - coverage_score: P(X <= k_obs | n, p_protein) — binomial CDF coverage score
#'   - light_intensity_score: 1 (no filter) or binary top-N indicator (per protein)
#'   - monotonicity_score: Kendall correlation (time vs response), floored at 0
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
#' # Inspect coverage diagnostics
#' ratios_weighted %>%
#'   group_by(Protein, BaseSequence) %>%
#'   slice(1) %>%
#'   select(Protein, BaseSequence, k_obs, p_protein, coverage_score, monotonicity_score, weight)
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
#' @importFrom dplyr group_by mutate ungroup across all_of distinct summarise left_join if_else dense_rank
#' @importFrom stats cor pbinom median
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
  # Total non-zero timepoints: t=0 excluded because it precedes label incorporation
  all_tps <- unique(data[[time_col]])
  non_zero_tps <- all_tps[all_tps != 0]
  n_coverage <- if (length(non_zero_tps) > 0) length(non_zero_tps) else length(all_tps)

  data_weighted <- data %>%
    group_by(across(all_of(c(protein_col, peptide_col)))) %>%
    arrange(across(all_of(time_col))) %>%
    mutate(
      n_obs = n(),
      k_obs = sum(.data[[time_col]] != 0),
      coverage_per_peptide = k_obs / n_coverage,
      peptide_median_light = median(.data[[light_intensity_col]], na.rm = TRUE)
    ) %>%
    ungroup()

  # Protein-level detection rate: mean coverage proportion across unique peptides
  protein_p <- data_weighted %>%
    distinct(across(all_of(c(protein_col, peptide_col))), .keep_all = TRUE) %>%
    group_by(across(all_of(protein_col))) %>%
    summarise(p_protein = mean(coverage_per_peptide), .groups = "drop")

  data_weighted <- data_weighted %>%
    left_join(protein_p, by = protein_col) %>%
    mutate(
      # Binomial CDF: low score = unusually low coverage relative to protein norm
      coverage_score = pbinom(k_obs, n_coverage, p_protein)
    ) %>%
    group_by(across(all_of(protein_col))) %>%
    mutate(
      peptide_rank = dense_rank(-peptide_median_light)
    ) %>%
    ungroup() %>%
    mutate(
      light_intensity_score = if (is.null(top_n_peptides)) {
        1
      } else {
        if_else(peptide_rank <= top_n_peptides, 1, 0)
      }
    ) %>%
    group_by(across(all_of(c(protein_col, peptide_col)))) %>%
    mutate(
      monotonicity_score = pmax(0,
                                cor(.data[[time_col]], .data[[response_col]],
                                    method = "kendall", use = "complete.obs")
      ),
      monotonicity_score = if_else(is.na(monotonicity_score), 0, monotonicity_score)
    ) %>%
    ungroup() %>%
    mutate(
      validity_flag = if_else(
        .data[[response_col]] > validity_threshold |
          .data[[response_col]] < 0 |
          is.na(.data[[response_col]]),
        0,
        1
      ),
      weight = coverage_score * light_intensity_score * monotonicity_score * validity_flag
    )

  return(data_weighted)
}


#' Calculate per-protein QC score from light-channel coverage
#'
#' Computes the fraction of expected (peptide x timepoint) light-channel
#' observations actually detected for each protein. This is a pure
#' measurement-quality metric -- it ignores the heavy channel entirely, so it
#' remains meaningful even for proteins where no heavy incorporation occurred
#' (e.g., very long-lived proteins, or peptides without heavy-label residues).
#'
#' For each protein with n distinct light peptides observed across T total
#' timepoints in the experiment, the maximum possible (peptide x timepoint)
#' light measurements is n * T. The qc_score is the fraction of those that
#' were actually observed, capped at 1.
#'
#' @param feature_data Data frame from MSstats dataProcess()$FeatureLevelData.
#' @param protein_col Character. Column containing protein identifiers. Default = "PROTEIN"
#' @param peptide_col Character. Column containing peptide sequences. Default = "PEPTIDE"
#' @param label_col Character. Column containing Heavy/Light labels. Default = "LABEL"
#' @param intensity_col Character. Column with intensity values. Default = "INTENSITY"
#' @param time_col Character. Column containing timepoint information. Default = "GROUP"
#' @param light_label Character. Value in label_col indicating light channel. Default = "L"
#'
#' @return Data frame with one row per protein containing:
#'   - protein identifier (column name per `protein_col`)
#'   - n_light_peptides: distinct peptides observed in light channel
#'   - observed_cells: distinct (peptide x timepoint) cells observed
#'   - n_max_possible: n_light_peptides x n_distinct_timepoints
#'   - qc_score: observed_cells / n_max_possible, capped at 1
#'
#' @examples
#' \dontrun{
#' qc <- calculateQCScore(quant_data$FeatureLevelData)
#' qc %>% arrange(desc(qc_score)) %>% head()
#' }
#'
#' @export
#' @importFrom dplyr filter group_by summarise n_distinct
calculateQCScore <- function(feature_data,
                             protein_col   = "PROTEIN",
                             peptide_col   = "PEPTIDE",
                             label_col     = "LABEL",
                             intensity_col = "INTENSITY",
                             time_col      = "GROUP",
                             light_label   = "L") {
  n_tp <- n_distinct(feature_data[[time_col]])

  feature_data %>%
    filter(.data[[label_col]] == light_label, is.finite(.data[[intensity_col]])) %>%
    group_by(.data[[protein_col]]) %>%
    summarise(
      n_light_peptides = n_distinct(.data[[peptide_col]]),
      observed_cells   = n_distinct(paste(.data[[peptide_col]], .data[[time_col]])),
      n_max_possible   = n_light_peptides * n_tp,
      qc_score         = pmin(1, observed_cells / n_max_possible),
      .groups = "drop"
    )
}


#' Calculate per-protein confidence score for turnover fits
#'
#' Combines peptide-quality weights, fit residuals, light-channel QC, and
#' a Bayesian shrinkage factor on heavy-peptide count into a single
#' per-protein confidence score in \[0, 1\]. Higher values indicate the
#' dose-response fit is supported by clean, complete, abundant data.
#'
#' Formula:
#'
#'   confidence = mean_weight * 1/(1 + SSE_Full) * qc_score * n_heavy / (n_heavy + k_shrinkage)
#'
#' The `n_heavy / (n_heavy + k)` term is a Bayesian / Laplace smoothing
#' factor that penalizes thin peptide support: with `k_shrinkage = 2`, a
#' 1-peptide protein is capped at 1/3 of its otherwise-achievable score.
#'
#' @param weights_df Output of `calculatePeptideWeights()`.
#' @param fit_df Output of `doseResponseFit()` (must contain an SSE column).
#' @param qc_df Output of `calculateQCScore()`.
#' @param feature_data Raw feature-level data (used to count heavy peptides per protein).
#' @param protein_col Character. Column in weights_df / fit_df identifying proteins. Default = "Protein"
#' @param weight_col Character. Column in weights_df containing per-observation weight. Default = "weight"
#' @param sse_col Character. Column in fit_df containing SSE. Default = "SSE_Full"
#' @param qc_protein_col Character. Protein column in qc_df. Default = "PROTEIN"
#' @param qc_score_col Character. QC score column in qc_df. Default = "qc_score"
#' @param feature_protein_col Character. Protein column in feature_data. Default = "PROTEIN"
#' @param feature_peptide_col Character. Peptide column in feature_data. Default = "PEPTIDE"
#' @param feature_label_col Character. Label column in feature_data. Default = "LABEL"
#' @param heavy_label Character. Value indicating heavy channel. Default = "H"
#' @param k_shrinkage Numeric. Bayesian shrinkage constant for the peptide-count factor.
#'   Larger values penalize low-peptide proteins more strongly. Default = 2
#'
#' @return The input fit_df with additional columns:
#'   - mean_weight: average per-observation weight
#'   - n_obs: number of observations used in the fit
#'   - qc_score: joined from qc_df
#'   - n_heavy_peptides: count from feature_data
#'   - pep_factor: n_heavy / (n_heavy + k_shrinkage)
#'   - confidence: combined score in \[0, 1\]
#'
#' @examples
#' \dontrun{
#' qc      <- calculateQCScore(df_feat)
#' weights <- calculatePeptideWeights(ratios)
#' fit     <- doseResponseFit(ratios, increasing = TRUE, precalculated_ratios = TRUE)
#' conf    <- calculateConfidence(weights, fit, qc, df_feat)
#'
#' conf %>% arrange(desc(confidence)) %>% head()
#' }
#'
#' @export
#' @importFrom dplyr filter group_by summarise n_distinct left_join mutate select rename all_of coalesce
calculateConfidence <- function(weights_df,
                                fit_df,
                                qc_df,
                                feature_data,
                                protein_col           = "Protein",
                                weight_col            = "weight",
                                sse_col               = "SSE_Full",
                                qc_protein_col        = "PROTEIN",
                                qc_score_col          = "qc_score",
                                feature_protein_col   = "PROTEIN",
                                feature_peptide_col   = "PEPTIDE",
                                feature_label_col     = "LABEL",
                                heavy_label           = "H",
                                k_shrinkage           = 2) {

  pep_counts <- feature_data %>%
    filter(.data[[feature_label_col]] == heavy_label) %>%
    group_by(.data[[feature_protein_col]]) %>%
    summarise(n_heavy_peptides = n_distinct(.data[[feature_peptide_col]]),
              .groups = "drop") %>%
    rename(!!protein_col := !!feature_protein_col)

  mean_weights <- weights_df %>%
    group_by(.data[[protein_col]]) %>%
    summarise(mean_weight = mean(.data[[weight_col]], na.rm = TRUE),
              n_obs       = n(),
              .groups     = "drop")

  qc_simple <- qc_df %>%
    select(all_of(c(qc_protein_col, qc_score_col))) %>%
    rename(!!protein_col := !!qc_protein_col)

  fit_df %>%
    left_join(mean_weights, by = protein_col) %>%
    left_join(qc_simple,    by = protein_col) %>%
    left_join(pep_counts,   by = protein_col) %>%
    mutate(
      n_heavy_peptides = coalesce(n_heavy_peptides, 0L),
      pep_factor       = n_heavy_peptides / (n_heavy_peptides + k_shrinkage),
      confidence       = mean_weight * (1 / (1 + .data[[sse_col]])) *
                           .data[[qc_score_col]] * pep_factor,
      confidence       = pmin(1, pmax(0, confidence))
    )
}


#' Classify proteins into turnover categories and confidence tiers
#'
#' Assigns each protein a biological `category` describing its turnover
#' behavior and a `tier` (HIGH / MEDIUM / LOW) summarizing scoring confidence.
#' Combines QC + confidence + IC50 predictions in one call.
#'
#' Categories:
#'   - `fit`: IC50 reached at the long-target response (default 0.50)
#'   - `medium_lived`: reached the short target (0.21) but not the long target
#'   - `long_lived`: failed both targets and max H_frac stayed below 0.5
#'   - `fast`: failed both targets but max H_frac exceeds 0.5 (IC50 below observed range)
#'   - `no_heavy`: no fit was possible (no paired heavy peptides)
#'
#' Tiers use percentile cutoffs computed from the input data. Proteins with a
#' fit are tiered on `confidence`; `no_heavy` proteins are tiered on `qc_score`.
#' HIGH tier additionally requires a minimum number of observations.
#'
#' @param weights_df Output of `calculatePeptideWeights()`.
#' @param fit_df Output of `doseResponseFit()`.
#' @param qc_df Output of `calculateQCScore()`.
#' @param conf_df Output of `calculateConfidence()`.
#' @param high_quantile Numeric. Upper percentile cutoff for HIGH tier. Default = 0.85 (top 15%).
#' @param low_quantile Numeric. Lower percentile cutoff for LOW tier. Default = 0.25 (bottom 25%).
#' @param min_obs Numeric. Minimum observations required for HIGH tier. Default = 3.
#' @param target_short Numeric. Lower response target passed to predictIC50() for lifetime
#'   classification. Default = 0.21.
#' @param target_long Numeric. Upper response target passed to predictIC50(). Default = 0.50.
#'
#' @return Data frame with one row per protein containing input columns plus:
#'   - max_h_frac: per-protein maximum H_frac
#'   - category: one of `fit`, `medium_lived`, `long_lived`, `fast`, `no_heavy`
#'   - tier: one of `HIGH`, `MEDIUM`, `LOW`
#'
#' @examples
#' \dontrun{
#' qc      <- calculateQCScore(df_feat)
#' weights <- calculatePeptideWeights(ratios)
#' fit     <- doseResponseFit(ratios, increasing = TRUE, precalculated_ratios = TRUE)
#' conf    <- calculateConfidence(weights, fit, qc, df_feat, k_shrinkage = 2)
#'
#' final   <- classifyTurnoverProteins(weights, fit, qc, conf)
#' final %>% count(category, tier)
#' }
#'
#' @export
#' @importFrom dplyr filter pull group_by summarise rename left_join select any_of mutate case_when
#' @importFrom stats quantile
classifyTurnoverProteins <- function(weights_df,
                                     fit_df,
                                     qc_df,
                                     conf_df,
                                     high_quantile = 0.85,
                                     low_quantile  = 0.25,
                                     min_obs       = 3,
                                     target_short  = 0.21,
                                     target_long   = 0.50) {

  predict_at <- function(target) {
    predictIC50(weights_df,
                increasing = TRUE, transform_dose = FALSE,
                precalculated_ratios = TRUE, bootstrap = FALSE,
                target_response = target)
  }
  short_NA <- predict_at(target_short) %>% filter(is.na(IC50)) %>% pull(Protein)
  long_NA  <- predict_at(target_long)  %>% filter(is.na(IC50)) %>% pull(Protein)
  long_lived <- intersect(short_NA, long_NA)

  max_hfrac <- weights_df %>%
    group_by(Protein) %>%
    summarise(max_h_frac = max(H_frac, na.rm = TRUE), .groups = "drop")

  out <- qc_df %>%
    rename(Protein = PROTEIN) %>%
    left_join(conf_df %>% select(-any_of("qc_score")), by = "Protein") %>%
    left_join(max_hfrac, by = "Protein") %>%
    mutate(
      category = case_when(
        is.na(confidence)                                       ~ "no_heavy",
        Protein %in% long_lived & max_h_frac > 0.5              ~ "fast",
        Protein %in% long_lived & max_h_frac < 0.35           ~ "long_lived",
        Protein %in% long_NA & !(Protein %in% short_NA)         ~ "medium_lived",
        TRUE                                                     ~ "fit"
      )
    )

  conf_high <- quantile(out$confidence, high_quantile, na.rm = TRUE)
  conf_low  <- quantile(out$confidence, low_quantile,  na.rm = TRUE)
  qc_high   <- quantile(out$qc_score[out$category == "no_heavy"], high_quantile, na.rm = TRUE)
  qc_low    <- quantile(out$qc_score[out$category == "no_heavy"], low_quantile,  na.rm = TRUE)

  out %>%
    mutate(
      tier = case_when(
        category == "no_heavy" & qc_score >= qc_high & observed_cells >= min_obs ~ "HIGH",
        category == "no_heavy" & qc_score >= qc_low                              ~ "MEDIUM",
        category == "no_heavy"                                                    ~ "LOW",
        confidence >= conf_high & n_obs >= min_obs                                ~ "HIGH",
        confidence >= conf_low                                                     ~ "MEDIUM",
        TRUE                                                                       ~ "LOW"
      )
    )
}
