#' Drug-protein interaction detection tested by F-test (fitted curve vs average response)
#'
#' Fits an isotonic regression model to protein abundance data.
#' Performs an F-test to assess the significance of the dose-response curve and applies FDR correction.
#'
#' @param data Protein-level data, formatted with MSstatsPreparedoseResponseFit().
#' @param weights Optional numeric vector of weights. Defaults to equal weights.
#' @param increasing Logical or character. If TRUE, fits a non-decreasing model. If FALSE, fits non-increasing.
#'   If "both", fits both increasing and decreasing models and selects the one with better fit (lower SSE).
#'   Recommended for exploratory analysis, but ~2x slower. Default is FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default is TRUE.
#' @param ratio_response Logical. If TRUE, converts log2 abundance to ratios relative to DMSO. Default is FALSE.
#' @param precalculated_ratios Logical. If TRUE, assumes the response column contains pre-calculated ratios
#'   (not log2 values) and skips internal ratio calculation. This is useful for protein turnover data
#'   or other datasets where ratios are computed externally. Default is FALSE.
#'
#' @return A data frame with protein-wise F-test results and BH-adjusted p-values.
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
#' # Subset for quick example
#' example_data <- prepared_data[prepared_data$protein %in%
#'                               unique(prepared_data$protein)[1:5], ]
#'
#' # Example 1: Basic interaction detection on log2 scale
#' interaction_results <- doseResponseFit(
#'   data = example_data,
#'   increasing = FALSE,
#'   transform_dose = TRUE,
#'   ratio_response = FALSE
#' )
#'
#' # Example 2: Fit both directions and select better
#' interaction_both <- doseResponseFit(
#'   data = example_data,
#'   increasing = "both",
#'   transform_dose = TRUE,
#'   ratio_response = FALSE
#' )
#'
#' # Check which direction was selected for each protein
#' table(interaction_both$direction)
#'
#' # Example 3: Using pre-calculated ratios (e.g., from protein turnover)
#' # Assume you have pre-calculated ratio data
#' turnover_data <- example_data
#' turnover_data$response <- 2^turnover_data$response /
#'                           mean(2^turnover_data$response[turnover_data$dose == 0])
#'
#' interaction_results_precalc <- doseResponseFit(
#'   data = turnover_data,
#'   increasing = FALSE,
#'   transform_dose = TRUE,
#'   ratio_response = FALSE,
#'   precalculated_ratios = TRUE
#' )
#'
#' \dontrun{
#' # Example 4: Full dataset analysis
#' full_results <- doseResponseFit(
#'   data = prepared_data,
#'   increasing = FALSE,
#'   transform_dose = TRUE,
#'   ratio_response = FALSE
#' )
#' }
#'
#' @export
#' @importFrom stats pf p.adjust
#' @importFrom dplyr filter select mutate group_by summarise arrange distinct
#' @importFrom data.table rbindlist
#' @importFrom utils setTxtProgressBar txtProgressBar
doseResponseFit = function(data, weights = NULL,
                           increasing = FALSE,
                           transform_dose = TRUE,
                           ratio_response = FALSE,
                           precalculated_ratios = FALSE) {

  drug_list = unique(data$drug[data$drug != "DMSO"])

  # If no drugs to test, return empty data frame with correct structure
  if (length(drug_list) == 0) {
    return(data.frame(
      Protein = character(),
      drug = character(),
      SSE_Full = numeric(),
      SSE_Null = numeric(),
      F_statistic = numeric(),
      P_value = numeric(),
      adj.pvalue = numeric(),
      direction = character(),
      log2FC = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  all_results = list()

  for (drug_type in drug_list) {
    drug_subset = data %>% dplyr::filter(drug %in% c("DMSO", drug_type))
    protein_list = unique(drug_subset$protein)
    results_list = list()

    progress_bar = txtProgressBar(min = 0, max = length(protein_list), style = 3)
    for (i in seq_along(protein_list)) {
      tryCatch({
        suppressWarnings({
          data_single_protein = drug_subset %>% dplyr::filter(protein == protein_list[i])

          x = data_single_protein$dose
          y = data_single_protein$response

          # Handle weights
          if (is.null(weights)) {
            w = rep(1, length(y))
          } else if (length(weights) == nrow(data)) {
            protein_drug_idx = which(data$protein == protein_list[i] &
                                       data$drug %in% c("DMSO", drug_type))
            w = weights[protein_drug_idx]
          } else if (length(weights) == length(y)) {
            w = weights
          } else {
            warning("Weights length doesn't match data, using equal weights")
            w = rep(1, length(y))
          }

          # Handle "both" option for increasing parameter
          if (is.character(increasing) && increasing == "both") {
            # Fit both directions - keep both results for FDR correction
            fit_dec = fitIsotonicRegression(
              x = x, y = y, w = w,
              increasing = FALSE,
              transform_x = transform_dose,
              ratio_y = ratio_response,
              precalculated_ratios = precalculated_ratios,
              test_significance = TRUE
            )

            fit_inc = fitIsotonicRegression(
              x = x, y = y, w = w,
              increasing = TRUE,
              transform_x = transform_dose,
              ratio_y = ratio_response,
              precalculated_ratios = precalculated_ratios,
              test_significance = TRUE
            )

            # Store both results temporarily for FDR correction
            results_dec = fit_dec$f_test
            results_dec$Protein = protein_list[i]
            results_dec$drug = drug_type
            results_dec$direction = "decreasing"
            results_dec$F_statistic_alt = fit_inc$f_test$F_statistic  # Store alternative F for later comparison

            results_inc = fit_inc$f_test
            results_inc$Protein = protein_list[i]
            results_inc$drug = drug_type
            results_inc$direction = "increasing"
            results_inc$F_statistic_alt = fit_dec$f_test$F_statistic

            # Fitted max fold change: y_pred is on ratio scale if ratio_response=TRUE or precalculated_ratios=TRUE
            on_ratio_scale_fit = ratio_response || precalculated_ratios

            if (on_ratio_scale_fit) {
              log2FC_dec = log2(max(fit_dec$y_pred, na.rm = TRUE) / max(min(fit_dec$y_pred, na.rm = TRUE), 1e-10))
              log2FC_inc = log2(max(fit_inc$y_pred, na.rm = TRUE) / max(min(fit_inc$y_pred, na.rm = TRUE), 1e-10))
            } else {
              log2FC_dec = max(fit_dec$y_pred, na.rm = TRUE) - min(fit_dec$y_pred, na.rm = TRUE)
              log2FC_inc = max(fit_inc$y_pred, na.rm = TRUE) - min(fit_inc$y_pred, na.rm = TRUE)
            }

            results_dec$log2FC = log2FC_dec
            results_inc$log2FC = log2FC_inc

            # Add both to results list
            results_list[[length(results_list) + 1]] = results_dec
            results_list[[length(results_list) + 1]] = results_inc

          } else {
            # Original behavior - single direction
            fit = fitIsotonicRegression(
              x = x, y = y, w = w,
              increasing = increasing,
              transform_x = transform_dose,
              ratio_y = ratio_response,
              precalculated_ratios = precalculated_ratios,
              test_significance = TRUE
            )
            chosen_direction = if(is.logical(increasing) && increasing) "increasing" else "decreasing"

            results_temp = fit$f_test
            results_temp$Protein = protein_list[i]
            results_temp$drug = drug_type
            results_temp$direction = chosen_direction

            on_ratio_scale_fit = ratio_response || precalculated_ratios

            if (on_ratio_scale_fit) {
              results_temp$log2FC = log2(max(fit$y_pred, na.rm = TRUE) / max(min(fit$y_pred, na.rm = TRUE), 1e-10))
            } else {
              results_temp$log2FC = max(fit$y_pred, na.rm = TRUE) - min(fit$y_pred, na.rm = TRUE)
            }

            results_temp = results_temp[, c("Protein", "drug", "direction",
                                            setdiff(names(results_temp), c("Protein", "drug", "direction")))]
            results_list[[i]] = results_temp
          }
          setTxtProgressBar(progress_bar, i)
        })
      }, error = function(e) {
        warning(paste("Error for drug:", drug_type, "protein:", protein_list[i], ":", conditionMessage(e)))
      })
    }
    close(progress_bar)

    # Combine and adjust p-values for this drug
    if (length(results_list) > 0) {
      results_df = do.call(rbind, results_list)

      # If using "both", we need to handle FDR correction and selection
      if (is.character(increasing) && increasing == "both") {
        # Separate FDR correction for increasing and decreasing
        inc_idx = which(results_df$direction == "increasing")
        dec_idx = which(results_df$direction == "decreasing")

        # Initialize adj.pvalue column
        results_df$adj.pvalue = NA

        # FDR correction within each direction
        if (length(inc_idx) > 0) {
          results_df$adj.pvalue[inc_idx] = p.adjust(results_df$P_value[inc_idx], method = "BH")
        }
        if (length(dec_idx) > 0) {
          results_df$adj.pvalue[dec_idx] = p.adjust(results_df$P_value[dec_idx], method = "BH")
        }

        # Now select the better direction for each protein
        proteins_unique = unique(results_df$Protein)
        final_results = list()

        for (prot in proteins_unique) {
          prot_rows = results_df[results_df$Protein == prot, ]

          if (nrow(prot_rows) == 2) {
            # Compare F-statistics to select better fit
            dec_row = prot_rows[prot_rows$direction == "decreasing", ]
            inc_row = prot_rows[prot_rows$direction == "increasing", ]

            if (dec_row$F_statistic >= inc_row$F_statistic) {
              selected_row = dec_row
            } else {
              selected_row = inc_row
            }

            # Remove the F_statistic_alt column
            selected_row$F_statistic_alt = NULL
            final_results[[length(final_results) + 1]] = selected_row
          } else {
            # Only one direction (shouldn't happen, but handle it)
            prot_rows$F_statistic_alt = NULL
            final_results[[length(final_results) + 1]] = prot_rows
          }
        }

        results_df = do.call(rbind, final_results)

      } else {
        # Standard FDR correction when direction is fixed
        results_df$adj.pvalue = p.adjust(results_df$P_value, method = "BH")
      }

      # Reorder columns to put direction after drug
      results_df = results_df[, c("Protein", "drug", "direction",
                                  setdiff(names(results_df), c("Protein", "drug", "direction")))]
      all_results[[drug_type]] = results_df
    }
  }

  # Return empty data frame if no results
  if (length(all_results) == 0) {
    return(data.frame(
      protein = character(),
      drug = character(),
      SSE_Full = numeric(),
      SSE_Null = numeric(),
      F_statistic = numeric(),
      P_value = numeric(),
      adj.pvalue = numeric(),
      direction = character(),
      stringsAsFactors = FALSE
    ))
  }

  final_results = do.call(rbind, all_results)
  rownames(final_results) = NULL
  return(final_results)
}

#' Fit Isotonic Regression Model
#'
#' Fits an isotonic regression model to protein intensity data with log-transformed drug doses.
#' Optionally performs an F-test to assess the significance of the dose-response curve.
#'
#' @param x Numeric vector of dose values.
#' @param y Numeric vector of response values.
#' @param w Optional numeric vector of weights. Defaults to equal weights.
#' @param increasing Logical. If TRUE, fits a non-decreasing model. If FALSE, fits non-increasing.
#' @param transform_x Logical. If TRUE, applies log10(x + 1) transformation. Default is TRUE.
#' @param ratio_y Logical. If TRUE, converts log2 abundance to ratios relative to DMSO. Default is FALSE.
#' @param precalculated_ratios Logical. If TRUE, assumes y contains pre-calculated ratios and skips
#'   internal ratio calculation. Default is FALSE.
#' @param test_significance Logical. If TRUE, performs F-test to assess significance.
#'
#' @return A list representing the isotonic regression model (class = "isotonic_model").
#'
#' @importFrom stats pf approx quantile
fitIsotonicRegression = function(x, y, w = rep(1, length(y)),
                                 increasing = FALSE,
                                 transform_x = TRUE,
                                 ratio_y = FALSE,
                                 precalculated_ratios = FALSE,
                                 test_significance = FALSE) {
  stopifnot(length(x) == length(y), length(y) == length(w))

  x_org = x

  if (transform_x) {
    x = log10(x + 1)
  }

  # Handle ratio calculation
  if (precalculated_ratios) {
    # y already contains ratios - use as-is
    # No transformation needed
  } else if (ratio_y) {
    # Original behavior: convert log2 to ratios
    y = 2^y
    baseline = mean(y[x_org == 0], na.rm = TRUE)
    y = y / baseline
  }
  # Otherwise y stays as-is (log2 values)

  order_idx = order(x)
  x_sorted = x[order_idx]
  y_sorted = y[order_idx]
  w_sorted = w[order_idx]

  if (!increasing) {
    y_sorted = -y_sorted
  }

  fit_y = y_sorted
  n_obs = length(fit_y)
  i = 1
  target = seq_len(n_obs)

  while (i <= n_obs) {
    k = target[i] + 1
    if (k > n_obs) break # stops if at end last value
    if (fit_y[i] < fit_y[k]) { # if breaks constraint , do not let line change
      i = k
      next
    } # if constraint is satisfied perform model fit
    sum_wy = w_sorted[i] * fit_y[i]
    sum_w = w_sorted[i]
    repeat {
      sum_wy = sum_wy + w_sorted[k] * fit_y[k]
      sum_w = sum_w + w_sorted[k]
      k = target[k] + 1
      if (k > n_obs || fit_y[k - 1] < fit_y[k]) {
        fit_y[i] = sum_wy / sum_w
        w_sorted[i] = sum_w
        target[i] = k - 1
        target[k - 1] = i
        if (i > 1) i = target[i - 1]
        break
      }
    }
  }

  i = 1
  while (i <= n_obs) {
    k = target[i] + 1
    if (k > i + 1) {
      fit_y[(i + 1):(k - 1)] = fit_y[i]
    }
    i = k
  }

  if (!increasing) {
    fit_y = -fit_y
  }

  y_final = numeric(n_obs)
  y_final[order_idx] = fit_y

  # fix to duplicate x values warnings from approx function
  if (length(unique(x)) < length(x)) {
    # Aggregate y values for duplicate x values
    x_unique = unique(x)
    y_aggregated = sapply(x_unique, function(xi) {
      mean(y_final[x == xi])
    })
    y_pred_new = stats::approx(x_unique, y_aggregated, xout = x, rule = 2)$y
  } else {
    y_pred_new = stats::approx(x, y_final, xout = x, rule = 2)$y
  }

  result = list(
    x = x,
    y_pred = y_pred_new,
    original_x = x_org[order_idx],
    original_y = y[order_idx],
    increasing = increasing,
    transform_x = transform_x,
    ratio_y = ratio_y,
    precalculated_ratios = precalculated_ratios
  )

  if (test_significance) {
    # Check for k* (when no pooling, k* = k, otherwise k* < k)
    k_star = length(unique(y_final))

    # F-test using original order
    null_model = mean(y)
    null_pred = rep(null_model, length(y))
    sse_full = sum((y - y_pred_new)^2)
    sse_null = sum((y - null_pred)^2)
    df_full = length(y) - length(unique(x))
    df_null = length(y) - 1
    f_stat = ((sse_null - sse_full) / (df_null - df_full)) / (sse_full / df_full)
    p_value = 1 - stats::pf(f_stat, df_null - df_full, df_full)
    f_res = data.frame(
      SSE_Full = sse_full,
      SSE_Null = sse_null,
      F_statistic = f_stat,
      P_value = p_value,
      stringsAsFactors = FALSE
    )
    result$f_test = f_res
  }

  class(result) = "isotonic_model"
  return(result)
}
