#' Drug-protein interaction detection tested by F-test (fitted curve vs average response)
#'
#' Fits an isotonic regression model to protein abundance data.
#' Performs an F-test to assess the significance of the dose-response curve and applies FDR correction.
#'
#' @param data Protein-level data, formatted with MSstatsPrepareDoseResponseFit().
#' @param weights Optional numeric vector of weights. Defaults to equal weights.
#' @param increasing Logical. If TRUE, fits a non-decreasing model. If FALSE, fits non-increasing.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default is TRUE.
#' @param ratio_response Logical. If TRUE, converts log2 abundance to ratios relative to DMSO. Default is FALSE.
#'
#' @return A data frame with protein-wise F-test results and BH-adjusted p-values.
#' @export
#' @importFrom stats pf p.adjust
#' @import dplyr

DoseResponseFit = function(data, weights = NULL,
                           increasing = FALSE,
                           transform_dose = TRUE,
                           ratio_response = FALSE) {

  drug_list = unique(data$drug[data$drug != "DMSO"])
  all_results = list()

  for (drug_type in drug_list) {
    drug_subset = data %>% filter(drug %in% c("DMSO", drug_type))
    protein_list = unique(drug_subset$protein)
    results_list = list()

    for (i in seq_along(protein_list)) {
      tryCatch({
        suppressWarnings({
          data_single_protein = drug_subset %>% filter(protein == protein_list[i])

          x = data_single_protein$dose
          y = data_single_protein$response
          w = if (is.null(weights)) rep(1, length(y)) else weights

          fit = fit_isotonic_regression(
            x = x, y = y, w = w,
            increasing = increasing,
            transform_x = transform_dose,
            ratio_y = ratio_response,
            test_significance = TRUE
          )

          results_temp = fit$f_test
          results_temp$protein = protein_list[i]
          results_temp$drug = drug_type
          results_temp = results_temp[, c("protein", "drug", setdiff(names(results_temp), c("protein", "drug")))]
          results_list[[i]] = results_temp
        })
      }, error = function(e) {
        cat("ERROR for drug:", drug_type, "protein:", protein_list[i], ":", conditionMessage(e), "\n")
      })
    }

    # Combine and adjust p-values for this drug
    results_df = do.call(rbind, results_list)
    results_df$adjust_pval = p.adjust(results_df$P_value, method = "BH")
    all_results[[drug_type]] = results_df
  }

  final_results = do.call(rbind, all_results)
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
#' @param test_significance Logical. If TRUE, performs F-test to assess significance.
#'
#' @return A list representing the isotonic regression model (class = "isotonic_model").
#'
#' @export
#' @importFrom stats pf approx quantile

fit_isotonic_regression = function(x, y, w = rep(1, length(y)),
                                   increasing = FALSE,
                                   transform_x = TRUE,
                                   ratio_y = FALSE,
                                   test_significance = FALSE) {
  stopifnot(length(x) == length(y), length(y) == length(w))

  x_org = x

  if (transform_x) {
    x = log10(x + 1)
  }

  if (ratio_y) {
    y = 2^y
    baseline = mean(y[x == 0], na.rm = TRUE)
    y = y / baseline
  }


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
    k = target[i] + 1 # compares to next value in sequence
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

  y_pred_new = approx(x, y_final, xout = x, rule = 2)$y

  result = list(
    x = x,
    #y = y_final,
    y_pred = y_pred_new,
    original_x = x_org[order_idx],
    original_y = y[order_idx],
    increasing = increasing,
    transform_x = transform_x,
    ratio_y = ratio_y
  )

  if (test_significance) {
    # F-test using original order
    null_model = mean(y)
    null_pred = rep(null_model, length(y))
    sse_full = sum((y - y_pred_new)^2)
    sse_null = sum((y - null_pred)^2)
    df_full = length(y) - length(unique(x))
    df_null = length(y) - 1
    f_stat = ((sse_null - sse_full) / (df_null - df_full)) / (sse_full / df_full)
    p_value = 1 - pf(f_stat, df_null - df_full, df_full)
    f_res = data.frame(
      SSE_Full = sse_full,
      SSE_Null = sse_null,
      F_statistic = f_stat,
      P_value = p_value
    )
    result$f_test = f_res
  }

  class(result) = "isotonic_model"
  return(result)
}
