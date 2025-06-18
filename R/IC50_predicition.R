#' Estimate IC50 and Confidence Intervals from Isotonic Regression
#'
#' Fits isotonic regression models to protein-level dose-response data
#' and estimates IC50 along with confidence intervals using bootstrapping.
#'
#' @param data A data.frame formatted with MSstatsPrepareDoseResponseFit().
#' @param n_samples Number of bootstrap samples for CI estimation (default = 1000).
#' @param alpha Significance level for confidence intervals (default = 0.10).
#' @param increasing Logical. If TRUE, fits a non-decreasing model. Default is FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default is TRUE.
#' @param ratio_response Logical. If TRUE, transforms y to ratios relative to DMSO. Default is TRUE.
#'
#' @return A data.frame with protein, IC50, IC50_lower_bound, IC50_upper_bound
#' @export
#' @importFrom dplyr filter

PredictIC50 = function(data,
                       n_samples = 1000,
                       alpha = 0.10,
                       increasing = FALSE,
                       transform_dose = TRUE,
                       ratio_response = TRUE) {

  protein_list = unique(data$protein)
  results_list = list()

  for (i in seq_along(protein_list)) {
    tryCatch({
      suppressWarnings({
        df = data %>% filter(protein == protein_list[i])
        x = df$dose
        y = df$response

        # Sort x and y
        order_idx = order(x)
        x = x[order_idx]
        y = y[order_idx]

        if (ratio_response) {
          # ratio-based response
          y_unlog = 2^y
          baseline = mean(y_unlog[x == 0], na.rm = TRUE)
          y_ratio = y_unlog / baseline

          fit_try = fit_isotonic_regression(x, y,
                                            increasing = increasing,
                                            transform_x = transform_dose,
                                            ratio_y = TRUE,
                                            test_significance = FALSE)

          ic50_est = predict_ic50(fit_try, target_response = 0.5)

          if (is.na(ic50_est)) {
            results_list[[i]] = data.frame(
              protein = protein_list[i],
              IC50 = NA,
              IC50_lower_bound = NA,
              IC50_upper_bound = NA
            )
            next
          }

          bootstrap_res = bootstrap_ic50(
            x = x, y = y_ratio,
            n_samples = n_samples, alpha = alpha,
            increasing = increasing
          )

          ic50 = 10^ic50_est
          lower = as.numeric(bootstrap_res$ci_lower_transform)
          upper = as.numeric(bootstrap_res$ci_upper_transform)

        } else {
          # log-based response
          dmso_mean = mean(y[x == 0], na.rm = TRUE)
          target_response = dmso_mean - 1

          fit_try = fit_isotonic_regression(x, y,
                                            increasing = increasing,
                                            transform_x = transform_dose,
                                            ratio_y = FALSE,
                                            test_significance = FALSE)

          ic50_est = predict_ic50(fit_try, target_response = target_response)

          if (is.na(ic50_est)) {
            results_list[[i]] = data.frame(
              protein = protein_list[i],
              IC50 = NA,
              IC50_lower_bound = NA,
              IC50_upper_bound = NA
            )
            next
          }

          bootstrap_res = bootstrap_ic50_logscale(
            x = x, y = y,
            n_samples = n_samples, alpha = alpha,
            increasing = increasing
          )

          ic50 = 10^ic50_est
          lower = as.numeric(bootstrap_res$ci_lower_transform)
          upper = as.numeric(bootstrap_res$ci_upper_transform)
        }

        # Store result
        results_list[[i]] = data.frame(
          protein = protein_list[i],
          IC50 = ic50,
          IC50_lower_bound = lower,
          IC50_upper_bound = upper
        )

      })  # end suppressWarnings
    }, error = function(e) {
      cat("ERROR for protein:", protein_list[i], ":", conditionMessage(e), "\n")
    })
  }

  final_df = do.call(rbind, results_list)
  return(final_df)
}




#' Predict IC50 (dose where response is 0.5 of DMSO)
#'
#' @param fit An object of class "isotonic_model" from fit_isotonic_regression().
#' @param target_response The target response value for prediction (default = 0.5).
#'
#' @return Estimated dose value corresponding to target response (IC50)
#' @export
predict_ic50 = function(fit, target_response = 0.5) {
  dose = fit$x
  response = fit$y_pred

  # Find all points where the response exactly matches the target
  above = which(response >= target_response)
  below = which(response <= target_response)
  common = intersect(above, below)

  if (length(common) > 0) {
    return(dose[common[1]])
  }

  # Linearly interpolate between nearest points where response crosses target
  for (i in seq_along(response)[-length(response)]) {
    if ((response[i] < target_response && response[i + 1] > target_response) ||
        (response[i] > target_response && response[i + 1] < target_response)) {
      dose1 = dose[i]; dose2 = dose[i + 1]
      response1 = response[i]; response2 = response[i + 1]
      slope = (response2 - response1) / (dose2 - dose1)
      return(dose1 + (target_response - response1) / slope)
    }
  }

  warning("Target response value not reached by fit")
  return(NA)
}

#' Bootstrap IC50 Estimates and Confidence Interval (ratio scale)
#'
#' @param dose Numeric vector of dose values.
#' @param response Numeric vector of response values (on log2 scale).
#' @param n_samples Number of bootstrap samples (default = 1000).
#' @param alpha Significance level for confidence interval (default = 0.10).
#' @param increasing Logical. Fit non-decreasing if TRUE.
#'
#' @return List with mean IC50, CI bounds, and transformed estimates.
#' @export
bootstrap_ic50 = function(dose, response, n_samples = 1000, alpha = 0.10,
                          increasing = FALSE) {
  n = length(dose)
  ic50_vals = numeric(n_samples)

  for (i in seq_len(n_samples)) {
    idx = sample(seq_len(n), n, replace = TRUE)
    dose_sample = dose[idx]
    response_sample = response[idx]

    tryCatch({
      fit_sample = fit_isotonic_regression(dose_sample, response_sample,
                                           increasing = increasing,
                                           transform_dose = TRUE,
                                           ratio_y = TRUE,
                                           test_significance = FALSE)
      ic50_est = predict_ic50(fit_sample)
      ic50_vals[i] = ifelse(is.na(ic50_est), NA, ic50_est)
    }, error = function(e) {
      ic50_vals[i] = NA
    })
  }

  ic50_values_clean = ic50_vals[!is.na(ic50_vals) & ic50_vals != 0]
  ci_lower = quantile(ic50_values_clean, alpha / 2, na.rm = TRUE)
  ci_upper = quantile(ic50_values_clean, 1 - alpha / 2, na.rm = TRUE)
  mean_ic50 = mean(ic50_values_clean, na.rm = TRUE)

  return(list(
    ic50_values = ic50_vals,
    mean_ic50 = mean_ic50,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_ic50_transform = 10^mean_ic50,
    ci_lower_transform = 10^ci_lower,
    ci_upper_transform = 10^ci_upper
  ))
}

#' Bootstrap IC50 Estimates and Confidence Interval (log scale)
#'
#' @param x Numeric vector of dose values.
#' @param y Numeric vector of log2 response values.
#' @param n_samples Number of bootstrap samples (default = 1000).
#' @param alpha Significance level for confidence interval (default = 0.05).
#' @param increasing Logical. Fit non-decreasing if TRUE.
#'
#' @return List with mean IC50, CI bounds, and transformed estimates.
#' @export
bootstrap_ic50_logscale = function(x, y, n_samples = 1000, alpha = 0.05,
                                   increasing = FALSE) {
  n = length(x)
  ic50_vals = numeric(n_samples)

  for (i in seq_len(n_samples)) {
    idx = sample(seq_len(n), n, replace = TRUE)
    x_sample = x[idx]
    y_sample = y[idx]

    dmso_idx = which(x_sample == 0)
    if (length(dmso_idx) == 0) {
      ic50_vals[i] = NA
      next
    }

    mean_dmso = mean(y_sample[dmso_idx], na.rm = TRUE)
    target_response = mean_dmso - 1

    tryCatch({
      fit_sample = fit_isotonic_regression(x_sample, y_sample,
                                           increasing = increasing,
                                           transform_x = TRUE,
                                           ratio_y = FALSE,
                                           test_significance = FALSE)
      ic50_est = predict_ic50(fit_sample, target_response = target_response)
      ic50_vals[i] = ifelse(is.na(ic50_est), NA, ic50_est)
    }, error = function(e) {
      ic50_vals[i] = NA
    })
  }

  ic50_values_clean = ic50_vals[!is.na(ic50_vals) & ic50_vals != 0]
  ci_lower = quantile(ic50_values_clean, alpha / 2, na.rm = TRUE)
  ci_upper = quantile(ic50_values_clean, 1 - alpha / 2, na.rm = TRUE)
  mean_ic50 = mean(ic50_values_clean, na.rm = TRUE)

  return(list(
    ic50_values = ic50_vals,
    mean_ic50 = mean_ic50,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_ic50_transform = 10^mean_ic50,
    ci_lower_transform = 10^ci_lower,
    ci_upper_transform = 10^ci_upper
  ))
}
