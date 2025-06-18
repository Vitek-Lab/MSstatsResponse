#' Bootstrap IC50 Estimates and Confidence Interval
#'
#' @param dose Numeric vector of dose values.
#' @param response Numeric vector of response values.
#' @param n_samples Number of bootstrap samples (default = 1000).
#' @param alpha Significance level for confidence interval (default = 0.10).
#' @param ... Additional arguments passed to fit_isotonic_regression().
#'
#' @return List with mean IC50, lower and upper CI bounds, and all bootstrapped estimates.
#'
#' @importFrom stats pf approx quantile

bootstrap_ic50 = function(dose, response, n_samples = 1000, alpha = 0.10,
                          increasing = FALSE) {
  n = length(dose)
  ic50_vals = numeric(n_samples)

  for (i in seq_len(n_samples)) {
    # Resample data with replacement
    idx = sample(seq_len(n), n, replace = TRUE)
    dose_sample = dose[idx]
    response_sample = response[idx]

    # Fit model and predict IC50 for bootstrap sample
    tryCatch({
      fit_sample = fit_isotonic_regression(dose_sample, response_sample,
                                           increasing = increasing,
                                           transform_dose = TRUE,
                                           ratio_response = TRUE,
                                           test_significance = FALSE)
      ic50_est = predict_ic50(fit_sample)

      # Only store IC50 if it's not NA
      if (!is.na(ic50_est)) {
        ic50_vals[i] = ic50_est
      } else {
        ic50_vals[i] = NA
      }
    }, error = function(e) {
      ic50_vals[i] = NA
    })
  }

  # Remove missing or invalid estimates
  ic50_values_clean = ic50_vals[!is.na(ic50_vals) & ic50_vals != 0]

  # Compute bootstrap summary statistics
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
