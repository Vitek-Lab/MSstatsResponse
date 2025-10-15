.calcSingleIC50 = function(df,
                           n_samples,
                           alpha,
                           increasing,
                           transform_dose,
                           ratio_response,
                           bootstrap,
                           prot,
                           drug_type,
                           target_response){


  x = df$dose
  y = df$response

  # Check for DMSO presence (dose = 0)
  if (!any(x == 0)) {
    warning(paste("No DMSO control for protein", prot, "and drug", drug_type))
    return(data.frame(
      protein = prot,
      drug = drug_type,
      IC50 = NA,
      IC50_lower_bound = NA,
      IC50_upper_bound = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Validate input data
  if (nrow(df) < 2) {
    warning(paste("Insufficient data for protein", prot, "and drug", drug_type))
    return(data.frame(
      protein = prot,
      drug = drug_type,
      IC50 = NA,
      IC50_lower_bound = NA,
      IC50_upper_bound = NA,
      stringsAsFactors = FALSE
    ))
  }

  order_idx = order(x)
  x = x[order_idx]
  y = y[order_idx]

  y_unlog = 2^y
  baseline = mean(y_unlog[x == 0], na.rm = TRUE)

  # Check if baseline is valid
  if (is.na(baseline) || !is.finite(baseline)) {
    warning(paste("Invalid baseline for protein", prot, "and drug", drug_type))
    return(data.frame(
      protein = prot,
      drug = drug_type,
      IC50 = NA,
      IC50_lower_bound = NA,
      IC50_upper_bound = NA,
      stringsAsFactors = FALSE
    ))
  }


  y_ratio = y_unlog / baseline

  if (!ratio_response){
    dmso_mean = mean(y[x == 0], na.rm = TRUE)
    adjusted_target_response = dmso_mean + log2(target_response)
  } else {
    adjusted_target_response = target_response
  }

  result = NULL
  ic50_est = NA

  tryCatch({
    fit_try = fitIsotonicRegression(x, y,
                                    increasing = increasing,
                                    transform_x = transform_dose,
                                    ratio_y = ratio_response,
                                    test_significance = FALSE)
    ic50_est = PredictIC50(fit_try, target_response = adjusted_target_response)
    ic50_est = ifelse(is.na(ic50_est), NA, ic50_est)
  }, error = function(e) {
    ic50_est <- NA
  })

  if (is.na(ic50_est)) {
    result = data.frame(
      protein = prot,
      drug = drug_type,
      IC50 = NA,
      IC50_lower_bound = NA,
      IC50_upper_bound = NA
    )
  } else {
    ic50 = -log10(ic50_est)

    if (bootstrap) {
      if(ratio_response){
        bootstrap_res = bootstrapIC50(
          dose = x, response = y,
          n_samples = n_samples, alpha = alpha,
          increasing = increasing,
          target_response = target_response
        )
      } else {
        bootstrap_res = bootstrapIC50LogScale(
          x = x, y = y,
          n_samples = n_samples, alpha = alpha,
          increasing = increasing,
          target_response = target_response
        )
      }
      lower = as.numeric(bootstrap_res$ci_lower_transform)
      upper = as.numeric(bootstrap_res$ci_upper_transform)
    } else {
      lower = NA
      upper = NA
    }
    result = data.frame(
      protein = prot,
      drug = drug_type,
      IC50 = ic50,
      IC50_lower_bound = lower,
      IC50_upper_bound = upper
    )
  }

  return(result)
}

#' Predict IC50 (dose where response is 0.5 of DMSO)
#'
#' @param fit An object of class "isotonic_model" from fitIsotonicRegression().
#' @param target_response The target response value for prediction (default = 0.5).
#'
#' @return Estimated dose value corresponding to target response (IC50)
#' @importFrom stats approx
PredictIC50 = function(fit, target_response = 0.5) {
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
#' @param target_response Numeric value for response level (default = 0.5).
#' @param transform_x Logical. If TRUE, applies log10(dose + 1) transformation. Default = TRUE.
#'
#' @return List with mean IC50, CI bounds, and transformed estimates.
#' @importFrom stats quantile
#' @importFrom dplyr filter select mutate group_by summarise arrange distinct

bootstrapIC50 = function(dose, response, n_samples = 1000, alpha = 0.10,
                         increasing = FALSE, target_response = 0.5,
                         transform_x = TRUE) {
  ic50_vals = numeric(n_samples)
  df = data.frame(dose = dose, response = response)

  for (i in seq_len(n_samples)) {
    boot_df = df %>%
      dplyr::group_by(dose) %>%
      dplyr::sample_frac(size = 1, replace = TRUE) %>%
      dplyr::ungroup()

    x_sample = boot_df$dose
    y_sample = boot_df$response

    tryCatch({
      fit_sample = fitIsotonicRegression(x_sample, y_sample,
                                         increasing = increasing,
                                         transform_x = transform_x,
                                         ratio_y = TRUE,
                                         test_significance = FALSE)
      ic50_est = PredictIC50(fit_sample, target_response = target_response)
      ic50_vals[i] = ifelse(is.na(ic50_est), NA, ic50_est)
    }, error = function(e) {
      ic50_vals[i] = NA
    })
  }

  ic50_values_clean = ic50_vals[!is.na(ic50_vals) & ic50_vals != 0]
  ci_lower = stats::quantile(ic50_values_clean, alpha / 2, na.rm = TRUE)
  ci_upper = stats::quantile(ic50_values_clean, 1 - alpha / 2, na.rm = TRUE)
  mean_ic50 = mean(ic50_values_clean, na.rm = TRUE)

  return(list(
    ic50_values = ic50_vals,
    mean_ic50 = mean_ic50,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_ic50_transform = -log10(mean_ic50),
    ci_lower_transform = -log10(ci_lower),
    ci_upper_transform = -log10(ci_upper)
  ))
}

#' Bootstrap IC50 Estimates and Confidence Interval (log scale)
#'
#' @param x Numeric vector of dose values.
#' @param y Numeric vector of log2 response values.
#' @param n_samples Number of bootstrap samples (default = 1000).
#' @param alpha Significance level for confidence interval (default = 0.05).
#' @param increasing Logical. Fit non-decreasing if TRUE.
#' @param target_response Numeric value for response level (default = 0.5).
#'
#' @return List with mean IC50, CI bounds, and transformed estimates.
#' @importFrom stats quantile
#' @importFrom dplyr filter select mutate group_by summarise arrange distinct
bootstrapIC50LogScale = function(x, y, n_samples = 1000, alpha = 0.05,
                                 increasing = FALSE, target_response = 0.5) {
  ic50_vals = numeric(n_samples)
  df = data.frame(dose = x, response = y)

  for (i in seq_len(n_samples)) {
    boot_df = df %>%
      dplyr::group_by(dose) %>%
      dplyr::sample_frac(size = 1, replace = TRUE) %>%
      dplyr::ungroup()

    x_sample = boot_df$dose
    y_sample = boot_df$response

    dmso_idx = which(x_sample == 0)
    if (length(dmso_idx) == 0) {
      ic50_vals[i] = NA
      next
    }

    mean_dmso = mean(y_sample[dmso_idx], na.rm = TRUE)
    adjusted_target_response = mean_dmso + log2(target_response)

    tryCatch({
      fit_sample = fitIsotonicRegression(x_sample, y_sample,
                                         increasing = increasing,
                                         transform_x = TRUE,
                                         ratio_y = FALSE,
                                         test_significance = FALSE)
      ic50_est = PredictIC50(fit_sample, target_response = adjusted_target_response)
      ic50_vals[i] = ifelse(is.na(ic50_est), NA, ic50_est)
    }, error = function(e) {
      ic50_vals[i] = NA
    })
  }

  ic50_values_clean = ic50_vals[!is.na(ic50_vals) & ic50_vals != 0]
  ci_lower = stats::quantile(ic50_values_clean, alpha / 2, na.rm = TRUE)
  ci_upper = stats::quantile(ic50_values_clean, 1 - alpha / 2, na.rm = TRUE)
  mean_ic50 = mean(ic50_values_clean, na.rm = TRUE)

  return(list(
    ic50_values = ic50_vals,
    mean_ic50 = mean_ic50,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_ic50_transform = -log10(mean_ic50),
    ci_lower_transform = -log10(ci_lower),
    ci_upper_transform = -log10(ci_upper)
  ))
}
