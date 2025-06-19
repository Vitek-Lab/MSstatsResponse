#' Predict IC50 (dose where response is 0.5) for each protein and drug
#'
#' @param data A data frame with columns: protein, drug, dose, response.
#' @param increasing Logical. If TRUE, fit a non-decreasing trend.
#' @param ratio_response Logical. If TRUE, use ratio response; else use log2 scale.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1) transformation. Default is TRUE.
#' @param bootstrap Logical. If FALSE, skip confidence interval bootstrap estimation and only return IC50.
#' @param n_samples Number of bootstrap samples. Default = 1000.
#' @param alpha Confidence level. Default = 0.05.
#' @return A data frame with columns: protein, drug, IC50, lower CI, upper CI.
#' @export
PredictIC50 = function(data,
                       n_samples = 1000,
                       alpha = 0.10,
                       increasing = FALSE,
                       transform_dose = TRUE,
                       ratio_response = TRUE,
                       bootstrap = TRUE) {

  protein_list = unique(data$protein)
  drug_list = unique(data$drug[data$drug != "DMSO"])
  results_list = list()

  for (drug_type in drug_list) {
    data_subset = data %>% dplyr::filter(drug %in% c("DMSO", drug_type))

    for (prot in protein_list) {
      tryCatch({
        suppressWarnings({
          df = data_subset %>% dplyr::filter(protein == prot)
          x = df$dose
          y = df$response

          order_idx = order(x)
          x = x[order_idx]
          y = y[order_idx]

          if (ratio_response) {
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
              results_list[[length(results_list) + 1]] = data.frame(
                protein = prot,
                drug = drug_type,
                IC50 = NA,
                IC50_lower_bound = NA,
                IC50_upper_bound = NA
              )
              next
            }

            ic50 = 10^ic50_est

            if (bootstrap) {
              bootstrap_res = bootstrap_ic50(
                dose = x, response = y_ratio,
                n_samples = n_samples, alpha = alpha,
                increasing = increasing
              )
              lower = as.numeric(bootstrap_res$ci_lower_transform)
              upper = as.numeric(bootstrap_res$ci_upper_transform)
            } else {
              lower = NA
              upper = NA
            }

          } else {
            dmso_mean = mean(y[x == 0], na.rm = TRUE)
            target_response = dmso_mean - 1

            fit_try = fit_isotonic_regression(x, y,
                                              increasing = increasing,
                                              transform_x = transform_dose,
                                              ratio_y = FALSE,
                                              test_significance = FALSE)
            ic50_est = predict_ic50(fit_try, target_response = target_response)

            if (is.na(ic50_est)) {
              results_list[[length(results_list) + 1]] = data.frame(
                protein = prot,
                drug = drug_type,
                IC50 = NA,
                IC50_lower_bound = NA,
                IC50_upper_bound = NA
              )
              next
            }

            ic50 = 10^ic50_est

            if (bootstrap) {
              bootstrap_res = bootstrap_ic50_logscale(
                x = x, y = y,
                n_samples = n_samples, alpha = alpha,
                increasing = increasing
              )
              lower = as.numeric(bootstrap_res$ci_lower_transform)
              upper = as.numeric(bootstrap_res$ci_upper_transform)
            } else {
              lower = NA
              upper = NA
            }
          }

          results_list[[length(results_list) + 1]] = data.frame(
            protein = prot,
            drug = drug_type,
            IC50 = ic50,
            IC50_lower_bound = lower,
            IC50_upper_bound = upper
          )
        })
      }, error = function(e) {
        cat("ERROR for", prot, "with", drug_type, ":", conditionMessage(e), "\n")
      })
    }
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
  ic50_vals = numeric(n_samples)
  df = data.frame(dose = dose, response = response)

  for (i in seq_len(n_samples)) {
    # Sample with replacement within each dose group
    boot_df = df %>%
      dplyr::group_by(dose) %>%
      dplyr::sample_frac(size = 1, replace = TRUE) %>%
      dplyr::ungroup()

    x_sample = boot_df$dose
    y_sample = boot_df$response

    tryCatch({
      fit_sample = fit_isotonic_regression(x_sample, y_sample,
                                           increasing = increasing,
                                           transform_x = TRUE,
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
  ic50_vals = numeric(n_samples)
  df = data.frame(dose = x, response = y)

  for (i in seq_len(n_samples)) {
    # Sample with replacement within each dose group
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
