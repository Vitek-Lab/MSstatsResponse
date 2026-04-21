#' Plot isotonic regression fit with optional IC50 for a single protein and drug
#'
#' @param data Protein-level dataset (e.g., output of MSstatsPrepareDoseResponseFit).
#' @param protein_name Character. Protein name to plot.
#' @param drug_name Character. Drug name to plot.
#' @param weights Optional numeric vector of weights matching nrow(data). Default is NULL (equal weights).
#' @param show_weights Logical. If TRUE and weights are provided, scale point size by weight. Default is FALSE.
#' @param increasing Logical or character. If TRUE, fits a non-decreasing model. If FALSE, fits non-increasing.
#'   If "both", fits both and selects better fit.
#' @param ratio_response Logical. If TRUE, compute IC50 on ratio scale; if FALSE, use log2 intensities.
#' @param precalculated_ratios Logical. If TRUE, assumes response contains pre-calculated ratios. Default is FALSE.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1). Default is TRUE.
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param target_response Numeric. Target response level for IC50/EC50 calculation. Default = 0.5.
#' @param add_ci Logical. Include IC50 95% confidence interval bands if TRUE. Default is FALSE.
#' @param n_samples Number of bootstrap samples if including confidence intervals. Default is 1000.
#' @param alpha Alpha level for confidence intervals. Default is 0.10.
#' @param color_by Character. Column name to color points by (e.g., "replicate", "peptide"). Default is NULL.
#' @param x_lab Character. Custom label for the x-axis. If NULL, uses default based on transform_dose.
#' @param y_lab Character. Custom label for the y-axis. If NULL, uses default based on ratio_response.
#' @param title Character. Custom plot title. If NULL, auto-generates from drug and protein names.
#'
#' @return A ggplot object.
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
#' # Example 1: Basic dose-response visualization
#' plot1 <- visualizeResponseProtein(
#'   data = prepared_data,
#'   protein_name = "PROTEIN_A",
#'   drug_name = "Drug1",
#'   ratio_response = TRUE,
#'   show_ic50 = FALSE,
#'   add_ci = FALSE
#' )
#' print(plot1)
#'
#' # Example 2: With weights (visualized as point size)
#' plot2 <- visualizeResponseProtein(
#'   data = prepared_data,
#'   protein_name = "PROTEIN_A",
#'   drug_name = "Drug1",
#'   weights = prepared_data$weight,
#'   show_weights = TRUE,
#'   ratio_response = TRUE,
#'   show_ic50 = TRUE
#' )
#' print(plot2)
#'
#' # Example 3: Color points by replicate
#' \dontrun{
#' plot3 <- visualizeResponseProtein(
#'   data = prepared_data,
#'   protein_name = "PROTEIN_A",
#'   drug_name = "Drug1",
#'   ratio_response = TRUE,
#'   show_ic50 = TRUE,
#'   color_by = "replicate"
#' )
#' print(plot3)
#' }
#'
#' @export
visualizeResponseProtein = function(data,
                                    protein_name,
                                    drug_name,
                                    weights = NULL,
                                    show_weights = FALSE,
                                    ratio_response = TRUE,
                                    precalculated_ratios = FALSE,
                                    transform_dose = TRUE,
                                    show_ic50 = TRUE,
                                    target_response = 0.5,
                                    add_ci = FALSE,
                                    n_samples = 1000,
                                    alpha = 0.10,
                                    increasing = FALSE,
                                    color_by = NULL,
                                    x_lab = NULL,
                                    y_lab = NULL,
                                    title = NULL) {
  # Subset to single protein + drug (+ DMSO)
  data_single = data %>%
    dplyr::filter(protein == protein_name & drug %in% c("DMSO", drug_name))
  
  x_org = data_single$dose
  y_org = data_single$response
  
  # Extract weights for this subset
  if (!is.null(weights)) {
    # Assume weights vector matches the full data
    subset_idx = which(data$protein == protein_name & 
                         data$drug %in% c("DMSO", drug_name))
    w = weights[subset_idx]
  } else {
    w = rep(1, length(x_org))
  }
  
  order_idx = order(x_org)
  x = x_org[order_idx]
  y = y_org[order_idx]
  w = w[order_idx]
  
  # Handle "both" option by fitting and selecting better model
  if (is.character(increasing) && increasing == "both") {
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
    
    # Select better fit
    if (fit_dec$f_test$SSE_Full <= fit_inc$f_test$SSE_Full) {
      fit = fit_dec
      increasing_final = FALSE
    } else {
      fit = fit_inc
      increasing_final = TRUE
    }
  } else {
    # Fit isotonic model
    fit = fitIsotonicRegression(
      x = x, y = y, w = w,
      increasing = increasing,
      transform_x = transform_dose,
      ratio_y = ratio_response,
      precalculated_ratios = precalculated_ratios,
      test_significance = FALSE
    )
    increasing_final = increasing
  }
  
  # Optional: bootstrap CIs
  ci_bounds = NULL
  if (add_ci) {
    suppressWarnings({
      ic50_est = predictIC50(data_single,
                             weights = w,  # Pass weights to IC50 calculation
                             ratio_response = ratio_response,
                             precalculated_ratios = precalculated_ratios,
                             transform_dose = transform_dose,
                             increasing = increasing_final,
                             target_response = target_response,
                             n_samples = n_samples,
                             alpha = alpha)
    })
    ci_bounds = list(
      ci_lower = ic50_est$IC50_lower_bound,
      ci_upper = ic50_est$IC50_upper_bound
    )
  }
  
  p = plotIsotonic(
    fit = fit,
    weights = w,
    show_weights = show_weights,
    ratio = ratio_response,
    precalculated_ratios = precalculated_ratios,
    show_ic50 = show_ic50,
    target_response = target_response,
    ci = ci_bounds,
    drug_name = drug_name,
    protein_name = protein_name,
    x_lab = x_lab,
    y_lab = y_lab,
    title = title,
    transform_x = transform_dose,
    color_by = color_by,
    original_data = data_single
  )
  
  return(p)
}

#' Plot Isotonic Regression Model
#'
#' @param fit A model object returned by fitIsotonicRegression().
#' @param weights Numeric vector of weights (same length as data points). Default is NULL.
#' @param show_weights Logical. If TRUE and weights provided, scale point size by weight. Default is FALSE.
#' @param ratio Logical. If TRUE, shows plot on the ratio scale relative to DMSO. Default is FALSE.
#' @param precalculated_ratios Logical. If TRUE, response values are pre-calculated ratios. Default is FALSE.
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param target_response Numeric. Target response level for IC50/EC50. Default = 0.5.
#' @param drug_name Drug name for plotting data.
#' @param protein_name Protein name for plot.
#' @param x_lab Character. Label for x-axis. If NULL, uses default based on transform_x.
#' @param y_lab Character. Label for y-axis. If NULL, uses default based on ratio/precalculated_ratios.
#' @param title Character. Title for the plot. If NULL, auto-generates.
#' @param ci List with ci_lower and ci_upper. Include IC50 confidence interval bands if provided. Default is NULL.
#' @param legend Logical. Show legend if TRUE.
#' @param theme_style ggplot2 theme name to apply (default = "classic").
#' @param original_label Logical. If TRUE, replace x-axis tick labels with original dose labels.
#' @param transform_x Logical. Whether doses were log-transformed. Default = TRUE.
#' @param color_by Character. Column name to color points by. Default is NULL.
#' @param original_data Data frame containing original data with grouping variables.
#'
#' @return A ggplot object.
#'
#' @importFrom stats approx
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment annotate
#' @importFrom ggplot2 scale_color_manual scale_x_continuous labs ylim theme_minimal
#' @importFrom ggplot2 theme element_rect element_blank element_text scale_size_continuous
plotIsotonic = function(fit,
                        weights = NULL,
                        show_weights = FALSE,
                        ratio = TRUE,
                        precalculated_ratios = FALSE,
                        show_ic50 = FALSE,
                        target_response = 0.5,
                        drug_name = NULL,
                        protein_name = NULL,
                        x_lab = NULL,
                        y_lab = NULL,
                        title = NULL,
                        ci = NULL,
                        legend = FALSE,
                        theme_style = "classic",
                        original_label = FALSE,
                        transform_x = TRUE,
                        color_by = NULL,
                        original_data = NULL) {
  
  # Construct title if not provided
  if (is.null(title)) {
    if (!is.null(drug_name) && !is.null(protein_name)) {
      title = paste0(drug_name, ": ", protein_name, " Response")
    } else {
      title = "Isotonic Regression Fit"
    }
  }
  
  # Set default labels if not provided
  if (is.null(y_lab)) {
    if (precalculated_ratios || ratio) {
      y_lab = "Ratio relative to DMSO"
    } else {
      y_lab = expression(Log[2] ~ "Intensity")
    }
  }
  
  if (is.null(x_lab)) {
    if (transform_x) {
      x_lab = expression(Log[10] ~ "[drug (M)]")
    } else {
      x_lab = "Dose"
    }
  }
  
  # Prepare data
  if (transform_x) {
    fit_df = data.frame(dose = log10(fit$x), y_pred = fit$y_pred)
    orig_df = data.frame(x = log10(fit$x), y = fit$original_y)
  } else {
    fit_df = data.frame(dose = fit$x, y_pred = fit$y_pred)
    orig_df = data.frame(x = fit$x, y = fit$original_y)
  }
  
  # Add weights if provided
  if (!is.null(weights) && show_weights) {
    orig_df$weight = weights
  }
  
  # Add color grouping variable if requested
  if (!is.null(color_by) && !is.null(original_data)) {
    if (!color_by %in% names(original_data)) {
      stop(paste("Column", color_by, "not found in data"))
    }
    order_idx = order(original_data$dose)
    orig_df[[color_by]] = original_data[[color_by]][order_idx]
  }
  
  # Base plot with conditional coloring and sizing
  if (!is.null(color_by) && !is.null(original_data)) {
    # Color points by grouping variable
    if (!is.null(weights) && show_weights) {
      model_fit_plot = ggplot2::ggplot() +
        ggplot2::geom_point(data = orig_df,
                            ggplot2::aes(x = x, y = y, color = .data[[color_by]], alpha = weight, size = weight)) +
  ggplot2::scale_alpha_continuous(range = c(0.35, 1), name = "Weight") +
  ggplot2::scale_size_continuous(range = c(1.5, 3), name = "Weight") +
  ggplot2::guides(alpha = "legend", size = "legend")
    } else {
      model_fit_plot = ggplot2::ggplot() +
        ggplot2::geom_point(data = orig_df,
                            ggplot2::aes(x = x, y = y, color = .data[[color_by]]),
                            size = 2.5)
    }
    model_fit_plot = model_fit_plot +
      ggplot2::geom_line(data = fit_df,
                         ggplot2::aes(x = dose, y = y_pred),
                         color = "orange", size = 1.2) +
      ggplot2::labs(x = x_lab, y = y_lab, title = title, color = color_by)
  } else {
    # Default: black points, orange line
    if (!is.null(weights) && show_weights) {
      model_fit_plot = ggplot2::ggplot() +
        ggplot2::geom_point(data = orig_df,
                            ggplot2::aes(x = x, y = y, size = weight),
                            color = "black") +
        ggplot2::scale_size_continuous(range = c(1, 5), name = "Weight")
    } else {
      model_fit_plot = ggplot2::ggplot() +
        ggplot2::geom_point(data = orig_df,
                            ggplot2::aes(x = x, y = y, color = "Observed Data"),
                            size = 2)
    }
    model_fit_plot = model_fit_plot +
      ggplot2::geom_line(data = fit_df,
                         ggplot2::aes(x = dose, y = y_pred, color = "Isotonic Regression Fit"),
                         size = 1.2)
    
    # Only add color scale if not showing weights
    if (!show_weights) {
      model_fit_plot = model_fit_plot +
        ggplot2::scale_color_manual(name = NULL,
                                    values = c("Observed Data" = "black",
                                               "Isotonic Regression Fit" = "orange"))
    }
    
    model_fit_plot = model_fit_plot +
      ggplot2::labs(x = x_lab, y = y_lab, title = title)
  }
  
  # Apply common theme
  model_fit_plot = model_fit_plot +
    ggplot2::ylim(0, NA) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  # Handle legend display
  if (!legend && is.null(color_by) && !show_weights) {
    model_fit_plot = model_fit_plot +
      ggplot2::theme(legend.position = "none")
  }
  
  # Optional IC50
  if (show_ic50) {
    if (precalculated_ratios || ratio) {
      adjusted_target = target_response
    } else {
      if (transform_x) {
        dmso_idx = which(is.infinite(orig_df$x) & orig_df$x < 0)
      } else {
        dmso_idx = which(fit$original_x == 0)
      }
      adjusted_target = mean(orig_df$y[dmso_idx], na.rm = TRUE) + log2(target_response)
    }
    
    ic50_pred = PredictIC50(fit, target_response = adjusted_target)
    
    # Handle duplicate x values for approx
    if (length(unique(fit$x)) < length(fit$x)) {
      x_unique = unique(fit$x)
      y_aggregated = sapply(x_unique, function(xi) {
        mean(fit$y_pred[fit$x == xi])
      })
      y_ic50 = stats::approx(x_unique, y_aggregated, xout = ic50_pred)$y
    } else {
      suppressWarnings({
        y_ic50 = stats::approx(fit$x, fit$y_pred, xout = ic50_pred)$y
      })
    }
    
    if (transform_x) {
      ic50_pred_transform = -log10(ic50_pred)
      ic50_label = paste("pIC50 =", round(ic50_pred_transform, 2))
      x_pos = log10(ic50_pred)
    } else {
      ic50_label = paste("IC50 =", round(ic50_pred, 2))
      x_pos = ic50_pred
    }
    
    model_fit_plot = model_fit_plot +
      ggplot2::geom_point(ggplot2::aes(x = x_pos, y = y_ic50),
                          shape = 23, size = 3.5, fill = "red", color = "black") +
      ggplot2::annotate("text", x = x_pos + 0.35, y = y_ic50,
                        label = ic50_label,
                        vjust = -0.5, hjust = 0.5, size = 4, fontface = "bold")
    
    # Optional display IC50 confidence interval
    if (show_ic50 && !is.null(ci) && !is.na(ci$ci_lower) && !is.na(ci$ci_upper)) {
      if (transform_x) {
        ci_lower_dose = 10^(-ci$ci_lower)
        ci_upper_dose = 10^(-ci$ci_upper)
        ci_lower_x = log10(ci_lower_dose)
        ci_upper_x = log10(ci_upper_dose)
      } else {
        ci_lower_x = ci$ci_lower
        ci_upper_x = ci$ci_upper
      }
      
      model_fit_plot = model_fit_plot +
        ggplot2::geom_segment(ggplot2::aes(x = ci_lower_x, xend = ci_lower_x,
                                           y = 0, yend = adjusted_target),
                              linetype = "dotted", color = "gray40", linewidth = 0.8) +
        ggplot2::geom_segment(ggplot2::aes(x = ci_upper_x, xend = ci_upper_x,
                                           y = 0, yend = adjusted_target),
                              linetype = "dotted", color = "gray40", linewidth = 0.8) +
        ggplot2::annotate("rect",
                          xmin = ci_lower_x, xmax = ci_upper_x,
                          ymin = 0, ymax = adjusted_target,
                          alpha = 0.1, fill = "red")
    }
  }
  
  return(model_fit_plot)
}
