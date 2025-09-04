#' Plot isotonic regression fit with optional IC50 for a single protein and drug
#'
#' @param data Protein-level dataset (e.g., output of MSstatsPrepareDoseResponseFit).
#' @param protein_name Character. Protein name to plot.
#' @param drug_name Character. Drug name to plot.
#' @param increasing Logical. If TRUE, fits a non-decreasing model. If FALSE, fits non-increasing.
#' @param ratio_response Logical. If TRUE, compute IC50 on ratio scale; if FALSE, use log2 intensities.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1). Default is TRUE.
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param add_ci Logical. Include IC50 95% confidence interval bands if TRUE. Default is FALSE.
#' @param n_samples Number of bootstrap samples if including confidence intervals. Default is 1000.
#' @param alpha Alpha level for confidence intervals. Default is 0.05.
#' @param y_lab Character. Label for the y-axis. Default is "Ratio Response".
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
#' dose_info <- ConvertGroupToNumericDose(dia_data$ProteinLevelData$GROUP)
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
#' plot1 <- VisualizeResponseProtein(
#'   data = prepared_data,
#'   protein_name = "P12931",
#'   drug_name = "Dasatinib",
#'   ratio_response = TRUE,
#'   show_ic50 = FALSE,
#'   add_ci = FALSE
#' )
#' print(plot1)
#'
#' # Example 2: Add IC50 annotation
#' plot2 <- VisualizeResponseProtein(
#'   data = prepared_data,
#'   protein_name = "P12931",
#'   drug_name = "Dasatinib",
#'   ratio_response = TRUE,
#'   show_ic50 = TRUE,
#'   add_ci = FALSE
#' )
#' print(plot2)
#'
#' @export
VisualizeResponseProtein = function(data,
                                    protein_name,
                                    drug_name,
                                    ratio_response = TRUE,
                                    transform_dose = TRUE,
                                    show_ic50 = TRUE,
                                    add_ci = FALSE,
                                    n_samples = 1000,
                                    alpha = 0.10,
                                    increasing = FALSE,
                                    y_lab = 'Ratio Response') {
  # Subset to single protein + drug (+ DMSO)
  data_single = data %>%
    dplyr::filter(protein == protein_name & drug %in% c("DMSO", drug_name))

  x_org = data_single$dose
  y_org = data_single$response

  order_idx = order(x_org)
  x = x_org[order_idx]
  y = y_org[order_idx]

  # Fit isotonic model using renamed function
  fit = fitIsotonicRegression(
    x = x, y = y,
    increasing = increasing,
    transform_x = transform_dose,
    ratio_y = ratio_response,
    test_significance = FALSE
  )

  # Optional: bootstrap CIs
  ci_bounds = NULL
  if (add_ci) {
    suppressWarnings({
    ic50_est = PredictIC50(data_single,
                           ratio_response = ratio_response,
                           transform_dose = transform_dose,
                           increasing = increasing,
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
    ratio = ratio_response,
    show_ic50 = show_ic50,
    ci = ci_bounds,
    drug_name = drug_name,
    protein_name = protein_name,
    y_lab = y_lab
  )

  return(p)
}

#' Plot Isotonic Regression Model
#'
#' @param fit A model object returned by fitIsotonicRegression().
#' @param ratio Logical. If TRUE, shows plot on the ratio scale relative to DMSO (i.e. 0-1 scale). Default is FALSE.
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param drug_name Drug name for plotting data.
#' @param protein_name Protein name for plot.
#' @param x_lab Label for x-axis.
#' @param y_lab Label for y-axis.
#' @param title Title for the plot.
#' @param ci Logical. Include IC50 95% confidence interval bands if TRUE. Default is FALSE.
#' @param legend Logical. Show legend if TRUE.
#' @param theme_style ggplot2 theme name to apply (default = "classic").
#' @param original_label Logical. If TRUE, replace x-axis tick labels with original dose labels.
#'
#' @return A ggplot object.
#'
#' @importFrom stats approx
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment annotate
#' @importFrom ggplot2 scale_color_manual scale_x_continuous labs ylim theme_minimal
#' @importFrom ggplot2 theme element_rect element_blank element_text
plotIsotonic = function(fit,
                        ratio = TRUE,
                        show_ic50 = FALSE,
                        drug_name = NULL,
                        protein_name = NULL,
                        x_lab = expression(Log[10] ~ "[drug (M)]"),
                        y_lab = "Log2 Intensity",
                        title = NULL,
                        ci = NULL,
                        legend = FALSE,
                        theme_style = "classic",
                        original_label = FALSE) {

  # Construct title if not provided
  if (is.null(title) && !is.null(drug_name) && !is.null(protein_name)) {
    title = paste0(drug_name, ": ", protein_name, " Response")
  } else if (is.null(title)) {
    title = "Isotonic Regression Fit"
  }

  if (ratio) {
    y_lab = "Ratio relative to DMSO"
  } else {
    y_lab = expression(Log[2] ~ "Intensity")
  }

  # Prepare data
  fit_df = data.frame(dose = log10(fit$x), y_pred = fit$y_pred)
  orig_df = data.frame(x = log10(fit$x), y = fit$original_y)

  # Base plot
  model_fit_plot = ggplot2::ggplot() +
    ggplot2::geom_point(data = orig_df,
                        ggplot2::aes(x = x, y = y, color = "Observed Data"), size = 2) +
    ggplot2::geom_line(data = fit_df,
                       ggplot2::aes(x = dose, y = y_pred, color = "Isotonic Regression Fit"), size = 1.2) +
    ggplot2::scale_color_manual(name = NULL,
                                values = c("Observed Data" = "black",
                                           "Isotonic Regression Fit" = "orange")) +
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    ggplot2::ylim(0, NA) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16)
    )

  if (!legend) {
    model_fit_plot = model_fit_plot +
      ggplot2::theme(legend.position = "none")
  }

  # Optional IC50
  if (show_ic50) {
    if (ratio) {
      target_response = 0.5
    } else {
      dmso_idx = which(orig_df$x == -Inf)
      target_response = mean(orig_df$y[dmso_idx], na.rm = TRUE) - 1
    }

    ic50_pred = predictIC50(fit, target_response = target_response)
    y_ic50 = stats::approx(fit$x, fit$y_pred, xout = ic50_pred)$y
    ic50_pred_transform = -log10(ic50_pred) #pIC50

    model_fit_plot = model_fit_plot +
      ggplot2::geom_point(ggplot2::aes(x = log10(ic50_pred), y = y_ic50),
                          shape = 23, size = 3.5, fill = "red", color = "black") +
      ggplot2::annotate("text", x = log10(ic50_pred) + 0.35, y = y_ic50,
                        label = paste("pIC50 =", round(ic50_pred_transform, 2)),
                        vjust = -0.5, hjust = 0.5, size = 4, fontface = "bold")

    # Show IC50 confidence interval if provided
    if (show_ic50 && !is.null(ci)) {
      ci_lower_M = 10^(-ci$ci_lower)
      ci_upper_M = 10^(-ci$ci_upper)

      model_fit_plot = model_fit_plot +
        ggplot2::geom_segment(ggplot2::aes(x = log10(ci_lower_M), xend = log10(ci_lower_M),
                                           y = 0, yend = target_response),
                              linetype = "dotted", color = "gray40", linewidth = 0.8) +
        ggplot2::geom_segment(ggplot2::aes(x = log10(ci_upper_M), xend = log10(ci_upper_M),
                                           y = 0, yend = target_response),
                              linetype = "dotted", color = "gray40", linewidth = 0.8) +
        ggplot2::annotate("rect",
                          xmin = log10(ci_lower_M), xmax = log10(ci_upper_M),
                          ymin = 0, ymax = target_response,
                          alpha = 0.1, fill = "red")
    }
  }

  return(model_fit_plot)
}
