#' Plot isotonic regression fit with optional IC50 for a single protein and drug
#'
#' @param data Protein-level dataset (e.g., output of MSstatsPrepareDoseResponseFit).
#' @param protein_name Character. Protein name to plot.
#' @param drug_name Character. Drug name to plot.
#' @param ratio_response Logical. If TRUE, compute IC50 on ratio scale; if FALSE, use log2 intensities.
#' @param transform_dose Logical. If TRUE, applies log10(dose + 1). Default is TRUE.
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param ... Additional arguments passed to internal `plot_isotonic()` function.
#'
#' @return A ggplot object.
#' @export
VisualizeResponseProtein = function(data,
                                    protein_name,
                                    drug_name,
                                    ratio_response = TRUE,
                                    transform_dose = TRUE,
                                    show_ic50 = TRUE,
                                    ...) {
  # Filter data for the selected protein and drug + DMSO
  data_single = data %>%
    dplyr::filter(protein == protein_name & drug %in% c("DMSO", drug_name))

  x_org = data_single$dose
  y_org = data_single$response

  order_idx = order(x_org)
  x = x_org[order_idx]
  y = y_org[order_idx]

  # Only use plotting args in plot_isotonic, not here
  fit = fit_isotonic_regression(
    x = x, y = y,
    increasing = FALSE,
    transform_x = transform_dose,
    ratio_y = ratio_response,
    test_significance = FALSE
  )

  # Now pass the extra args into the plotting function
  p = plot_isotonic(
    fit = fit,
    ratio = ratio_response,
    show_ic50 = show_ic50,
    drug_name = drug_name,
    protein_name = protein_name,
    ...
  )

  return(p)
}




#' Plot Isotonic Regression Model
#'
#' @param fit A model object returned by fit_isotonic_regression().
#' @param show_ic50 Logical. If TRUE, adds vertical line and annotation for IC50.
#' @param x_lab Label for x-axis.
#' @param y_lab Label for y-axis.
#' @param title Title for the plot.
#' @param legend Logical. Show legend if TRUE.
#' @param theme_style ggplot2 theme name to apply (default = "classic").
#' @param original_label Logical. If TRUE, replace x-axis tick labels with original dose labels.
#'
#' @return A ggplot object.
#'
#' @importFrom stats pf
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment annotate
#' @importFrom ggplot2 scale_color_manual scale_x_continuous labs ylim theme_minimal
#' @importFrom ggplot2 theme element_rect element_blank element_text

plot_isotonic = function(fit,
                         ratio = TRUE,
                         show_ic50 = FALSE,
                         drug_name = NULL,
                         protein_name = NULL,
                         x_lab = "dose (nM)",
                         y_lab = "log2 Intensity",
                         title = NULL,
                         legend = FALSE,
                         theme_style = "classic",
                         original_label = TRUE) {
  library(ggplot2)

  # Construct title if not provided
  if (is.null(title) && !is.null(drug_name) && !is.null(protein_name)) {
    title = paste0(drug_name, ": ", protein_name, " Response")
  } else if (is.null(title)) {
    title = "Isotonic Regression Fit"
  }

  # Prepare data
  fit_df = data.frame(dose = fit$x, y_pred = fit$y_pred)
  orig_df = data.frame(x = fit$x, y = fit$original_y)

  dose_list = log10(c(0, 1, 3, 10, 30, 100, 300, 1000, 3000) + 1)

  # Base plot
  model_fit_plot = ggplot() +
    geom_point(data = orig_df,
               aes(x = x, y = y, color = "Observed Data"), size = 2) +
    geom_line(data = fit_df,
              aes(x = dose, y = y_pred, color = "Isotonic Regression Fit"), size = 1.2) +
    scale_color_manual(name = NULL,
                       values = c("Observed Data" = "black",
                                  "Isotonic Regression Fit" = "orange")) +
    labs(x = x_lab, y = y_lab, title = title) +
    ylim(0, NA) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_blank(),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

  if (!legend) {
    model_fit_plot = model_fit_plot +
      theme(legend.position = "none")
  }

  if (original_label) {
    model_fit_plot = model_fit_plot +
      scale_x_continuous(
        breaks = dose_list,
        labels = c("DMSO", "1", "3", "10", "30", "100", "300", "1000", "3000")
      )
  }

  # Optional IC50
  if (show_ic50) {
    if (ratio) {
      target_response = 0.5
    } else {
      dmso_idx = which(orig_df$x == 0)
      target_response = mean(orig_df$y[dmso_idx], na.rm = TRUE) - 1
    }

    ic50_pred = predict_ic50(fit, target_response = target_response)
    y_ic50 = approx(fit$x, fit$y_pred, xout = ic50_pred)$y
    ic50_pred_transform = 10^ic50_pred

    model_fit_plot = model_fit_plot +
      geom_point(aes(x = ic50_pred, y = y_ic50),
                 shape = 23, size = 3.5, fill = "red", color = "black") +
      annotate("text", x = ic50_pred + 0.35, y = y_ic50,
               label = paste("IC50 =", round(ic50_pred_transform, 2)),
               vjust = -0.5, hjust = 0.5, size = 4, fontface = "bold")
  }

  return(model_fit_plot)
}

