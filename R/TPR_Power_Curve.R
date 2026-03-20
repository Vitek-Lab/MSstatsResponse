# Hardcoded concentration ladders: control (0) and highest dose always present,
# intermediate doses filled from a biologically motivated set.
CONC_MAP <- list(
  "2" = c(0, 3000),
  "3" = c(0, 1000, 3000),
  "4" = c(0, 1, 1000, 3000),
  "5" = c(0, 1, 100, 1000, 3000),
  "6" = c(0, 1, 100, 300, 1000, 3000),
  "7" = c(0, 1, 10, 100, 300, 1000, 3000),
  "8" = c(0, 1, 10, 30, 100, 300, 1000, 3000),
  "9" = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
)

#' Run TPR simulation across a grid of concentration counts and replicate counts
#'
#' Sweeps over combinations of dose counts (2-9) and replicate counts,
#' calling \code{futureExperimentSimulation()} for each combination.
#'
#' @param rep_range Integer vector of length 2, c(min, max) for replicate sweep.
#' @param n_proteins Integer. Number of proteins to simulate. Default: 1000.
#'
#' @return A data.frame with columns: Interaction, TPR, N_rep, NumConcs.
#'
#' @importFrom dplyr filter mutate select if_else
#' @export
run_tpr_simulation <- function(rep_range, n_proteins = 1000) {
  k_grid <- as.integer(names(CONC_MAP))
  rep_grid <- seq(rep_range[1], rep_range[2])

  run_one <- function(n_rep, k_conc, seed = 123) {
    set.seed(seed + n_rep * 100 + k_conc)
    concs_k <- CONC_MAP[[as.character(k_conc)]]

    sim_args <- list(
      N_proteins = n_proteins,
      N_rep = n_rep,
      Concentrations = concs_k,
      IC50_Prediction = FALSE
    )

    temp_res <- futureExperimentSimulation(
      N_proteins = sim_args$N_proteins,
      N_rep = sim_args$N_rep,
      Concentrations = sim_args$Concentrations,
      IC50_Prediction = sim_args$IC50_Prediction
    )
    temp_res$Hit_Rates_Data |>
      dplyr::filter(Category %in% c("TPR (Strong)", "TPR (Weak)")) |>
      dplyr::mutate(
        N_rep = n_rep,
        NumConcs = length(concs_k),
        Interaction = dplyr::if_else(Category == "TPR (Strong)", "Strong", "Weak")
      ) |>
      dplyr::select(Interaction, TPR = Percent, N_rep, NumConcs)
  }

  grid_df <- expand.grid(N_rep = rep_grid, k_conc = k_grid)
  results <- do.call(rbind, lapply(seq_len(nrow(grid_df)), function(i) {
    tryCatch(
      run_one(n_rep = grid_df$N_rep[i], k_conc = grid_df$k_conc[i]),
      error = function(e) {
        warning(paste("Simulation failed for N_rep=", grid_df$N_rep[i],
            ", k_conc=", grid_df$k_conc[i], ":", conditionMessage(e)))
        NULL
      }
    )
  }))

  if (is.null(results) || nrow(results) == 0) {
    stop("All simulation runs failed. Please check your input parameters.")
  }
  return(results)
}

#' Plot TPR power curves for dose response proteomics experiments
#'
#' Creates a two-panel interactive plot showing True Positive Rate
#' for Strong and Weak interactions across concentration and replicate counts.
#'
#' @param simulation_results A data.frame from \code{run_tpr_simulation()}.
#'
#' @return A plotly object with two panels (Strong and Weak interaction).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous
#'   scale_y_continuous scale_linetype_manual labs theme_bw theme element_text
#' @export
plot_tpr_power_curve <- function(simulation_results) {
  k_grid <- sort(unique(simulation_results$NumConcs))
  rep_levels <- sort(unique(simulation_results$N_rep))

  ltypes <- c("dotted", "dotdash", "dashed", "longdash", "solid")
  ltype_values <- setNames(ltypes[seq_along(rep_levels)], as.character(rep_levels))

  make_panel <- function(data, color, show_legend = FALSE) {
    p <- ggplot2::ggplot(data,
      ggplot2::aes(x = NumConcs, y = TPR, linetype = factor(N_rep))) +
      ggplot2::geom_line(linewidth = 1.2, color = color) +
      ggplot2::geom_point(size = 2, color = color) +
      ggplot2::scale_x_continuous(breaks = k_grid) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      ggplot2::scale_linetype_manual(values = ltype_values) +
      ggplot2::labs(
        x = "Number of concentrations",
        y = "True Positive Rate (%)",
        linetype = "Replicates"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )

    if (!show_legend) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    p
  }

  results_strong <- simulation_results[simulation_results$Interaction == "Strong", ]
  results_weak <- simulation_results[simulation_results$Interaction == "Weak", ]

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots. Please install it.")
  }

  p_strong <- plotly::ggplotly(make_panel(results_strong, "#1b9e77", FALSE))
  p_weak <- plotly::ggplotly(make_panel(results_weak, "#d95f02", TRUE))

  plotly::subplot(p_strong, p_weak,
    nrows = 1, shareY = TRUE, titleX = TRUE, titleY = TRUE,
    margin = 0.12
  ) |> plotly::layout(
    margin = list(t = 60),
    annotations = list(
      list(text = "<b>Strong interaction detection power</b>", x = 0.22, y = 1.08,
           xref = "paper", yref = "paper", showarrow = FALSE,
           font = list(size = 12), xanchor = "center"),
      list(text = "<b>Weak interaction detection power</b>", x = 0.78, y = 1.08,
           xref = "paper", yref = "paper", showarrow = FALSE,
           font = list(size = 12), xanchor = "center")
    )
  )
}