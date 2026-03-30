# ============================================================================
# TPR Power Curve - Simulation and Visualization
# ============================================================================

#' Build concentration subsets for dose-response power analysis
#'
#' Given a user's full set of experimental concentrations, generates subsets
#' of size 2 through N by greedily selecting doses that are farthest from
#' the log-median of the already-selected set. Control (0) and the highest
#' dose are always included.
#'
#' @param concentrations Numeric vector of all available concentrations.
#' @param dose_range Integer vector of length 2, c(min, max) specifying the
#'   range of subset sizes to generate.
#'
#' @return A named list where each element is a numeric vector of concentrations.
#'   Names are the subset sizes (e.g., "2", "3", ...).
#' @noRd
.build_concentration_ladders <- function(concentrations, dose_range) {
  concentrations <- sort(unique(concentrations))
  if (length(concentrations) < 2) {
    stop("At least 2 unique concentrations are required.")
  }
  if (!(0 %in% concentrations)) {
    stop("Concentrations must include 0 (control).")
  }

  control <- min(concentrations)
  highest <- max(concentrations)
  candidates <- setdiff(concentrations, c(control, highest))

  # Greedily build up from {control, highest} by picking the dose
  # farthest from log(median) of the current set at each step
  selected <- c(control, highest)
  selection_order <- list()
  selection_order[["2"]] <- sort(selected)

  remaining <- candidates
  max_k <- min(dose_range[2], length(concentrations))

  for (k in 3:max_k) {
    if (length(remaining) == 0) break

    log_selected <- log10(selected + 1)
    log_median <- median(log_selected)
    log_remaining <- log10(remaining + 1)

    distances <- abs(log_remaining - log_median)
    best_idx <- which.max(distances)
    pick <- remaining[best_idx]

    selected <- c(selected, pick)
    remaining <- remaining[-best_idx]
    selection_order[[as.character(k)]] <- sort(selected)
  }

  # Filter to requested dose_range
  k_range <- seq(dose_range[1], min(dose_range[2], length(concentrations)))
  result <- selection_order[as.character(k_range)]
  result <- result[!sapply(result, is.null)]

  if (length(result) == 0) {
    stop("No valid concentration subsets could be generated for the given dose_range.")
  }
  return(result)
}

#' Run a single TPR simulation for one replicate/concentration combination
#'
#' @param n_rep Integer. Number of replicates per dose.
#' @param concs Numeric vector. Concentrations to use.
#' @param n_proteins Integer. Number of proteins to simulate.
#' @param data Optional prepared dose-response data for user-based templates.
#' @param protein Optional protein ID for strong interaction template.
#' @param seed Integer. Base seed for reproducibility.
#'
#' @return A data.frame with columns: Interaction, TPR, N_rep, NumConcs.
#' @noRd
.run_one_tpr_simulation <- function(n_rep, concs, n_proteins, data = NULL,
                                     protein = NULL, seed = 123) {
  set.seed(seed + n_rep * 100 + length(concs))

  sim_args <- list(
    N_proteins = n_proteins,
    N_rep = n_rep,
    Concentrations = concs,
    IC50_Prediction = FALSE
  )
  if (!is.null(data)) {
    sim_args$data <- data
  }
  if (!is.null(protein) && !is.null(data)) {
    sim_args$strong_proteins <- protein
  }

  temp_res <- do.call(futureExperimentSimulation, sim_args)
  temp_res$Hit_Rates_Data |>
    dplyr::filter(Category %in% c("TPR (Strong)", "TPR (Weak)")) |>
    dplyr::mutate(
      N_rep = n_rep,
      NumConcs = length(concs),
      Interaction = dplyr::if_else(Category == "TPR (Strong)", "Strong", "Weak")
    ) |>
    dplyr::select(Interaction, TPR = Percent, N_rep, NumConcs)
}

#' Simulate detection power across experimental design configurations
#'
#' Evaluates how well a dose-response proteomics experiment can detect
#' drug-protein interactions under different experimental designs. For each
#' combination of replicate count and number of doses, the function simulates
#' a full experiment and measures the true positive rate (TPR); the
#' percentage of known interacting proteins that are correctly identified.
#'
#' This helps experimentalists answer practical questions like: "If I can
#' only afford 3 replicates per dose, how many dose levels do I need to
#' reliably detect weak interactions?"
#'
#' Concentration subsets are selected automatically from the user's available
#' doses: control (0) and the highest dose are always included, and
#' intermediate doses are added by choosing the one farthest from the
#' log-median of the already-selected set.
#'
#' @param rep_range Integer vector of length 2, \code{c(min, max)}. The range
#'   of biological replicates per dose to evaluate. Each value in the sequence
#'   \code{min:max} produces one line in the resulting power curve.
#' @param concentrations Numeric vector. The full set of drug concentrations
#'   available in the experiment (e.g., \code{c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)}).
#' @param dose_range Integer vector of length 2, \code{c(min, max)}. The range
#'   of dose counts to sweep. For example, \code{c(2, 9)} evaluates designs
#'   with 2 doses, 3 doses, ..., up to 9 doses.
#' @param data Optional. User's prepared dose-response data (e.g., from
#'   \code{MSstatsPrepareDoseResponseFit}). If provided, protein-specific
#'   templates are extracted from this data instead of using defaults.
#' @param protein Optional. Character string specifying a protein ID to use as
#'   the strong interaction template. Only used when \code{data} is provided.
#' @param n_proteins Integer. Number of proteins to simulate per run. Larger
#'   values give more stable estimates but increase runtime. Default: 1000.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{Interaction}{Character. "Strong" or "Weak".}
#'     \item{TPR}{Numeric. True positive rate as a percentage (0-100).}
#'     \item{N_rep}{Integer. Number of replicates used in this simulation.}
#'     \item{NumConcs}{Integer. Number of concentrations used in this simulation.}
#'   }
#'
#' @examples
#' # Quick example with small simulation (for speed)
#' results <- run_tpr_simulation(
#'   rep_range = c(1, 2),
#'   concentrations = c(0, 10, 100, 1000),
#'   dose_range = c(2, 4),
#'   n_proteins = 50
#' )
#' head(results)
#'
#' \dontrun{
#' # Full power analysis with standard chemoproteomic doses
#' results <- run_tpr_simulation(
#'   rep_range = c(1, 5),
#'   concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
#'   dose_range = c(2, 9),
#'   n_proteins = 1000
#' )
#'
#' # Visualize results
#' plot_tpr_power_curve(results)
#' }
#'
#' @importFrom dplyr filter mutate select if_else
#' @export
run_tpr_simulation <- function(rep_range, concentrations, dose_range,
                                data = NULL, protein = NULL, n_proteins = 1000) {
  if (!is.numeric(rep_range) || length(rep_range) != 2 || rep_range[1] > rep_range[2]) {
    stop("rep_range must be a numeric vector of length 2 with c(min, max) where min <= max.")
  }
  if (!is.numeric(dose_range) || length(dose_range) != 2 || dose_range[1] > dose_range[2]) {
    stop("dose_range must be a numeric vector of length 2 with c(min, max) where min <= max.")
  }
  if (dose_range[1] < 2) {
    stop("dose_range minimum must be at least 2 (control + one treatment).")
  }
  if (!is.null(data) && is.null(protein)) {
    stop("protein must be specified when data is provided.")
  }

  conc_subsets <- .build_concentration_ladders(concentrations, dose_range)
  k_grid <- as.integer(names(conc_subsets))
  rep_grid <- seq(rep_range[1], rep_range[2])

  grid_df <- expand.grid(N_rep = rep_grid, k_conc = k_grid)
  results <- do.call(rbind, lapply(seq_len(nrow(grid_df)), function(i) {
    k <- as.character(grid_df$k_conc[i])
    concs <- conc_subsets[[k]]
    tryCatch(
      .run_one_tpr_simulation(
        n_rep = grid_df$N_rep[i],
        concs = concs,
        n_proteins = n_proteins,
        data = data,
        protein = protein
      ),
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

#' Create a single TPR panel plot with color gradient by replicate count
#'
#' @param data Data frame with TPR results for one interaction type.
#' @param color_low Character. Light shade for fewest replicates.
#' @param color_high Character. Dark shade for most replicates.
#' @param k_grid Integer vector. X-axis breaks.
#' @param show_legend Logical. Whether to display the legend.
#'
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous
#'   scale_y_continuous scale_color_manual labs theme_bw theme element_text
#' @noRd
.make_tpr_panel <- function(data, color_low, color_high, k_grid, show_legend = FALSE) {
  rep_levels <- sort(unique(data$N_rep))
  n_levels <- length(rep_levels)

  # Generate gradient from light to dark
  color_palette <- grDevices::colorRampPalette(c(color_low, color_high))(n_levels)
  names(color_palette) <- as.character(rep_levels)

  p <- ggplot2::ggplot(data,
    ggplot2::aes(x = NumConcs, y = TPR,
                 color = factor(N_rep), group = factor(N_rep))) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = k_grid) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(
      x = "Number of concentrations",
      y = "True Positive Rate (%)",
      color = "Replicates"
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

#' Visualize detection power across experimental designs
#'
#' Creates an interactive plot showing how the true positive rate
#' (detection power) changes with the number of doses and replicates.
#' Only the results for the user-selected protein template (passed as
#' the strong interaction category) are displayed. Line color shading goes
#' from light (fewest replicates) to dark (most replicates).
#'
#' @param simulation_results A \code{data.frame} returned by
#'   \code{\link{run_tpr_simulation}}.
#'
#' @return An interactive \code{plotly} object.
#'
#' @examples
#' \dontrun{
#' results <- run_tpr_simulation(
#'   rep_range = c(1, 3),
#'   concentrations = c(0, 10, 100, 1000),
#'   dose_range = c(2, 4),
#'   n_proteins = 300
#' )
#' plot_tpr_power_curve(results)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous
#'   scale_y_continuous scale_color_manual labs theme_bw theme element_text
#' @importFrom plotly ggplotly layout
#' @importFrom grDevices colorRampPalette
#' @export
plot_tpr_power_curve <- function(simulation_results) {
  k_grid <- sort(unique(simulation_results$NumConcs))

  # Use only the Strong interaction results (user-selected protein template)
  results_protein <- simulation_results[simulation_results$Interaction == "Strong", ]

  p <- .make_tpr_panel(results_protein, "#a6dba0", "#1b7837", k_grid, TRUE)

  plotly::ggplotly(p) |> plotly::layout(
    margin = list(t = 60),
    annotations = list(
      list(text = "<b>Interaction detection power</b>", x = 0.5, y = 1.08,
           xref = "paper", yref = "paper", showarrow = FALSE,
           font = list(size = 14), xanchor = "center")
    )
  )
}