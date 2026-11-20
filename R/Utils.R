# Suppress R CMD check NOTE about undefined global variables used in NSE
utils::globalVariables(c(
  "drug", "protein", "dose", "x", "y", "y_pred"
))

# Suppress R CMD check notes about undefined globals used in NSE
utils::globalVariables(c(
  # Variables from .extractTemplatesFromData
  "dose_nM", "response", "LogIntensities", "Intensity", "ratio",

  # Variables from PlotExperimentPowerCurve
  "Category", "Interaction", "Percent", "N_rep", "NumConcs", "TPR",

  # Variables from plotHitRateMSstatsResponse
  "adj.pvalue",

  # Variables from calculateTurnoverRatios (data.table NSE)
  ".SD", "Protein", "BaseSequence", "Label", "TimeVal", "Run",
  "Heavy", "Light", "Total", "H_frac", "L_frac", "tracer_factor",

  # Variables from the other protein turnover helpers (dplyr NSE)
  "PROTEIN", "IC50", "k_obs", "coverage_per_peptide", "p_protein",
  "peptide_median_light", "peptide_rank", "coverage_score",
  "light_intensity_score", "monotonicity_score", "validity_flag",
  "n_light_peptides", "observed_cells", "n_max_possible",
  "n_heavy_peptides", "mean_weight", "pep_factor", "confidence",

  # Variables from plotIsotonic and .getMeanProfile
  "weight", "dose_str",

  # Add any other NSE variables used in your package
  "drug", "protein", "dose", "x", "y", "y_pred"
))
# Utility functions (if needed in future)

# custom step interpolation (left-continuous)
step_interpolate = function(x, y, new_x) {
  order_idx = order(x)
  x = x[order_idx]
  y = y[order_idx]
  sapply(new_x, function(xi) {
    idx = max(which(x <= xi), na.rm = TRUE)
    return(y[idx])
  })
}
