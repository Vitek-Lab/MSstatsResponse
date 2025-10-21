#' Test future experimental design using simulated data with user-defined or default templates
#'
#' @param N_proteins Number of proteins to simulate. Default = 300.
#' @param N_rep Number of replicates for each drug concentration. Default = 3.
#' @param N_Control_Rep Number of control replicates. If NULL, uses N_rep.
#' @param Concentrations Numeric vector of drug concentrations (in nM scale).
#'   Default = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000).
#' @param IC50_Prediction Logical. If TRUE, perform IC50 prediction. Default = FALSE.
#' @param data Optional. User's prepared dose-response data (e.g., from MSstatsPrepareDoseResponseFit).
#'   If provided, will extract templates from this data instead of using defaults.
#' @param strong_proteins Character vector of protein IDs to use as strong interaction templates.
#'   Only used if data is provided.
#' @param weak_proteins Character vector of protein IDs to use as weak interaction templates.
#'   Only used if data is provided.
#' @param no_interaction_proteins Character vector of protein IDs to use as no interaction templates.
#'   Only used if data is provided.
#' @param drug_name Character. Name of drug to extract templates for. Default = first non-DMSO drug in data.
#'
#' @return A list containing simulated data, MSstats formatted data, dose-response
#'   fit results, hit rate plots, and optionally IC50 predictions.
#'
#' @examples
#' # Example 1: Quick simulation with default templates (small scale for speed)
#' sim_results <- futureExperimentSimulation(
#'   N_proteins = 50,  # Small number for quick example
#'   N_rep = 2,
#'   N_Control_Rep = 3,
#'   Concentrations = c(0, 10, 100, 1000),  # Fewer doses for speed
#'   IC50_Prediction = FALSE
#' )
#'
#' # View hit rates
#' print(sim_results$Hit_Rates_Data)
#'
#' # Check simulation results
#' print(paste("Simulated", nrow(sim_results$Simulated_Data), "data points"))
#'
#' \dontrun{
#' # Example 2: Full simulation with standard parameters
#' full_sim <- futureExperimentSimulation(
#'   N_proteins = 3000,
#'   N_rep = 3,
#'   N_Control_Rep = 6,
#'   Concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
#'   IC50_Prediction = TRUE
#' )
#'
#' # Display power analysis plot
#' print(full_sim$Hit_Rates_Plot)
#'
#' # Example 3: Using custom templates from your own data
#' # Load and prepare your data
#' data_path <- system.file("extdata", "DIA_MSstats_Normalized.RDS",
#'                          package = "MSstatsResponse")
#' dia_data <- readRDS(data_path)
#'
#' dose_info <- convertGroupToNumericDose(dia_data$ProteinLevelData$GROUP)
#' dia_data$ProteinLevelData$dose <- dose_info$dose_nM * 1e-9
#' dia_data$ProteinLevelData$drug <- dose_info$drug
#'
#' prepared_data <- MSstatsPrepareDoseResponseFit(
#'   dia_data$ProteinLevelData,
#'   dose_column = "dose",
#'   drug_column = "drug",
#'   protein_column = "Protein",
#'   log_abundance_column = "LogIntensities"
#' )
#'
#' # Run simulation with custom templates
#' custom_sim <- futureExperimentSimulation(
#'   N_proteins = 1000,
#'   N_rep = 3,
#'   data = prepared_data,
#'   strong_proteins = c("PROTEIN_A"),
#'   weak_proteins = c("PROTEIN_B"),
#'   no_interaction_proteins = c("PROTEIN_C"),
#'   drug_name = "Drug1",
#'   Concentrations = c(0, 1, 10, 100, 1000, 3000)
#' )
#'}
#'
#' @export
futureExperimentSimulation = function(N_proteins = 300,
                                      N_rep = 3,
                                      N_Control_Rep = NULL,
                                      Concentrations = c(0, 1, 3, 10, 30, 100,
                                                         300, 1000, 3000),
                                      IC50_Prediction = FALSE,
                                      data = NULL,
                                      strong_proteins = NULL,
                                      weak_proteins = NULL,
                                      no_interaction_proteins = NULL,
                                      drug_name = NULL) {

  # Determine whether to use user templates or defaults
  if (!is.null(data)) {
    # Extract templates from user data
    if (is.null(drug_name)) {
      # Use first non-DMSO drug if not specified
      drug_name = unique(data$drug[data$drug != "DMSO"])[1]
      if (is.na(drug_name)) {
        stop("No drug found in data. Please specify drug_name.")
      }
      message(paste("Using drug:", drug_name))
    }

    # Check that at least one protein group is specified
    if (is.null(strong_proteins) && is.null(weak_proteins) && is.null(no_interaction_proteins)) {
      stop("When providing data, you must specify at least one of: strong_proteins, weak_proteins, or no_interaction_proteins")
    }

    # Create template from user-specified proteins
    template = .extractTemplatesFromData(
      data = data,
      strong_proteins = strong_proteins,
      weak_proteins = weak_proteins,
      no_interaction_proteins = no_interaction_proteins,
      drug_name = drug_name,
      concentrations = Concentrations
    )

  } else {
    # Use default templates from package
    template1_path = system.file("extdata", "template1.RDS",
                                 package = "MSstatsResponse")
    template3_path = system.file("extdata", "template3.RDS",
                                 package = "MSstatsResponse")

    if (!file.exists(template1_path)) {
      stop("Template file 1 not found. Please ensure template1.RDS is in inst/extdata/")
    }
    if (!file.exists(template3_path)) {
      stop("Template file 3 not found. Please ensure template3.RDS is in inst/extdata/")
    }

    template = readRDS(template1_path)  # Use template1 as default
    # Note: You could also merge template1 and template3 or choose based on parameters
  }

  # Run simulation
  simulated_data = simulateChemoProteinLevelNonParametric(
    N_proteins = N_proteins,
    template = template,
    rep = N_rep,
    control_rep = N_Control_Rep,
    concentrations = Concentrations,
    outlier_prob = 0
  )

  msstats_format_simulated_data = MSstatsPrepareDoseResponseFit(
    data = simulated_data,
    log_abundance_column = "log_intensity",
    protein_column = "protein",
    drug_column = "drug",
    dose_column = "dose",
    transform_nM_to_M = TRUE
  )

  Interaction_Res = doseResponseFit(msstats_format_simulated_data,
                                    increasing = FALSE,
                                    transform_dose = TRUE,
                                    ratio_response = FALSE)

  Hit_Rates_Plot = plotHitRateMSstatsResponse(
    Interaction_Res,
    rep_count = N_rep,
    concentration_count = length(Concentrations))

  if (IC50_Prediction) {
    IC50_Res = predictIC50(msstats_format_simulated_data,
                           increasing = FALSE,
                           transform_dose = TRUE,
                           ratio_response = FALSE,
                           bootstrap = TRUE,
                           BPPARAM = bpparam())
  } else {
    IC50_Res = NULL
  }

  return(list(
    Simulated_Data = simulated_data,
    MSstats_Simulated_Data = msstats_format_simulated_data,
    DoseResponseFit_Results = Interaction_Res,
    Hit_Rates_Plot = Hit_Rates_Plot$plot,
    Hit_Rates_Data = Hit_Rates_Plot$plot_data,
    PredictIC50_Results = IC50_Res,
    Template_Used = template  # Include template for verification
  ))
}

#' Helper function to extract template profiles from user data
#'
#' @param data Prepared dose-response data
#' @param strong_proteins Vector of protein IDs for strong responders
#' @param weak_proteins Vector of protein IDs for weak responders
#' @param no_interaction_proteins Vector of protein IDs for non-responders
#' @param drug_name Drug to extract templates for
#' @param concentrations Concentrations to include in template (in nM)
#'
#' @return List with template data for each interaction type
.extractTemplatesFromData = function(data,
                                     strong_proteins,
                                     weak_proteins,
                                     no_interaction_proteins,
                                     drug_name,
                                     concentrations) {

  # Filter for specified drug
  drug_data = data %>%
    filter(drug %in% c("DMSO", drug_name))

  # Helper to extract mean profile for a set of proteins
  .getMeanProfile = function(protein_ids, drug_data, concentrations) {
    if (is.null(protein_ids) || length(protein_ids) == 0) {
      # Return flat profile if no proteins specified
      baseline_intensity = 20
      return(data.frame(
        dose = concentrations,
        Intensity = rep(2^baseline_intensity, length(concentrations)),
        LogIntensities = rep(baseline_intensity, length(concentrations)),
        ratio = rep(1, length(concentrations))
      ))
    }

    # Get mean profile across proteins
    profile = drug_data %>%
      filter(protein %in% protein_ids) %>%
      mutate(dose_nM = round(dose * 1e9)) %>%  # Convert M to nM and round to handle precision
      group_by(dose_nM) %>%
      summarise(
        LogIntensities = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      rename(dose = dose_nM)

    # Create complete profile with all concentrations
    complete_profile = data.frame(dose = concentrations)

    # Merge with existing data
    complete_profile = complete_profile %>%
      left_join(profile, by = "dose")

    # Fill missing values with interpolation or nearest neighbor
    for (i in which(is.na(complete_profile$LogIntensities))) {
      curr_dose = complete_profile$dose[i]

      if (curr_dose == 0 && any(!is.na(complete_profile$LogIntensities))) {
        # If DMSO is missing, use the highest observed value (least affected)
        complete_profile$LogIntensities[i] = max(complete_profile$LogIntensities, na.rm = TRUE)
      } else if (any(!is.na(complete_profile$LogIntensities))) {
        # Find nearest dose with data on log scale
        doses_with_data = complete_profile$dose[!is.na(complete_profile$LogIntensities)]
        log_distances = abs(log10(curr_dose + 1) - log10(doses_with_data + 1))
        nearest_idx = which.min(log_distances)
        nearest_dose = doses_with_data[nearest_idx]
        nearest_value = complete_profile$LogIntensities[complete_profile$dose == nearest_dose]

        # Apply simple interpolation based on direction
        if (length(doses_with_data) >= 2) {
          # Linear interpolation on log scale if we have enough points
          lower_doses = doses_with_data[doses_with_data < curr_dose]
          upper_doses = doses_with_data[doses_with_data > curr_dose]

          if (length(lower_doses) > 0 && length(upper_doses) > 0) {
            # Interpolate between nearest lower and upper
            lower_dose = max(lower_doses)
            upper_dose = min(upper_doses)
            lower_val = complete_profile$LogIntensities[complete_profile$dose == lower_dose]
            upper_val = complete_profile$LogIntensities[complete_profile$dose == upper_dose]

            # Linear interpolation on log-dose scale
            log_lower = log10(lower_dose + 1)
            log_upper = log10(upper_dose + 1)
            log_curr = log10(curr_dose + 1)

            weight = (log_curr - log_lower) / (log_upper - log_lower)
            complete_profile$LogIntensities[i] = lower_val * (1 - weight) + upper_val * weight
          } else {
            # Use nearest neighbor
            complete_profile$LogIntensities[i] = nearest_value
          }
        } else {
          # Only one point available
          complete_profile$LogIntensities[i] = nearest_value
        }
      } else {
        # No data at all - use baseline
        complete_profile$LogIntensities[i] = 20
      }
    }

    # Calculate Intensity and ratio
    complete_profile = complete_profile %>%
      mutate(Intensity = 2^LogIntensities)

    # Calculate ratio relative to DMSO
    dmso_intensity = complete_profile$Intensity[complete_profile$dose == 0]
    if (length(dmso_intensity) == 0 || is.na(dmso_intensity)) {
      dmso_intensity = complete_profile$Intensity[1]
    }

    complete_profile = complete_profile %>%
      mutate(ratio = Intensity / dmso_intensity) %>%
      select(dose, Intensity, LogIntensities, ratio)

    return(complete_profile)
  }

  # Extract templates for each category
  template = list(
    strong_interaction = .getMeanProfile(strong_proteins, drug_data, concentrations),
    weak_interaction = .getMeanProfile(weak_proteins, drug_data, concentrations),
    no_interaction = .getMeanProfile(no_interaction_proteins, drug_data, concentrations)
  )

  # Validate template format
  for (category in names(template)) {
    if (!all(c("dose", "Intensity", "LogIntensities", "ratio") %in% names(template[[category]]))) {
      stop(paste("Template for", category, "missing required columns"))
    }
    if (nrow(template[[category]]) != length(concentrations)) {
      stop(paste("Template for", category, "doesn't match number of concentrations"))
    }
  }

  message("Templates successfully extracted from data")
  message(paste("  Strong interaction proteins:", paste(strong_proteins, collapse = ", ")))
  message(paste("  Weak interaction proteins:", paste(weak_proteins, collapse = ", ")))
  message(paste("  No interaction proteins:", paste(no_interaction_proteins, collapse = ", ")))

  return(template)
}

#' Simulate chemoproteomics data at the protein level - non-parametric approach
#'
#' @param N_proteins Number of proteins in simulation. Default = 3000.
#' @param TP Proportion of strong interacting proteins. Default = 0.333.
#' @param TW Proportion of weak interacting proteins. Default = 0.333.
#' @param TN Proportion of non-interacting proteins. Default = 0.333.
#' @param concentrations Numeric vector of drug concentrations in simulation experiment.
#'   Default = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000).
#' @param rep Number of replicates for each drug concentration. Default = 3.
#' @param seed Simulation seed for reproducibility. Default = 3.
#' @param var_tech Combined technical and biological variation. Default = 0.4.
#' @param control_rep Number of control replicates. If NULL, uses rep.
#' @param template List containing real protein level data representing different
#'   levels of interactions.
#' @param outlier_prob Probability of sample outlier. Default = 0.05.
#'
#' @return A data.frame with simulated chemoproteomics data.
#' @importFrom stats rnorm rlnorm

simulateChemoProteinLevelNonParametric = function(N_proteins = 3000,
                                                     TP = 0.333,
                                                     TW = 0.333,
                                                     TN = 0.333,
                                                     concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
                                                     rep = 3,
                                                     seed = NULL,
                                                     var_tech = 0.4,
                                                     control_rep = NULL,
                                                     template = list(
                                                       strong_interaction = data.frame(dose = c(), LogIntensities = c()),
                                                       weak_interaction = data.frame(dose = c(), LogIntensities = c()),
                                                       no_interaction = data.frame(dose = c(), LogIntensities = c())
                                                     ),
                                                     outlier_prob = 0.05) {

  if (is.null(control_rep)) {
    control_rep = rep
  }

  # Validate concentrations
  valid_concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
  if (!all(concentrations %in% valid_concentrations)) {
    stop("Error: `concentrations` must be a subset of c(0, 1, 3, 10, 30, 100, 300, 1000, 3000).")
  }

  # Subset template data frames by concentrations
  subset_template = function(data, concentrations) {
    data[data$dose %in% concentrations, ]
  }
  template = lapply(template, subset_template, concentrations = concentrations)

  # Set random seed for reproducibility
  #if (!is.null(seed)) {
  #  set.seed(seed)
 # }

  # Calculate the number of proteins for each interaction type
  tp = ceiling(N_proteins * TP)
  tw = ceiling(N_proteins * TW)
  tn = ceiling(N_proteins * TN)

  # Function to simulate data for a protein group with replicates
  simulate_proteins = function(num_proteins, interaction_type, var_tech, template_data, rep, control_rep) {
    if (num_proteins == 0) {
      return(data.frame())
    }
    data_list = vector("list", num_proteins)
    for (i in seq_len(num_proteins)) {
      # Determine replicates for each dose
      dose_vector = template_data$dose
      replicate_counts = ifelse(dose_vector == 0, control_rep, rep)
      df_list = vector("list", length(dose_vector))

      for (j in seq_along(dose_vector)) {
        dose_value = dose_vector[j]
        n_reps = replicate_counts[j]

        df_tmp = data.frame(
          dose = dose_value,
          replicate = seq_len(n_reps),
          log_intensity = rep(template_data$LogIntensities[j], n_reps) + rnorm(n_reps, 0, var_tech)
        )
        df_list[[j]] = df_tmp
      }

      df = do.call(rbind, df_list)
      df$protein = paste0("p_", interaction_type, "_", i)
      df$protein_group = interaction_type
      data_list[[i]] = df
    }
    do.call(rbind, data_list)
  }

  # Simulate data for each protein group
  SimData_TP = simulate_proteins(tp, "strong_interaction", var_tech, template$strong_interaction, rep, control_rep)
  SimData_TW = simulate_proteins(tw, "weak_interaction", var_tech, template$weak_interaction, rep, control_rep)
  SimData_TN = simulate_proteins(tn, "no_interaction", var_tech, template$no_interaction, rep, control_rep)

  # Combine data from all protein groups
  intensity_mu = rbind(SimData_TP, SimData_TW, SimData_TN)

  # Add DMSO labels for dose = 0
  intensity_mu$sample_group = ifelse(intensity_mu$dose == 0, "DMSO", "treated")

  # Introduce outliers and mark them
  add_outliers = function(data, outlier_prob) {
    num_samples = nrow(data)
    num_outliers = ceiling(outlier_prob * num_samples)

    # Sample outlier indices
    outlier_indices = sample(seq_len(num_samples), num_outliers, replace = FALSE)

    # Ensure at least one outlier in strong and weak interaction groups
    strong_idx = which(data$protein_group == "strong_interaction")
    weak_idx = which(data$protein_group == "weak_interaction")
    if (length(strong_idx) > 0) {
      outlier_indices = unique(c(outlier_indices, sample(strong_idx, 1)))
    }
    if (length(weak_idx) > 0) {
      outlier_indices = unique(c(outlier_indices, sample(weak_idx, 1)))
    }

    # Add outliers using a log-normal distribution for realistic variance
    data$outlier = FALSE
    lognorm_offset = rlnorm(length(outlier_indices), meanlog = 0, sdlog = 1)
    sign_flip = sample(c(-1, 1), length(outlier_indices), replace = TRUE)

    data$log_intensity[outlier_indices] = data$log_intensity[outlier_indices] +
      (sign_flip * lognorm_offset)

    data$outlier[outlier_indices] = TRUE

    # Add simulation as drug
    data$drug = 'simulation'
    return(data)
  }

  # Apply function to introduce outliers
  intensity_mu = add_outliers(intensity_mu, outlier_prob)

  return(intensity_mu)
}

#' Plot hit rates by category
#'
#' @param results Output of interaction test from doseResponseFit()
#' @param rep_count Number of replicates per concentration in simulation
#' @param concentration_count Number of concentrations in simulation
#'
#' @return A list containing the plot and plot data
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs scale_fill_manual
#'   scale_y_continuous theme_classic theme element_text
#' @import dplyr
plotHitRateMSstatsResponse = function(results, rep_count, concentration_count) {

  dt = data.frame(results)

  # Use correct column for filtering
  filtered_dt = dt %>%
    dplyr::filter(adjust_pval < 0.05)

  # Identify protein groups
  strong_proteins = dt %>%
    dplyr::filter(grepl("strong", protein)) %>%
    dplyr::pull(protein)

  weak_proteins = dt %>%
    dplyr::filter(grepl("weak", protein)) %>%
    dplyr::pull(protein)

  no_proteins = dt %>%
    dplyr::filter(grepl("no", protein)) %>%
    dplyr::pull(protein)

  # Function to calculate percentage passing per group
  calculate_percentage = function(group_proteins) {
    total = length(group_proteins)
    if (total == 0) return(NA_real_)
    passed = sum(group_proteins %in% filtered_dt$protein)
    round((passed / total) * 100, 1)
  }

  # Compute percentages
  percent_strong = calculate_percentage(strong_proteins)
  percent_weak = calculate_percentage(weak_proteins)
  percent_no = calculate_percentage(no_proteins)

  # Prepare data table for plotting
  plot_data = data.frame(
    Category = c("TPR (Strong)", "TPR (Weak)", "FPR"),
    Percent = c(percent_strong, percent_weak, percent_no)
  )

  # Plot
  p = ggplot2::ggplot(plot_data, ggplot2::aes(x = Category, y = Percent, fill = Category)) +
    ggplot2::geom_bar(stat = "identity", color = "black", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(Percent, 1), "%")),
                       vjust = -0.3, size = 6, fontface = "bold") +
    ggplot2::labs(title = "Future Experiment Performance",
                  subtitle = paste("with", concentration_count, "doses and", rep_count, "replicates per dose"),
                  x = "Performance Metrics",
                  y = "Rate (%)") +
    ggplot2::scale_fill_manual(values = c("TPR (Strong)" = "#1b9e77",
                                          "TPR (Weak)" = "#d95f02",
                                          "FPR" = "#7570b3")) +
    ggplot2::scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 100, 20)) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(color = "black", size = 16),
      axis.title = ggplot2::element_text(face = "bold", size = 18),
      axis.ticks = ggplot2::element_line(color = "black"),
      legend.position = "none"
    )

  return(list(plot = p, plot_data = plot_data))
}

