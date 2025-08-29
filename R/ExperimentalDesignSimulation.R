#' Test future experimental design using simulated data
#'
#' @param concentrations list of drug concentrations (in nM scale) in sim experiment
#' @param rep number of replicates for each drug concentration
#' @param template real protein level data representing different level of interactions
#' @param outlier_prob probability of sample outlier
#'
#'

FutureExperimentSimulation = function(N_proteins = 300,
                                      N_rep = 3,
                                      N_Control_Rep = NULL,
                                      Concentrations = c(0, 1, 3, 10, 30, 100,
                                                         300, 1000, 3000),
                                      IC50_Prediction = FALSE
                                      ){
  # Load template data ---
  template1 = readRDS("data/template1_SRC_TEC_EIF2AK4.RDS")
  template3 = readRDS("data/template3_Q9H0K1_P09693.RDS")

  simulated_data = simulate_chemo_proteinlevel_nonparametric(
    N_proteins = N_proteins,
    template = template1,
    rep = N_rep,
    control_rep = N_Control_Rep,
    concentrations = Concentrations,
    outlier_prob = 0
  )

  msstats_format_simulated_data = MSstatsPrepareDoseResponseFit(
    data = simulated_data,
    log_abundance_column = "log_intensity",
    protein_column = "protein",
    transform_nM_to_M = TRUE
    )

  Interaction_Res = DoseResponseFit(msstats_format_simulated_data,
                  increasing = FALSE,
                  transform_dose = TRUE,
                  ratio_response = FALSE)

  Hit_Rates_Plot = plot_hit_rate_MSstatsResponse(
    Interaction_Res,
    rep_count = N_rep,
    concentration_count = length(Concentrations))

  if(IC50_Prediction){
   IC50_Res = PredictIC50(msstats_format_simulated_data,
                             increasing = FALSE,
                             transform_dose = TRUE,
                             ratio_response = FALSE,
                             bootstrap = TRUE,
                             numberOfCores = 4)
   } else {IC50_Res = NULL}

  return(list(
    Simulated_Data = simulated_data,
    MSstats_Simulated_Data = msstats_format_simulated_data,
    DoseResponseFit_Results = Interaction_Res,
    Hit_Rates_Plot = Hit_Rates_Plot$plot,
    Hit_Rates_Data = Hit_Rates_Plot$plot_data,
    PredictIC50_Results =   IC50_Res))

}



#' Simulate chemoproteomics data at the protein level - non-parametric approach
#'
#' @param N_proteins number of proteins in sim
#' @param TP percent of strong interacting proteins in sim
#' @param TW percent of weak interacting proteins in sim
#' @param TN percent of non interacting proteins in sim
#' @param concentrations list of drug concentrations in sim experiment
#' @param rep number of replicates for each drug concentration
#' @param seed simulation seed
#' @param var_teach combined technical and biological variation
#' @param template real protein level data representing different level of interactions
#' @param outlier_prob probability of sample outlier
#'


simulate_chemo_proteinlevel_nonparametric = function(N_proteins = 3000,
                                                      TP = 0.333,
                                                      TW = 0.333,
                                                      TN = 0.333,
                                                      concentrations = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000),
                                                      rep = 3,
                                                      seed = 3,
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
  valid_concentrations <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
  if (!all(concentrations %in% valid_concentrations)) {
    stop("Error: `concentrations` must be a subset of c(0, 1, 3, 10, 30, 100, 300, 1000, 3000).")
  }

  # Subset template data frames by concentrations
  subset_template <- function(data, concentrations) {
    data[data$dose %in% concentrations, ]
  }
  template <- lapply(template, subset_template, concentrations = concentrations)

  # Set random seed for reproducibility
  set.seed(seed)

  # Calculate the number of proteins for each interaction type
  tp <- ceiling(N_proteins * TP)
  tw <- ceiling(N_proteins * TW)
  tn <- ceiling(N_proteins * TN)

  # Function to simulate data for a protein group with replicates
  simulate_proteins <- function(num_proteins, interaction_type, var_tech, template_data, rep, control_rep) {
    if (num_proteins == 0) {
      return(data.frame())
    }
    data_list <- vector("list", num_proteins)
    for (i in seq_len(num_proteins)) {
      # Determine replicates for each dose
      dose_vector <- template_data$dose
      replicate_counts <- ifelse(dose_vector == 0, control_rep, rep)
      df_list <- vector("list", length(dose_vector))

      for (j in seq_along(dose_vector)) {
        dose_value <- dose_vector[j]
        n_reps <- replicate_counts[j]

        df_tmp <- data.frame(
          dose = dose_value,
          replicate = seq_len(n_reps),
          log_intensity = rep(template_data$LogIntensities[j], n_reps) + rnorm(n_reps, 0, var_tech)
        )
        df_list[[j]] <- df_tmp
      }

      df <- do.call(rbind, df_list)
      df$protein <- paste0("p_", interaction_type, "_", i)
      df$protein_group <- interaction_type
      data_list[[i]] <- df
    }
    do.call(rbind, data_list)
  }


  # Simulate data for each protein group
  SimData_TP <- simulate_proteins(tp, "strong_interaction", var_tech, template$strong_interaction, rep, control_rep)
  SimData_TW <- simulate_proteins(tw, "weak_interaction", var_tech, template$weak_interaction, rep, control_rep)
  SimData_TN <- simulate_proteins(tn, "no_interaction", var_tech, template$no_interaction, rep, control_rep)

  # Combine data from all protein groups
  intensity_mu <- rbind(SimData_TP, SimData_TW, SimData_TN)

  # Add DMSO labels for dose = 0
  intensity_mu$sample_group <- ifelse(intensity_mu$dose == 0, "DMSO", "treated")

  # Introduce outliers and mark them
  add_outliers <- function(data, outlier_prob) {
    num_samples <- nrow(data)
    num_outliers <- ceiling(outlier_prob * num_samples)  # 1-2 outliers per protein

    # Sample outlier indices
    outlier_indices <- sample(seq_len(num_samples), num_outliers, replace = FALSE)

    # Ensure at least one outlier in strong and weak interaction groups
    strong_idx <- which(data$protein_group == "strong_interaction")
    weak_idx <- which(data$protein_group == "weak_interaction")
    if (length(strong_idx) > 0) {
      outlier_indices <- unique(c(outlier_indices, sample(strong_idx, 1)))
    }
    if (length(weak_idx) > 0) {
      outlier_indices <- unique(c(outlier_indices, sample(weak_idx, 1)))
    }

    # Add outliers using a log-normal distribution for realistic variance
    data$outlier <- FALSE
    lognorm_offset <- rlnorm(length(outlier_indices), meanlog = 0, sdlog = 1) # Log-normal deviation
    sign_flip <- sample(c(-1, 1), length(outlier_indices), replace = TRUE) # Random direction (up/down)

    data$log_intensity[outlier_indices] <- data$log_intensity[outlier_indices] +
      (sign_flip * lognorm_offset)

    data$outlier[outlier_indices] <- TRUE

    # add simulation as drug
    data$drug = 'simulation'
    return(data)
  }

  # Apply function to introduce outliers
  intensity_mu <- add_outliers(intensity_mu, outlier_prob)

  return(intensity_mu)
}



#' Plot hit rates by category
#'
#' @param results output of interaction test
#' @param rep_count number of replicates per concentration in simulation
#' @param concentration_count number of concentrations in simulation

plot_hit_rate_MSstatsResponse = function(results, rep_count,concentration_count ) {

  dt = data.frame(results)

  # Use correct column for filtering
  filtered_dt = dt %>%
    filter(adjust_pval < 0.05)

  # Identify protein groups
  strong_proteins = dt %>%
    filter(grepl("strong", protein)) %>% pull(protein)
  #dt[grepl("strong", protein), protein]
  weak_proteins = dt %>%
    filter(grepl("weak", protein)) %>% pull(protein)
  #weak_proteins   = dt[grepl("weak", protein), protein]
  no_proteins = dt %>%
    filter(grepl("no", protein)) %>% pull(protein)
  #no_proteins     = dt[grepl("no", protein), protein]

  # Function to calculate percentage passing per group
  calculate_percentage = function(group_proteins) {
    total = length(group_proteins)
    if (total == 0) return(NA_real_)
    passed = sum(group_proteins %in% filtered_dt$protein)
    round((passed / total) * 100, 1)
  }

  # Compute percentages
  percent_strong = calculate_percentage(strong_proteins)
  percent_weak   = calculate_percentage(weak_proteins)
  percent_no     = calculate_percentage(no_proteins)

  # Prepare data table for plotting
  plot_data = data.table(
    Category = c("TPR (Strong)", "TPR (Weak)", "FPR"),
    Percent = c(percent_strong, percent_weak, percent_no)
  )

  # Plot
  p = ggplot(plot_data, aes(x = Category, y = Percent, fill = Category)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(Percent, 1), "%")), vjust = -0.3, size = 6, fontface = "bold") +
    labs(title = "Future Experiment Performance",
         subtitle = paste("with", concentration_count, "doses and", rep_count, "replicates per dose"),
         x = "Performance Metrics",
         y = "Rate (%)") +
    scale_fill_manual(values = c("TPR (Strong)" = "#1b9e77", "TPR (Weak)" = "#d95f02", "FPR" = "#7570b3")) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), breaks = seq(0, 100, 20)) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(face = "bold", size = 18),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )

  #print(p)
  return(list(plot = p, plot_data = plot_data))
}


