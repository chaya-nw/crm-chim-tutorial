# Functions for summarising rule-based dose-finding simulations


# Find the candidate dose nearest to a specified dose on the log10 scale
nearest_grid <- function(x, grid) {
  grid[which.min(abs(log10(grid) - log10(x)))]
}


# Summarise operating characteristics of the rule-based design
#
# Arguments:
# sims_file:          path to the simulated rule-based trial data
# pars_file:          path to the true parameter values
# generation_method: "random" or "scenarios"
# dose_grid:          candidate dose grid
# tol_log:            maximum log10 distance used to define allocation near the true ID70


sum_rb <- function(
    sims_file,
    pars_file,
    generation_method,
    dose_grid,
    tol_log = 0.5
) {
  
  generation_method <- match.arg(
    generation_method,
    choices = c("random", "scenarios")
  )
  
  # -----------------------------------------------------------------------
  # Read simulation results and true parameter values
  # -----------------------------------------------------------------------
  sims <- readr::read_csv(sims_file, show_col_types = FALSE)
  
  scen <- readr::read_csv(pars_file, show_col_types = FALSE) %>%
    mutate(scenario = row_number())
  
  # Standardise the log10(ID70) column name used by the random-parameter
  # and scenario files
  if (
    !"ID70_log10" %in% names(scen) &&
    "log10_id70" %in% names(scen)
  ) {
    scen <- scen %>%
      dplyr::rename(
        ID70_log10 = log10_id70
      )
  }
  
  # =======================================================================
  # Construct true ID70 information
  # =======================================================================
  scen_truth <- scen %>%
    mutate(
      true_ID70 = ID70,
      true_log10 = ID70_log10,
      nearest_dose = purrr::map_dbl(true_ID70, nearest_grid, grid = dose_grid),
      nearest_log10 = log10(nearest_dose),
      
      # Grid points immediately below and above the true ID70
      lower_log10 = pmax(floor(true_log10), log10(min(dose_grid))),
      upper_log10 = pmin(ceiling(true_log10), log10(max(dose_grid)))
      ) %>%
        select(
          scenario,
          true_ID70,
          true_log10,
          nearest_dose,
          nearest_log10,
          lower_log10,
          upper_log10
        )
  
  # =======================================================================
  # 1. Final selected dose for each simulated trial
  # =======================================================================
  final_dose <- sims %>%
    group_by(scenario, sim) %>%
    summarise(
      final_dose = dplyr::last(dose),
      final_log10 = log10(final_dose),
      .groups = "drop"
    )
  
  # =======================================================================
  # 2. Correct selection and under/over-selection (grid-based)
  # =======================================================================
  res <- final_dose %>%
    left_join(scen_truth,by = "scenario") %>%
    mutate(
      under_true_grid = as.integer(final_log10 < nearest_log10),
      equal_true_grid = as.integer(abs(final_log10 - nearest_log10) <= 1e-12),
      over_true_grid  = as.integer(final_log10 > nearest_log10),
      
      # Final dose is one of the two grid doses bracketing the true ID70
      correct_two = as.integer(
        abs(final_log10 - lower_log10) <= 1e-12 |
          abs(final_log10 - upper_log10) <= 1e-12
      )
    )
  
  pcs_by_scenario <- res %>%
    group_by(scenario) %>%
    summarise(
      n_trials         = n(),
      true_log10_ID70  = first(true_log10),
      pcs_exact        = mean(equal_true_grid),  # final dose == nearest grid
      pcs_two          = mean(correct_two),      # final dose is lower OR upper grid
      p_select_below   = mean(under_true_grid),
      p_select_above   = mean(over_true_grid),
      .groups = "drop"
    )
  
  pcs_overall <- res %>%
    summarise(
      n_trials         = n(),
      pcs_exact        = mean(equal_true_grid),
      pcs_two          = mean(correct_two),
      p_select_below   = mean(under_true_grid),
      p_select_above   = mean(over_true_grid)
    )
  
  # =======================================================================
  # 3. Allocation summaries
  # =======================================================================
  # Number of participants treated at the final selected dose
  n_at_final <- sims %>%
    group_by(scenario, sim, dose) %>%
    summarise(n_at_dose = sum(total), .groups = "drop") %>%
    left_join(final_dose, by = c("scenario", "sim")) %>%
    filter(dose == final_dose) %>%
    group_by(scenario,sim) %>%
    summarise(
      n_at_final = sum(n_at_dose),
      .groups = "drop"
    )
  
  # Total sample size and allocation relative to the true (true ID70) and final doses
  n_near_true <- sims %>%
    left_join(scen_truth, by = "scenario") %>%
    left_join(final_dose, by = c("scenario", "sim")) %>%
    group_by(scenario, sim) %>%
    summarise(
      total_N       = sum(total),
      n_near_true   = sum(total[abs(log10(dose) - true_log10) <= tol_log]),
      n_above_true  = sum(total[dose > true_ID70]),
      n_above_final = sum(total[dose > final_dose]),      # above final dose
      n_cohorts     = n(),
      .groups       = "drop"
    )
  
  # Combine allocation and final-dose information
  alloc_combined <- n_near_true %>%
    left_join(n_at_final, by = c("scenario", "sim")) %>%
    left_join(
      res %>%
        select(
          scenario,
          sim,
          final_dose,
          final_log10
        ),
      by = c("scenario", "sim")
    )
  
  alloc_by_scenario <- alloc_combined %>%
    group_by(scenario) %>%
    summarise(
      avg_total_N        = mean(total_N),
      q95_total_N        = quantile(total_N, 0.95),
      mean_at_final      = mean(n_at_final),
      mean_at_true       = mean(n_near_true),
      mean_above_true    = mean(n_above_true),
      mean_above_final   = mean(n_above_final),
      median_cohorts     = median(n_cohorts),
      median_total_N     = median(total_N),
      median_at_final    = median(n_at_final),
      median_above_final = median(n_above_final),
      .groups = "drop"
    )
  
  alloc_overall <- alloc_combined %>%
    summarise(
      avg_total_N        = mean(total_N),
      q95_total_N        = quantile(total_N, 0.95),
      mean_at_final      = mean(n_at_final),
      mean_at_true       = mean(n_near_true),
      mean_above_true    = mean(n_above_true),
      mean_above_final   = mean(n_above_final),
      median_cohorts     = median(n_cohorts),
      median_total_N     = median(total_N),
      median_at_final    = median(n_at_final),
      median_above_final = median(n_above_final)
    )
  
  
  
  # =======================================================================
  # 4. Extremes: Trials ending at the lowest or highest candidate dose
  # =======================================================================
  extremes_by_scenario <- final_dose %>%
    group_by(scenario) %>%
    summarise(
      n_trials   = n(),
      p_end_low  = mean(final_dose == min(dose_grid)),
      p_end_high = mean(final_dose == max(dose_grid)),
      .groups    = "drop"
    )
  
  extremes_overall <- final_dose %>%
    summarise(
      p_end_low  = mean(final_dose == min(dose_grid)),
      p_end_high = mean(final_dose == max(dose_grid))
    )
  
  
  # =======================================================================
  # Return the summaries required for the selected generation method
  # =======================================================================
  
  if (generation_method == "random") {
    
    return(
      list(
        pcs_overall          = pcs_overall,
        alloc_overall        = alloc_overall,
        extremes_overall     = extremes_overall
      )
    )
    
  } else {
    
    return(
      list(
        pcs_by_scenario      = pcs_by_scenario,
        alloc_by_scenario    = alloc_by_scenario,
        extremes_by_scenario = extremes_by_scenario
      )
    )
  }
}