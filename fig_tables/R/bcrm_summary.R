# Summarise Bayesian CRM simulation results


# Find the candidate dose nearest to a specified dose on the log10 scale
nearest_grid <- function(x, grid) {
  grid[which.min(abs(log10(grid) - log10(x)))]
}


# Summarise operating characteristics of the Bayesian CRM design
#
# Arguments:
# sims_file:          path to cohort-level CRM simulation results
# sumstats_file:      path to posterior infection-probability summaries
# est_id_file:        path to posterior log10(ID70) summaries
# pars_file:          path to the true parameter or scenario file
# generation_method:  either "random" or "scenarios"
# dose_grid:          candidate dose grid
# tol_log:            maximum log10 distance used to define allocation near the true ID70

sum_bcrm <- function(
    sims_file,
    sumstats_file,
    est_id_file,
    pars_file,
    generation_method,
    dose_grid = 10^(3:7),
    tol_log = 0.5
) {
  
  
  
  generation_method <- match.arg(
    generation_method,
    choices = c("random", "scenarios")
  )
  
  
  # =======================================================================
  # Read and prepare data
  # =======================================================================
  pars <- readr::read_csv(pars_file, show_col_types = FALSE)
  
  # Standardize the scenario variable
  if (!"scenario" %in% names(pars)) {
    
    if (
      generation_method == "random" &&
      "sim" %in% names(pars)
    ) {
      
      pars <- pars %>%
        dplyr::rename(
          scenario = sim
        )
      
    } else {
      
      pars <- pars %>%
        dplyr::mutate(
          scenario = dplyr::row_number()
        )
    }
  }
  
  # Standardize the log10(ID70) variable
  if (
    !"ID70_log10" %in% names(pars) &&
    "log10_id70" %in% names(pars)
  ) {
    
    pars <- pars %>%
      dplyr::rename(
        ID70_log10 = log10_id70
      )
  }
  
  
  # Cohort-level CRM simulation results
  sims <- readr::read_csv(sims_file, col_names = FALSE, show_col_types = FALSE)
  
  colnames(sims) <- c(
    "ds",
    "obs",
    "nn",
    "scenario",
    "design",
    "sim",
    "kk"
  )
  
  
  # Posterior infection-probability summaries
  sumstats <- readr::read_csv(sumstats_file, col_names = FALSE, show_col_types = FALSE)
  
  colnames(sumstats) <- c(
    "ds",
    "q5",
    "q25",
    "q50",
    "q75",
    "q95",
    "model_est_id70_idx",
    "scenario",
    "design",
    "sim",
    "kk"
  )
  
  
  # Posterior summaries of log10(ID70)
  est_id <- readr::read_csv(est_id_file, col_names = FALSE, show_col_types = FALSE)
  
  colnames(est_id) <- c(
    "q5",
    "q25",
    "q50",
    "q75",
    "q95",
    "sim",
    "n_doses",
    "design",
    "scenario",
    "rn"
  )
  
  
  # ===============================================================================
  # 1. Correct Selection (PCS): PCS = how often the final selected dose equals the 
  # nearest grid dose to the true ID70
  # ===============================================================================
  
  scen <- pars %>%
    mutate(
      nearest_dose   = purrr::map_dbl(ID70, nearest_grid, grid = dose_grid),
      nearest_log10  = log10(nearest_dose),
      out_of_grid    = ID70_log10 < log10(min(dose_grid)) | ID70_log10 > log10(max(dose_grid)),
      
      #Lower and upper grid points around true ID70
      lower_log10   = pmax(floor(ID70_log10), log10(min(dose_grid))),
      upper_log10   = pmin(ceiling(ID70_log10), log10(max(dose_grid)))
    )
  
  
  # -------------------------------------------------------------------------
  # Final selected dose and probability of correct selection
  # -------------------------------------------------------------------------
  final_dose <- sims %>%
    group_by(scenario, design, sim) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    transmute(
      scenario,
      design,
      sim,
      final_dose = ds,
      final_log10 = log10(ds)
    )
  
  
  res <- final_dose %>%
    left_join(
      scen %>%
        select(
          scenario,
          ID70,
          ID70_log10,
          nearest_dose,
          nearest_log10,
          out_of_grid,
          lower_log10,
          upper_log10
        ),
      by = "scenario"
    ) %>%
    mutate(
      # 1) PCS: exactly nearest grid dose
      correct_exact   = as.integer(final_dose == nearest_dose),
      
      # 2) PCS: either lower or upper grid around true ID70
      correct_two     = as.integer(abs(final_log10 - lower_log10) <= 1e-12 |
                                 abs(final_log10 - upper_log10) <= 1e-12),
      
      # 3) PCS: within +-1 log10 of nearest grid dose
      correct_within1 = as.integer(abs(final_log10 - nearest_log10) <= 1 + 1e-12),
      
      # under_true_grid = as.integer(final_log10 < nearest_log10),
      # equal_true_grid = as.integer(abs(final_log10 - nearest_log10) <= 1e-12),
      # over_true_grid = as.integer(final_log10 > nearest_log10)
    )
  
  
  pcs_by_scenario <- res %>%
    group_by(scenario) %>%
    summarise(
      n_trials           = n(),
      true_log10_ID70    = first(ID70_log10),
      nearest_log10      = first(nearest_log10),
      out_of_grid        = first(out_of_grid),
      pcs_exact          = mean(correct_exact),      # final dose == nearest grid only
      pcs_two            = mean(correct_two),        # final dose is one of the two nearest grid points
      pcs_within_1log    = mean(correct_within1),    # final dose within +-1 log of nearest grid
      .groups            = "drop"
    )
  
  overall_pcs <- res %>%
    summarise(
      n_trials          = dplyr::n(),
      pcs_exact         = mean(correct_exact),
      pcs_two           = mean(correct_two),
      pcs_within_1log   = mean(correct_within1)
    )
  
  
  
  # =======================================================================================
  # 2. How often the model-estimated ID70 (from CRM posterior mode bin) hits the true ID70
  # =======================================================================================
  
  # keep final cohort per trial
  est_last_category <- sumstats %>%
    group_by(scenario, design, sim) %>%
    mutate(cohort = row_number()) %>%
    slice_max(order_by = cohort, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      est_log10_mode = pmin(pmax(model_est_id70_idx,min(log10(dose_grid))), max(log10(dose_grid))),
      est_mode_dose = 10^est_log10_mode) %>%
      select(
      scenario,
      design,
      sim,
      est_log10_mode,
      est_mode_dose
    )
  
  
  pcs_est_trial <- est_last_category %>%
    dplyr::left_join(
      scen,
      by = "scenario"
    ) %>%
    dplyr::mutate(
      correct_est_exact = as.integer(
        est_mode_dose == nearest_dose
      ),
      
      correct_est_within1 = as.integer(
        abs(est_log10_mode - nearest_log10) <= 1 + 1e-12
      )
    )
  
  
  pcs_by_estimate_each <- pcs_est_trial %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      true_log10_ID70 = dplyr::first(ID70_log10),
      nearest_log10 = dplyr::first(nearest_log10),
      pcs_est_exact = mean(correct_est_exact),
      pcs_est_within_1log = mean(correct_est_within1),
      .groups = "drop"
    )
  
  
  overall_pcs_est <- pcs_est_trial %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      pcs_est_exact = mean(correct_est_exact),
      pcs_est_within_1log = mean(correct_est_within1)
    )
  
  
  # =======================================================================================
  # 3. Bias, RMSE, interval width and coverage for log10(ID70)
  # =======================================================================================
  est_last <- est_id %>%
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    dplyr::slice_max(
      order_by = n_doses,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      scen %>%
        dplyr::select(
          scenario,
          ID70,
          ID70_log10
        ),
      by = "scenario"
    ) %>%
    dplyr::mutate(
      est_log10 = q50,
      bias_log10 = est_log10 - ID70_log10,
      bias_abs = abs(bias_log10),
      err_log10 = est_log10 - ID70_log10,
      width_90 = q95 - q5,
      width_50 = q75 - q25,
      covered_90 =
        ID70_log10 >= q5 &
        ID70_log10 <= q95
    )
  
  
  bias_summary <- est_last %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      true_log10 = dplyr::first(ID70_log10),
      mean_bias = mean(bias_log10),
      median_bias = median(bias_log10),
      mean_abs_bias = mean(bias_abs),
      .groups = "drop"
    )
  
  
  overall_bias <- est_last %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      mean_bias = mean(bias_log10),
      median_bias = median(bias_log10),
      mean_abs_bias = mean(bias_abs)
    )
  
  
  rmse_by_scenario <- est_last %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      true_log10 = dplyr::first(ID70_log10),
      rmse_log10 = sqrt(mean(err_log10^2)),
      rmse_fold = 10^rmse_log10,
      .groups = "drop"
    )
  
  
  overall_rmse <- est_last %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      rmse_log10 = sqrt(mean(err_log10^2)),
      rmse_fold = 10^rmse_log10
    )
  
  
  width_by_scenario <- est_last %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      mean_width_90 = mean(width_90),
      mean_width_50 = mean(width_50),
      .groups = "drop"
    )
  
  
  overall_width <- est_last %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      mean_width_90 = mean(width_90),
      mean_width_50 = mean(width_50)
    )
  
  
  coverage_by_scenario <- est_last %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      coverage_90 = mean(covered_90),
      .groups = "drop"
    )
  
  
  overall_coverage <- est_last %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      coverage_90 = mean(covered_90)
    )
  
  
  # =======================================================================================
  # 4. Sample size and participant allocation
  # =======================================================================================
  
  n_target <- sims %>%
    dplyr::group_by(
      scenario,
      design,
      sim,
      ds
    ) %>%
    dplyr::summarise(
      n_at_dose = sum(nn),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      final_dose,
      by = c(
        "scenario",
        "design",
        "sim"
      )
    ) %>%
    dplyr::filter(
      ds == final_dose
    ) %>%
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    dplyr::summarise(
      n_at_target = sum(n_at_dose),
      .groups = "drop"
    )
  
  
  n_near_true <- sims %>%
    dplyr::left_join(
      scen %>%
        dplyr::select(
          scenario,
          ID70,
          ID70_log10
        ),
      by = "scenario"
    ) %>%
    dplyr::left_join(
      final_dose,
      by = c(
        "scenario",
        "design",
        "sim"
      )
    ) %>%
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    dplyr::summarise(
      total_N = sum(nn),
      
      n_near_true = sum(
        nn[
          abs(log10(ds) - ID70_log10) <= tol_log
        ]
      ),
      
      n_above_true = sum(
        nn[
          ds > ID70
        ]
      ),
      
      n_above_final = sum(
        nn[
          ds > final_dose
        ]
      ),
      
      n_cohorts = dplyr::n(),
      .groups = "drop"
    )
  
  
  combined <- n_near_true %>%
    dplyr::left_join(
      n_target,
      by = c(
        "scenario",
        "design",
        "sim"
      )
    )
  
  
  summary_by_scen <- combined %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      avg_total_N = mean(total_N),
      q_95_N = quantile(total_N, 0.95),
      mean_at_final = mean(n_at_target),
      mean_at_true = mean(n_near_true),
      mean_above_true = mean(n_above_true),
      mean_above_final = mean(n_above_final),
      median_cohorts = median(n_cohorts),
      mean_cohorts = mean(n_cohorts),
      median_total_N = median(total_N),
      median_at_final = median(n_at_target),
      median_above_final = median(n_above_final),
      .groups = "drop"
    )
  
  
  overall_N <- combined %>%
    dplyr::summarise(
      avg_total_N = mean(total_N),
      q_95_N = quantile(total_N, 0.95),
      mean_at_final = mean(n_at_target),
      mean_at_true = mean(n_near_true),
      mean_above_true = mean(n_above_true),
      mean_above_final = mean(n_above_final),
      median_cohorts = median(n_cohorts),
      mean_cohorts = mean(n_cohorts),
      median_total_N = median(total_N),
      median_at_final = median(n_at_target),
      median_above_final = median(n_above_final)
    )
  
  
  # =======================================================================================
  # 5. Trials ending at the lowest or highest candidate dose
  # =======================================================================================
  
  pcs_extremes <- final_dose %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      pct_lowest = mean(
        final_dose == min(dose_grid)
      ) * 100,
      pct_highest = mean(
        final_dose == max(dose_grid)
      ) * 100,
      .groups = "drop"
    )
  
  
  overall_extremes <- final_dose %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      pct_lowest = mean(
        final_dose == min(dose_grid)
      ) * 100,
      pct_highest = mean(
        final_dose == max(dose_grid)
      ) * 100
    )
  
  
  # =======================================================================
  # Return summaries for the selected generation method
  # =======================================================================
  
  if (generation_method == "random") {
    
    list(
      overall_pcs = overall_pcs,
      overall_pcs_est = overall_pcs_est,
      overall_bias = overall_bias,
      overall_rmse = overall_rmse,
      overall_width = overall_width,
      overall_coverage = overall_coverage,
      overall_N = overall_N,
      overall_extremes = overall_extremes
    )
    
  } else {
    
    list(
      pcs_by_scenario = pcs_by_scenario,
      pcs_by_estimate_each = pcs_by_estimate_each,
      bias_summary = bias_summary,
      rmse_by_scenario = rmse_by_scenario,
      width_by_scenario = width_by_scenario,
      coverage_by_scenario = coverage_by_scenario,
      summary_by_scen = summary_by_scen,
      pcs_extremes = pcs_extremes
    )
  }
}