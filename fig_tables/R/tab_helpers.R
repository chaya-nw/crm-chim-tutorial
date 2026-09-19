# Create the sample-size summary table from the four simulation summaries - Tab2

create_sample_size_tab <- function(simulation_summaries) {
  
  # -----------------------------------------------------------------------
  # Extract sample-size summaries
  # -----------------------------------------------------------------------
  
  rb_uncapped <-
    simulation_summaries$uncapped$rule_based$alloc_overall
  
  bcrm_uncapped <-
    simulation_summaries$uncapped$bcrm$overall_N
  
  rb_capped <-
    simulation_summaries$capped$rule_based$alloc_overall
  
  bcrm_capped <-
    simulation_summaries$capped$bcrm$overall_N
  
  
  # -----------------------------------------------------------------------
  # Format one column of the table
  # -----------------------------------------------------------------------
  
  format_sample_size_results <- function(
    results,
    q95_variable
  ) {
    
    q95_value <- results[[q95_variable]]
    
    c(
      # Average total N [95th percentile]
      paste0(
        round(results$avg_total_N, 1),
        " [",
        round(q95_value, 1),
        "]"
      ),
      
      # Mean N at final selected dose
      as.character(
        round(results$mean_at_final, 1)
      ),
      
      # Mean N above final selected dose
      as.character(
        round(results$mean_above_final, 1)
      ),
      
      # Median total N
      as.character(
        round(results$median_total_N, 1)
      ),
      
      # Median N at final selected dose
      as.character(
        round(results$median_at_final, 1)
      ),
      
      # Median N above final selected dose
      as.character(
        round(results$median_above_final, 1)
      ),
      
      # Median number of cohorts
      as.character(
        round(results$median_cohorts, 1)
      )
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Construct the table dataset
  # -----------------------------------------------------------------------
  
  tibble::tibble(
    Statistic = c(
      "Average total $N$ [95th percentile]",
      "Mean $N$ at target",
      "Mean $N$ above target",
      "Median total $N$",
      "Median $N$ at target",
      "Median $N$ above target",
      "Median number of cohorts"
    ),
    
    rb_uncapped = format_sample_size_results(
      results = rb_uncapped,
      q95_variable = "q95_total_N"
    ),
    
    bcrm_uncapped = format_sample_size_results(
      results = bcrm_uncapped,
      q95_variable = "q_95_N"
    ),
    
    rb_capped = format_sample_size_results(
      results = rb_capped,
      q95_variable = "q95_total_N"
    ),
    
    bcrm_capped = format_sample_size_results(
      results = bcrm_capped,
      q95_variable = "q_95_N"
    )
  )
}








# Create the scenario-specific operating-characteristics table - Table 3

create_scenario_oc_tab <- function(simulation_summaries) {
  
  # -----------------------------------------------------------------------
  # Extract capped scenario-specific summaries
  # -----------------------------------------------------------------------
  
  rb_allocation <-
    simulation_summaries$capped$rule_based$alloc_by_scenario %>%
    dplyr::arrange(scenario)
  
  bcrm_allocation <-
    simulation_summaries$capped$bcrm$summary_by_scen %>%
    dplyr::arrange(scenario)
  
  bcrm_width <-
    simulation_summaries$capped$bcrm$width_by_scenario %>%
    dplyr::arrange(scenario)
  
  scenario_truth <-
    simulation_summaries$capped$bcrm$pcs_by_scenario %>%
    dplyr::arrange(scenario)
  
  
  # -----------------------------------------------------------------------
  # Formatting helpers
  # -----------------------------------------------------------------------
  
  format_number <- function(x) {
    format(
      round(x, 1),
      trim = TRUE,
      scientific = FALSE
    )
  }
  
  format_mean_q95 <- function(mean_value, q95_value) {
    paste0(
      format_number(mean_value),
      " [",
      format_number(q95_value),
      "]"
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Construct the table dataset
  # -----------------------------------------------------------------------
  
  tibble::tibble(
    `Operating characteristic` = c(
      "True $ID_{70}$ ($\\log_{10}$ CFU)",
      
      # Panel A: Rule-based design
      "Average total $N$ [95th percentile]",
      "Average width of 90% CrI ($\\log_{10}$ CFU)",
      "Average width of 50% CrI ($\\log_{10}$ CFU)",
      "Mean $N$ at target",
      "Median number of cohorts",
      
      # Panel B: Bayesian CRM
      "Average total $N$ [95th percentile]",
      "Average width of 90% CrI ($\\log_{10}$ CFU)",
      "Average width of 50% CrI ($\\log_{10}$ CFU)",
      "Mean $N$ at target",
      "Median number of cohorts"
    ),
    
    `1` = c(
      format_number(scenario_truth$true_log10_ID70[1]),
      format_mean_q95(
        rb_allocation$avg_total_N[1],
        rb_allocation$q95_total_N[1]
      ),
      "-",
      "-",
      format_number(rb_allocation$mean_at_final[1]),
      format_number(rb_allocation$median_cohorts[1]),
      format_mean_q95(
        bcrm_allocation$avg_total_N[1],
        bcrm_allocation$q_95_N[1]
      ),
      format_number(bcrm_width$mean_width_90[1]),
      format_number(bcrm_width$mean_width_50[1]),
      format_number(bcrm_allocation$mean_at_final[1]),
      format_number(bcrm_allocation$median_cohorts[1])
    ),
    
    `2` = c(
      format_number(scenario_truth$true_log10_ID70[2]),
      format_mean_q95(
        rb_allocation$avg_total_N[2],
        rb_allocation$q95_total_N[2]
      ),
      "-",
      "-",
      format_number(rb_allocation$mean_at_final[2]),
      format_number(rb_allocation$median_cohorts[2]),
      format_mean_q95(
        bcrm_allocation$avg_total_N[2],
        bcrm_allocation$q_95_N[2]
      ),
      format_number(bcrm_width$mean_width_90[2]),
      format_number(bcrm_width$mean_width_50[2]),
      format_number(bcrm_allocation$mean_at_final[2]),
      format_number(bcrm_allocation$median_cohorts[2])
    ),
    
    `3` = c(
      format_number(scenario_truth$true_log10_ID70[3]),
      format_mean_q95(
        rb_allocation$avg_total_N[3],
        rb_allocation$q95_total_N[3]
      ),
      "-",
      "-",
      format_number(rb_allocation$mean_at_final[3]),
      format_number(rb_allocation$median_cohorts[3]),
      format_mean_q95(
        bcrm_allocation$avg_total_N[3],
        bcrm_allocation$q_95_N[3]
      ),
      format_number(bcrm_width$mean_width_90[3]),
      format_number(bcrm_width$mean_width_50[3]),
      format_number(bcrm_allocation$mean_at_final[3]),
      format_number(bcrm_allocation$median_cohorts[3])
    ),
    
    `4` = c(
      format_number(scenario_truth$true_log10_ID70[4]),
      format_mean_q95(
        rb_allocation$avg_total_N[4],
        rb_allocation$q95_total_N[4]
      ),
      "-",
      "-",
      format_number(rb_allocation$mean_at_final[4]),
      format_number(rb_allocation$median_cohorts[4]),
      format_mean_q95(
        bcrm_allocation$avg_total_N[4],
        bcrm_allocation$q_95_N[4]
      ),
      format_number(bcrm_width$mean_width_90[4]),
      format_number(bcrm_width$mean_width_50[4]),
      format_number(bcrm_allocation$mean_at_final[4]),
      format_number(bcrm_allocation$median_cohorts[4])
    ),
    
    `5` = c(
      format_number(scenario_truth$true_log10_ID70[5]),
      format_mean_q95(
        rb_allocation$avg_total_N[5],
        rb_allocation$q95_total_N[5]
      ),
      "-",
      "-",
      format_number(rb_allocation$mean_at_final[5]),
      format_number(rb_allocation$median_cohorts[5]),
      format_mean_q95(
        bcrm_allocation$avg_total_N[5],
        bcrm_allocation$q_95_N[5]
      ),
      format_number(bcrm_width$mean_width_90[5]),
      format_number(bcrm_width$mean_width_50[5]),
      format_number(bcrm_allocation$mean_at_final[5]),
      format_number(bcrm_allocation$median_cohorts[5])
    )
  )
}











# Create the scenario-specific performance table (Bias, RMSE, Cov, )

create_scenario_performance_tab <- function(
    simulation_summaries
) {
  
  # -----------------------------------------------------------------------
  # Extract capped scenario-specific summaries
  # -----------------------------------------------------------------------
  
  rb_allocation <-
    simulation_summaries$capped$rule_based$alloc_by_scenario %>%
    dplyr::arrange(scenario)
  
  rb_extremes <-
    simulation_summaries$capped$rule_based$extremes_by_scenario %>%
    dplyr::arrange(scenario)
  
  bcrm_allocation <-
    simulation_summaries$capped$bcrm$summary_by_scen %>%
    dplyr::arrange(scenario)
  
  bcrm_bias <-
    simulation_summaries$capped$bcrm$bias_summary %>%
    dplyr::arrange(scenario)
  
  bcrm_rmse <-
    simulation_summaries$capped$bcrm$rmse_by_scenario %>%
    dplyr::arrange(scenario)
  
  bcrm_coverage <-
    simulation_summaries$capped$bcrm$coverage_by_scenario %>%
    dplyr::arrange(scenario)
  
  bcrm_extremes <-
    simulation_summaries$capped$bcrm$pcs_extremes %>%
    dplyr::arrange(scenario)
  
  scenario_truth <-
    simulation_summaries$capped$bcrm$pcs_by_scenario %>%
    dplyr::arrange(scenario)
  
  
  # -----------------------------------------------------------------------
  # Format numeric results
  # -----------------------------------------------------------------------
  
  format_number <- function(x) {
    format(
      round(x, 2),
      trim = TRUE,
      scientific = FALSE
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Construct the table
  # -----------------------------------------------------------------------
  
  table_values <- rbind(
    # True scenario value
    format_number(
      scenario_truth$true_log10_ID70
    ),
    
    # Panel A: Rule-based design
    format_number(
      100 * rb_extremes$p_end_low
    ),
    format_number(
      100 * rb_extremes$p_end_high
    ),
    format_number(
      rb_allocation$mean_above_final
    ),
    
    # Panel B: Bayesian CRM
    format_number(
      bcrm_bias$mean_abs_bias
    ),
    format_number(
      bcrm_rmse$rmse_log10
    ),
    format_number(
      bcrm_coverage$coverage_90
    ),
    format_number(
      bcrm_extremes$pct_lowest
    ),
    format_number(
      bcrm_extremes$pct_highest
    ),
    format_number(
      bcrm_allocation$mean_above_final
    )
  )
  
  colnames(table_values) <- as.character(
    scenario_truth$scenario
  )
  
  dplyr::bind_cols(
    tibble::tibble(
      `Operating characteristic` = c(
        "True $ID_{70}$ ($\\log_{10}$ CFU)",
        
        # Panel A: Rule-based design
        "End at lowest dose (%)",
        "End at highest dose (%)",
        "Mean $N$ above target",
        
        # Panel B: Bayesian CRM
        "Mean absolute bias",
        "RMSE of $ID_{70}$",
        "Coverage of true $ID_{70}$",
        "End at lowest dose (%)",
        "End at highest dose (%)",
        "Mean $N$ above target"
      )
    ),
    
    tibble::as_tibble(
      table_values,
      .name_repair = "minimal"
    )
  )
}






# Calculate the prior-implied skeleton
summarise_prior_implied_skeleton <- function(
    nsamp,
    id70_settings,
    ab_settings,
    candidate_doses = 10^(3:7),
    seed = 2025
) {
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # Primary prior: prior on a and log10(ID70)
  # ------------------------------------------------------------
  
  a_primary <- exp(
    runif(
      nsamp,
      min = log(id70_settings$a_lower),
      max = log(id70_settings$a_upper)
    )
  )
  
  log10_id70 <- truncnorm::rtruncnorm(
    nsamp,
    a = id70_settings$log10_id70_lower,
    b = id70_settings$log10_id70_upper,
    mean = id70_settings$log10_id70_mean,
    sd = id70_settings$log10_id70_sd
  )
  
  id70_primary <- 10^log10_id70
  
  b_primary <- id70_primary /
    (0.30^(-1 / a_primary) - 1)
  
  primary_draws <- tibble::tibble(
    draw = seq_len(nsamp),
    a = a_primary,
    b = b_primary
  )
  
  # ------------------------------------------------------------
  # Alternative prior: prior on a and b
  # ------------------------------------------------------------
  
  a_alternative <- exp(
    runif(
      nsamp,
      min = log(ab_settings$a_lower),
      max = log(ab_settings$a_upper)
    )
  )
  
  b_alternative <- exp(
    runif(
      nsamp,
      min = log(ab_settings$b_lower),
      max = log(ab_settings$b_upper)
    )
  )
  
  alternative_draws <- tibble::tibble(
    draw = seq_len(nsamp),
    a = a_alternative,
    b = b_alternative
  )
  
  # ------------------------------------------------------------
  # Function to calculate summaries at candidate doses
  # ------------------------------------------------------------
  
  calculate_summaries <- function(parameter_draws, prior_name) {
    
    tidyr::crossing(
      parameter_draws,
      dose_level = seq_along(candidate_doses)
    ) %>%
      dplyr::mutate(
        challenge_dose = candidate_doses[dose_level],
        
        infection_probability =
          1 - (1 + challenge_dose / b)^(-a)
      ) %>%
      dplyr::group_by(dose_level, challenge_dose) %>%
      dplyr::summarise(
        prior_median = stats::median(infection_probability),
        lower_95 = stats::quantile(
          infection_probability,
          probs = 0.025
        ),
        upper_95 = stats::quantile(
          infection_probability,
          probs = 0.975
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(prior_specification = prior_name)
  }
  
  primary_summary <- calculate_summaries(
    primary_draws,
    "Primary"
  )
  
  alternative_summary <- calculate_summaries(
    alternative_draws,
    "Alternative"
  )
  
  # ------------------------------------------------------------
  # Combine results into the structure needed for the table
  # ------------------------------------------------------------
  
  dplyr::bind_rows(
    primary_summary,
    alternative_summary
  ) %>%
    tidyr::pivot_wider(
      id_cols = c(dose_level, challenge_dose),
      names_from = prior_specification,
      values_from = c(prior_median, lower_95, upper_95),
      names_glue = "{prior_specification}_{.value}"
    ) %>%
    dplyr::arrange(dose_level)
}