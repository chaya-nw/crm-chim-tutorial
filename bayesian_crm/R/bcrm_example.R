# Run one illustrative Bayesian CRM trial

run_bcrm_example <- function(
    scenario_no,
    sim_no,
    prior_spec,
    target_prob,
    id70_settings,
    ab_settings,
    cap_method,
    N_cap,
    seed.no,
    start.dose,
    group.size,
    target.doses,
    compiled_model,
    chains,
    iter,
    warmup,
    save.output.to
) {
  
  # -----------------------------------------------------------------------
  # Check selected options
  # -----------------------------------------------------------------------
  
  prior_spec <- match.arg(
    prior_spec,
    choices = c("id70_p", "ab_p")
  )
  
  cap_method <- match.arg(
    cap_method,
    choices = c("capped", "uncapped")
  )
  
  
  # -----------------------------------------------------------------------
  # Read fixed-scenario parameter values
  # -----------------------------------------------------------------------
  
  parameter_file <- paste0(
    "pars_scenarios_",
    prior_spec,
    ".csv"
  )
  
  scenarios <- readr::read_csv(
    here::here(
      "data",
      parameter_file
    ),
    show_col_types = FALSE
  )
  
  if (
    scenario_no < 1 ||
    scenario_no > nrow(scenarios)
  ) {
    stop(
      "`scenario_no` must be between 1 and ",
      nrow(scenarios),
      "."
    )
  }
  
  true_scenario <- scenarios %>%
    dplyr::slice(scenario_no) %>%
    dplyr::select(a, b)
  
  
  # -----------------------------------------------------------------------
  # Create output folder and file paths
  # -----------------------------------------------------------------------
  
  dir.create(
    save.output.to,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  learning_curve_file <- file.path(
    save.output.to,
    "learning_post_curves.csv"
  )
  
  trial_data_file <- file.path(
    save.output.to,
    "learning_trial_dat.csv"
  )
  
  est_probs_file <- file.path(
    save.output.to,
    "est.probs.csv"
  )
  
  est_id70_file <- file.path(
    save.output.to,
    "est.id70s.csv"
  )
  
  
  # Prevent accidental appending to an existing example
  files_to_check <- c(
    learning_curve_file,
    trial_data_file,
    est_probs_file,
    est_id70_file
  )
  
  if (any(file.exists(files_to_check))) {
    stop(
      "Illustrative-trial output files already exist in: ",
      save.output.to,
      ". Remove or rename them before rerunning."
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Initialise the trial
  # -----------------------------------------------------------------------
  
  my.dat <- NULL
  my.summary.stats <- NULL
  
  new.dose <- start.dose
  n.at.target.dose <- 0
  cohort <- 1
  
  design_no <- 1
  row_no <- 1
  
  
  # Use the same scenario-specific seed as the simulation study
  replicate_seed <-
    seed.no +
    100 * scenario_no +
    sim_no
  
  set.seed(replicate_seed)
  
  
  # -----------------------------------------------------------------------
  # Define stopping criterion
  # -----------------------------------------------------------------------
  
  continue_trial <- function() {
    
    if (cap_method == "capped") {
      
      n.at.target.dose < 20 &&
        cohort * group.size <= N_cap
      
    } else {
      
      n.at.target.dose < 20
    }
  }
  
  
  # -----------------------------------------------------------------------
  # Run the sequential illustrative trial
  # -----------------------------------------------------------------------
  
  while (continue_trial()) {
    
    # Simulate the next cohort
    my.dat <- my.sim.fn(
      dat = my.dat,
      scenario = true_scenario,
      new.dose = new.dose,
      new.group.size = group.size
    )
    
    
    # Update the Bayesian CRM and save the posterior learning curve
    summary.stats <- update_crm_model(
      dat = my.dat,
      scenario.no = scenario_no,
      des.no = design_no,
      sim.no = sim_no,
      target.doses = target.doses,
      save.output.to = save.output.to,
      row.numb = row_no,
      compiled_model = compiled_model,
      chains = chains,
      iter = iter,
      warmup = warmup,
      prior_spec = prior_spec,
      target_prob = target_prob,
      id70_settings = id70_settings,
      ab_settings = ab_settings,
      save_learning_curve = TRUE
    )
    
    
    # # Store the interim summary
    # my.summary.stats <- dplyr::bind_rows(
    #   my.summary.stats,
    #   summary.stats
    # )
    
    
    # Apply the standard dose-allocation criteria
    decision <- decision.criteria(
      dat = my.dat,
      summary.stats = summary.stats
    )
    
    n.at.target.dose <- as.numeric(
      decision["target.n"]
    )
    
    proposed_log10_dose <- as.numeric(
      decision["new.dose"]
    )
    
    current_log10_dose <- log10(new.dose)
    
    
    # Apply the maximum two-log10-level escalation rule
    if (
      summary.stats$model.est.id70 >
      proposed_log10_dose &&
      proposed_log10_dose >
      current_log10_dose
    ) {
      
      next_log10_dose <- dplyr::case_when(
        
        summary.stats$model.est.id70 -
          current_log10_dose > 2 ~
          current_log10_dose + 2,
        
        TRUE ~
          summary.stats$model.est.id70
      )
      
    } else {
      
      next_log10_dose <- proposed_log10_dose
    }
    
    
    # Convert back to the original CFU scale
    new.dose <- 10^next_log10_dose
    
    # Store the posterior summary and the dose recommended after this cohort
    summary.stats <- summary.stats %>%
      dplyr::mutate(
        cohort = cohort,
        recommended_next_dose = new.dose
      )
    
    my.summary.stats <- dplyr::bind_rows(
      my.summary.stats,
      summary.stats
    )
    
    # Move to the next cohort
    cohort <- cohort + 1
  }
  
  
  # -----------------------------------------------------------------------
  # Add trial identifiers and true parameter information
  # -----------------------------------------------------------------------
  
  true_a <- true_scenario$a[[1]]
  true_b <- true_scenario$b[[1]]
  
  true_id70 <- true_b * (
    (1 - target_prob)^(-1 / true_a) - 1
  )
  
  trial_data <- my.dat %>%
    dplyr::mutate(
      cohort = dplyr::row_number(),
      scenario = scenario_no,
      design = design_no,
      sim = sim_no,
      true_a = true_a,
      true_b = true_b,
      true_id70 = true_id70
    ) %>%
    dplyr::left_join(
      my.summary.stats %>%
        dplyr::select(
          cohort,
          recommended_next_dose
        ),
      by = "cohort"
    )
  
  
  # -----------------------------------------------------------------------
  # Save trial-level data
  # -----------------------------------------------------------------------
  
  readr::write_csv(
    trial_data,
    trial_data_file
  )
  
  
  # -----------------------------------------------------------------------
  # Read posterior curves produced during the interim updates
  # -----------------------------------------------------------------------
  
  posterior_curves <- readr::read_csv(
    learning_curve_file,
    show_col_types = FALSE
  )
  
  
  # -----------------------------------------------------------------------
  # Return results
  # -----------------------------------------------------------------------
  
  list(
    trial_data        = trial_data,
    posterior_curves  = posterior_curves,
    summary_stats     = my.summary.stats,
    true_scenario     = tibble(
      scenario   = scenario_no,
      a          = true_a,
      b          = true_b,
      ID70       = true_id70,
      ID70_log10 = log10(true_id70)
    ),
    final_administered_dose = last(trial_data$ds),
    final_recommended_dose  = last(trial_data$recommended_next_dose),
    n_cohorts           = nrow(trial_data),
    total_N             = sum(trial_data$nn),
    output_folder       = save.output.to
  )
}