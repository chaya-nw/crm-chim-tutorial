# Function for running rule-based dose-escalation simulations

run_rule_based <- function(
    cap_method,
    generation_method,
    prior_spec,
    n_sim_scenario,
    simulation_seed,
    n_cap
) {
  
  # -----------------------------------------------------------------------
  # Check selected options
  # -----------------------------------------------------------------------
  
  cap_method <- match.arg(
    cap_method,
    choices = c("capped", "uncapped")
  )
  
  generation_method <- match.arg(
    generation_method,
    choices = c("random", "scenarios")
  )
  
  prior_spec <- match.arg(
    prior_spec,
    choices = c("id70_p", "ab_p")
  )
  
  
  # -----------------------------------------------------------------------
  # Read true parameter values
  # -----------------------------------------------------------------------
  
  parameter_file <- paste0(
    "pars_",
    generation_method,
    "_",
    prior_spec,
    ".csv"
  )
  
  parameter_path <- here::here(
    "data",
    parameter_file
  )
  
  if (!file.exists(parameter_path)) {
    stop(
      "The parameter file does not exist: ",
      parameter_path
    )
  }
  
  parameter_values <- readr::read_csv(
    parameter_path,
    show_col_types = FALSE
  )
  
  
  # Create a list in which to store all simulated trials
  simulation_results <- list()
  
  
  # -----------------------------------------------------------------------
  # Random parameter sets
  # -----------------------------------------------------------------------
  
  if (generation_method == "random") {
    
    n_trials <- nrow(parameter_values)
    
    for (j in seq_len(n_trials)) {
      
      set.seed(simulation_seed + j)
      
      simulation_results[[j]] <- rulebased_path(
        b = parameter_values$b[j],
        a = parameter_values$a[j],
        scenario = j,
        sim = j,
        cap_method = cap_method,
        N_cap = n_cap
      )
    }
    
    
    # -----------------------------------------------------------------------
    # Fixed scenarios
    # -----------------------------------------------------------------------
    
  } else {
    
    parameter_values <- parameter_values %>%
      dplyr::mutate(
        scenario = dplyr::row_number()
      )
    
    result_index <- 1
    
    for (s in parameter_values$scenario) {
      
      a_s <- parameter_values$a[s]
      b_s <- parameter_values$b[s]
      
      for (i in seq_len(n_sim_scenario)) {
        
        # Gives each scenario-replicate combination a unique seed
        trial_seed <-
          simulation_seed +
          (s - 1) * n_sim_scenario +
          i
        
        set.seed(trial_seed)
        
        simulation_results[[result_index]] <- rulebased_path(
          b = b_s,
          a = a_s,
          scenario = s,
          sim = i,
          cap_method = cap_method,
          N_cap = n_cap
        )
        
        result_index <- result_index + 1
      }
    }
  }
  
  
  # -----------------------------------------------------------------------
  # Combine all simulated trials
  # -----------------------------------------------------------------------
  
  simulation_results <- dplyr::bind_rows(
    simulation_results
  )
  
  return(simulation_results)
}