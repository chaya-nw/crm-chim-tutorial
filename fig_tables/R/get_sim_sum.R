# Obtain operating-characteristic summaries for one selected design

get_simulation_summary <- function(
    design_method,
    cap_method,
    generation_method,
    prior_spec,
    dose_grid
) {



# -------------------------------------------------------------------------
# Check selected options
# -------------------------------------------------------------------------
design_method <- match.arg(
  design_method,
  choices = c("rule_based", "bcrm")
)

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


# -------------------------------------------------------------------------
# Construct parameter-file path
# -------------------------------------------------------------------------
if (generation_method == "random") {
  
  pars_file <- here::here(
    "data",
    paste0(
      "pars_random_",
      prior_spec,
      ".csv"
    )
  )
  
} else {
  
  pars_file <- here::here(
    "data",
    paste0(
      "pars_scenarios_",
      prior_spec,
      ".csv"
    )
  )
}


# -------------------------------------------------------------------------
# Calculate summaries for the selected design
# -------------------------------------------------------------------------
if (design_method == "rule_based") {
  
  sims_file <- here::here(
    "rule_based",
    "results",
    paste0(
      "rb_",
      prior_spec,
      "_",
      generation_method,
      "_",
      cap_method,
      ".csv"
    )
  )
  
  simulation_summary <- sum_rb(
    sims_file = sims_file,
    pars_file = pars_file,
    generation_method = generation_method,
    dose_grid = dose_grid
  )
  
} else if (design_method == "bcrm") {
  
  results_folder <- here::here(
    "bayesian_crm",
    "results",
    paste0(
      "bcrm_",
      prior_spec,
      "_",
      generation_method,
      "_",
      cap_method
    )
  )
  
  # results_folder <- here::here(
  #   "extra_RESULTS",
  #   "bcrm_results",
  #   paste0(
  #     "bcrm_",
  #     prior_spec,
  #     "_",
  #     generation_method,
  #     "_",
  #     cap_method
  #   )
  # )
  
  simulation_summary <- sum_bcrm(
    sims_file = file.path(
      results_folder,
      "sims.csv"
    ),
    sumstats_file = file.path(
      results_folder,
      "summary_stats.csv"
    ),
    est_id_file = file.path(
      results_folder,
      "est.id70s.csv"
    ),
    pars_file = pars_file,
    generation_method = generation_method,
    dose_grid = dose_grid
  )
}


# Display summary results
simulation_summary


}





# Obtain summaries for both designs under capped and uncapped settings

get_all_simulation_summaries <- function(
    generation_method,
    prior_spec,
    dose_grid,
    cap_methods = c("uncapped", "capped")
) {
  
  # Check the requested cap methods
  cap_methods <- match.arg(
    cap_methods,
    choices = c("uncapped", "capped"),
    several.ok = TRUE
  )
  
  # Calculate rule-based and Bayesian CRM summaries
  # for each requested cap method
  summaries <- lapply(
    cap_methods,
    function(current_cap_method) {
      
      list(
        rule_based = get_simulation_summary(
          design_method = "rule_based",
          cap_method = current_cap_method,
          generation_method = generation_method,
          prior_spec = prior_spec,
          dose_grid = dose_grid
        ),
        
        bcrm = get_simulation_summary(
          design_method = "bcrm",
          cap_method = current_cap_method,
          generation_method = generation_method,
          prior_spec = prior_spec,
          dose_grid = dose_grid
        )
      )
    }
  )
  
  # Name each list element "capped" or "uncapped"
  names(summaries) <- cap_methods
  
  return(summaries)
}