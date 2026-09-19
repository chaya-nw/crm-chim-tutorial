# Helper functions for the Bayesian CRM dose-finding simulations



# Evaluate the probability of infection under the dose-response model
#
# Arguments:
# DD: administered dose, which may be a single value or a vector
# b:  dose-response scale parameter
# a:  dose-response shape parameter
#
# Returns the infection probability:
# p(DD) = 1 - (1 + DD / b)^(-a)
dr.model <- function(DD, b, a){
  pp <- 1 - (1 + DD/b)^-a  
  names(pp) <- DD
  return(pp)
}



# Simulate one cohort at the current dose
#
# The function:
# 1. calculates the true probability of infection at the current dose using the dose-response parameters a and b;
# 2. simulates the number of infected participants in the cohort from a binomial distribution; and
# 3. appends the new cohort results to the accumulated trial data.
#
# Arguments:
# dat:            accumulated trial data from previously simulated cohorts
# scenario:       true dose-response parameter values a and b
# new.dose:       dose administered to the new cohort
# new.group.size: number of participants in the new cohort
#
# Returns a tibble with:
# ds:  administered dose
# obs: number of infected participants
# nn:  total number of participants in the cohort
my.sim.fn <- function(dat, scenario, new.dose, new.group.size){
  
  p.inf <- dr.model(new.dose, scenario$b, scenario$a)                            # calculate the true probability of infection at the current dose
  new.dat <- rbinom(1, size = new.group.size, prob = p.inf)                      # Simulate the number of infected participants in the new cohort
  out <- rbind(dat, c(new.dose, new.dat, new.group.size))                        # Append the new cohort results to the accumulated trial data
  colnames(out) <- c("ds","obs","nn")                                            # Assign variable names to the accumulated trial data
  return(as_tibble(out))
}




# Apply the CRM decision rules for dose allocation
#
# The function uses posterior infection-probability quantiles at the most
# recently administered dose to determine whether the next cohort should
# receive a higher, lower, or unchanged dose.
#
# The proposed dose is restricted to the allowable range from 10^3 to
# 10^7 CFU.
#
# Arguments:
# dat:           accumulated cohort-level trial data
# summary.stats: posterior infection-probability quantiles at the current dose
#
# Returns:
# target.n: total number of participants treated at the current dose
# new.dose: proposed next dose on the log10 scale

decision.criteria <- function(dat, summary.stats){
  
  target.dose <- dat %>% dplyr::slice(nrow(dat)) %>% pull(ds)                    # Identify the most recently administered dose
  target.n <- dat %>% filter(ds == target.dose) %>% summarise(nn = sum(nn))      # Calculate the total number of participants treated at this dose
  
  #Reads posterior quantiles of infection probability at current dose, and determine the proposed dose for the next cohort
  new.dose <- case_when(
    #If even the 25th percentile < 0.6 --> likely under-dosing --> escalate +1 log10 (×10 dose)
    summary.stats$`25%` < 0.6 ~ log10(target.dose) + 1,
    
    #Else if even the 5th percentile > 0.8 --> very likely over-dosing --> de-escalate −1 log10 (÷10 dose)
    summary.stats$`5%` > 0.8 ~ log10(target.dose) - 1,
    
    #Else → stay at the same log10 dose.
    TRUE ~ log10(target.dose)
  )
  
  # Restrict the proposed dose to the allowable log10 dose range
  if (new.dose < 3){
    new.dose <- 3
  } else if (new.dose > 7){
    new.dose <- 7
  } 
  
  return(c("target.n" = as.numeric(target.n), "new.dose" = new.dose))
}





# Summarise two-dose-level escalations in Bayesian CRM simulations

summarise_two_level_esc <- function(
    sims_file,
    tolerance = 1e-8
) {
  
  # -----------------------------------------------------------------------
  # Read cohort-level simulation results
  # -----------------------------------------------------------------------
  
  sims <- readr::read_csv(
    sims_file,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  colnames(sims) <- c(
    "ds",
    "obs",
    "nn",
    "scenario",
    "design",
    "sim",
    "kk"
  )
  
  
  # -----------------------------------------------------------------------
  # Identify dose changes between consecutive cohorts
  # -----------------------------------------------------------------------
  
  escalation_data <- sims %>%
    
    # Preserve the original order of rows in the results file
    dplyr::mutate(
      file_row = dplyr::row_number(),
      dose_log10 = log10(ds)
    ) %>%
    
    # Examine each simulated trial separately
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    
    # Arrange cohorts in the order in which they occurred
    dplyr::arrange(
      file_row,
      .by_group = TRUE
    ) %>%
    
    dplyr::mutate(
      cohort = dplyr::row_number(),
      
      # Dose given to the preceding cohort
      previous_dose_log10 = dplyr::lag(dose_log10),
      
      # Change from the preceding cohort's dose
      dose_change_log10 =
        dose_log10 - previous_dose_log10,
      
      # TRUE when the dose increased by exactly two log10 levels
      two_level_escalation =
        !is.na(dose_change_log10) &
        abs(dose_change_log10 - 2) < tolerance
    ) %>%
    
    dplyr::ungroup()
  
  
  # -----------------------------------------------------------------------
  # Trial-level indicators
  # -----------------------------------------------------------------------
  
  escalation_by_trial <- escalation_data %>%
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    
    dplyr::summarise(
      # Whether the trial contained at least one two-level escalation
      had_two_level_escalation = any(
        two_level_escalation
      ),
      
      # Number of two-level escalations within the trial
      n_two_level_escalations = sum(
        two_level_escalation
      ),
      
      .groups = "drop"
    )
  
  
  # -----------------------------------------------------------------------
  # Overall trial-level percentage
  # -----------------------------------------------------------------------
  
  overall_summary <- escalation_by_trial %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      
      n_trials_with_two_level_escalation =
        sum(had_two_level_escalation),
      
      pct_trials_with_two_level_escalation =
        100 * mean(had_two_level_escalation),
      
      n_two_level_escalation_occasions =
        sum(n_two_level_escalations)
    )
  
  
  # -----------------------------------------------------------------------
  # Trial-level percentage by scenario
  # -----------------------------------------------------------------------
  
  summary_by_scenario <- escalation_by_trial %>%
    dplyr::group_by(scenario) %>%
    
    dplyr::summarise(
      n_trials = dplyr::n(),
      
      n_trials_with_two_level_escalation =
        sum(had_two_level_escalation),
      
      pct_trials_with_two_level_escalation =
        100 * mean(had_two_level_escalation),
      
      n_two_level_escalation_occasions =
        sum(n_two_level_escalations),
      
      .groups = "drop"
    )
  
  
  # -----------------------------------------------------------------------
  # Extract the individual escalation occasions for checking
  # -----------------------------------------------------------------------
  
  escalation_occasions <- escalation_data %>%
    dplyr::filter(
      two_level_escalation
    ) %>%
    
    dplyr::transmute(
      scenario,
      design,
      sim,
      cohort,
      previous_dose_log10,
      new_dose_log10 = dose_log10,
      dose_change_log10
    ) %>%
    
    dplyr::arrange(
      scenario,
      design,
      sim,
      cohort
    )
  
  
  # -----------------------------------------------------------------------
  # Return results
  # -----------------------------------------------------------------------
  
  list(
    overall = overall_summary,
    by_scenario = summary_by_scenario,
    by_trial = escalation_by_trial,
    occasions = escalation_occasions
  )
}









# Assess cohort-level coherence of Bayesian CRM dose decisions

summarise_crm_coherence <- function(
    sims_file,
    target_prob = 0.70
) {
  
  # -----------------------------------------------------------------------
  # 1. Read cohort-level simulation results
  # -----------------------------------------------------------------------
  
  sims <- readr::read_csv(
    sims_file,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  colnames(sims) <- c(
    "ds",
    "obs",
    "nn",
    "scenario",
    "design",
    "sim",
    "kk"
  )
  
  
  # -----------------------------------------------------------------------
  # 2. Identify the next dose following each completed cohort
  # -----------------------------------------------------------------------
  
  cohort_transitions <- sims %>%
    
    # Preserve the order in which rows were written to the file
    dplyr::mutate(
      file_row = dplyr::row_number()
    ) %>%
    
    # Examine each simulated trial separately
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    
    dplyr::arrange(
      file_row,
      .by_group = TRUE
    ) %>%
    
    dplyr::mutate(
      # Cohort number within the simulated trial
      cohort = dplyr::row_number(),
      
      # Infection proportion in the current cohort
      observed_rate = obs / nn,
      
      # Dose administered to the following cohort
      next_dose = dplyr::lead(ds),
      
      # Current and next doses on the log10 scale
      current_log10_dose = log10(ds),
      next_log10_dose = log10(next_dose),
      
      # Direction of the implemented dose decision
      dose_direction = dplyr::case_when(
        is.na(next_dose) ~ NA_character_,
        next_dose > ds ~ "Escalation",
        next_dose < ds ~ "De-escalation",
        TRUE ~ "Stay"
      )
    ) %>%
    
    dplyr::ungroup() %>%
    
    # Retain only cohorts for which the next administered dose is observable
    dplyr::filter(
      !is.na(next_dose)
    )
  
  
  # -----------------------------------------------------------------------
  # 3. Identify formally incoherent dose transitions
  # -----------------------------------------------------------------------
  
  cohort_transitions <- cohort_transitions %>%
    dplyr::mutate(
      
      # De-escalation after a below-target outcome
      incoherent_deescalation =
        observed_rate < target_prob &
        next_dose < ds,
      
      # Escalation after an above-target outcome
      incoherent_escalation =
        observed_rate > target_prob &
        next_dose > ds,
      
      # Either type of incoherent transition
      incoherent_decision =
        incoherent_deescalation |
        incoherent_escalation,
      
      # Description of the type of incoherence
      incoherence_type = dplyr::case_when(
        incoherent_deescalation ~
          "De-escalation after below-target outcome",
        
        incoherent_escalation ~
          "Escalation after above-target outcome",
        
        TRUE ~ NA_character_
      )
    )
  
  
  # -----------------------------------------------------------------------
  # 4. Trial-level coherence indicators
  # -----------------------------------------------------------------------
  
  all_trials <- sims %>%
    dplyr::distinct(
      scenario,
      design,
      sim
    )
  
  
  coherence_by_trial <- cohort_transitions %>%
    dplyr::group_by(
      scenario,
      design,
      sim
    ) %>%
    
    dplyr::summarise(
      # Number of observable dose transitions in the trial
      n_transitions = dplyr::n(),
      
      # Number of incoherent transitions in the trial
      n_incoherent_decisions =
        sum(incoherent_decision),
      
      # Whether the trial had at least one incoherent de-escalation
      had_incoherent_deescalation =
        any(incoherent_deescalation),
      
      # Whether the trial had at least one incoherent escalation
      had_incoherent_escalation =
        any(incoherent_escalation),
      
      # Whether the trial had either type
      had_any_incoherent_decision =
        any(incoherent_decision),
      
      .groups = "drop"
    ) %>%
    
    # Retain all simulated trials
    dplyr::right_join(
      all_trials,
      by = c(
        "scenario",
        "design",
        "sim"
      )
    ) %>%
    
    # Trials with no observable transitions are assigned zero/false
    dplyr::mutate(
      n_transitions = tidyr::replace_na(
        n_transitions,
        0L
      ),
      
      n_incoherent_decisions = tidyr::replace_na(
        n_incoherent_decisions,
        0L
      ),
      
      had_incoherent_deescalation = tidyr::replace_na(
        had_incoherent_deescalation,
        FALSE
      ),
      
      had_incoherent_escalation = tidyr::replace_na(
        had_incoherent_escalation,
        FALSE
      ),
      
      had_any_incoherent_decision = tidyr::replace_na(
        had_any_incoherent_decision,
        FALSE
      )
    ) %>%
    
    dplyr::arrange(
      scenario,
      design,
      sim
    )
  
  
  # -----------------------------------------------------------------------
  # 5. Overall trial-level summary
  # -----------------------------------------------------------------------
  
  overall_by_trial <- coherence_by_trial %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      
      n_trials_with_incoherent_deescalation =
        sum(had_incoherent_deescalation),
      
      pct_trials_with_incoherent_deescalation =
        100 * mean(had_incoherent_deescalation),
      
      n_trials_with_incoherent_escalation =
        sum(had_incoherent_escalation),
      
      pct_trials_with_incoherent_escalation =
        100 * mean(had_incoherent_escalation),
      
      n_trials_with_any_incoherent_decision =
        sum(had_any_incoherent_decision),
      
      pct_trials_with_any_incoherent_decision =
        100 * mean(had_any_incoherent_decision)
    )
  
  
  # -----------------------------------------------------------------------
  # 6. Trial-level summary by fixed scenario
  # -----------------------------------------------------------------------
  
  summary_by_scenario <- coherence_by_trial %>%
    dplyr::group_by(scenario) %>%
    
    dplyr::summarise(
      n_trials = dplyr::n(),
      
      n_trials_with_incoherent_deescalation =
        sum(had_incoherent_deescalation),
      
      pct_trials_with_incoherent_deescalation =
        100 * mean(had_incoherent_deescalation),
      
      n_trials_with_incoherent_escalation =
        sum(had_incoherent_escalation),
      
      pct_trials_with_incoherent_escalation =
        100 * mean(had_incoherent_escalation),
      
      n_trials_with_any_incoherent_decision =
        sum(had_any_incoherent_decision),
      
      pct_trials_with_any_incoherent_decision =
        100 * mean(had_any_incoherent_decision),
      
      .groups = "drop"
    )
  
  
  # -----------------------------------------------------------------------
  # 7. Occasion-level summary
  # -----------------------------------------------------------------------
  
  overall_by_occasion <- cohort_transitions %>%
    dplyr::summarise(
      n_observable_transitions = dplyr::n(),
      
      n_incoherent_deescalations =
        sum(incoherent_deescalation),
      
      n_incoherent_escalations =
        sum(incoherent_escalation),
      
      n_incoherent_decisions =
        sum(incoherent_decision),
      
      pct_incoherent_transitions =
        100 * mean(incoherent_decision)
    )
  
  
  # -----------------------------------------------------------------------
  # 8. Extract individual incoherent decisions for checking
  # -----------------------------------------------------------------------
  
  incoherent_decisions <- cohort_transitions %>%
    dplyr::filter(
      incoherent_decision
    ) %>%
    
    dplyr::select(
      scenario,
      design,
      sim,
      cohort,
      ds,
      obs,
      nn,
      observed_rate,
      next_dose,
      current_log10_dose,
      next_log10_dose,
      dose_direction,
      incoherence_type
    ) %>%
    
    dplyr::arrange(
      scenario,
      design,
      sim,
      cohort
    )
  
  
  # -----------------------------------------------------------------------
  # 9. Return results
  # -----------------------------------------------------------------------
  
  list(
    overall_by_trial = overall_by_trial,
    summary_by_scenario = summary_by_scenario,
    overall_by_occasion = overall_by_occasion,
    coherence_by_trial = coherence_by_trial,
    cohort_transitions = cohort_transitions,
    incoherent_decisions = incoherent_decisions
  )
}

