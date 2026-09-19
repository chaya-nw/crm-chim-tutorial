# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(here)
library(tidyverse)
library(rstan)
library(truncnorm)


# -------------------------------------------------------------------------
# Load supporting functions
# -------------------------------------------------------------------------
source(here::here("bayesian_crm","R","helpers_crm.R"))



# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

# Select: "capped" or "uncapped"
# cap_method <- "capped"

# Maximum total sample size for the capped design
# This setting is ignored when cap_method = "uncapped"
# N_cap <- 35

# Select: "random" or "scenarios"
# generation_method <- "scenarios"

# Select prior: "id70_p" or "ab_p"; id70_p = Prior on a and log10(ID70); ab_p = Prior on a and b
prior_spec <- "id70_p"



#Prior-predictive simulations
prior_escalation_results <-
  summarise_two_level_esc(
    sims_file = here::here(
      "bayesian_crm",
      "results",
      paste0(
        "bcrm_",
        prior_spec,
        "_random_capped"
      ),
      "sims.csv"
    )
  )

prior_escalation_results$overall



#Fixed-scenario simulations
scenario_escalation_results <-
  summarise_two_level_esc(
    sims_file = here::here(
      "bayesian_crm",
      "results",
      paste0(
        "bcrm_",
        prior_spec,
        "_scenarios_capped"
      ),
      "sims.csv"
    )
  )

scenario_escalation_results$by_scenario


# scenario_5_escalation <-
#   scenario_escalation_results$by_scenario %>%
#   dplyr::filter(
#     scenario == 5
#   )
# 
# scenario_5_escalation

#prior_escalation_results$occasions





#####################################random
prior_spec <- "id70_p"

prior_coherence_results <- summarise_crm_coherence(
  sims_file = here::here(
    "bayesian_crm",
    "results",
    paste0(
      "bcrm_",
      prior_spec,
      "_random_capped"
    ),
    "sims.csv"
  ),
  target_prob = 0.70
)

prior_coherence_results$overall_by_trial
prior_coherence_results$overall_by_occasion
prior_coherence_results$incoherent_decisions

#####################################fiveS
scenario_coherence_results <- summarise_crm_coherence(
  sims_file = here::here(
    "bayesian_crm",
    "results",
    paste0(
      "bcrm_",
      prior_spec,
      "_scenarios_capped"
    ),
    "sims.csv"
  ),
  target_prob = 0.70
)

scenario_coherence_results$summary_by_scenario
prior_coherence_results$overall_by_trial

