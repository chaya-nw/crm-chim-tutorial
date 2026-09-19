# Generate one illustrative Bayesian CRM trial

library(here)
library(tidyverse)
library(rstan)
library(truncnorm)


# -----------------------------------------------------------------------
# Load supporting functions
# -----------------------------------------------------------------------

source(here::here("bayesian_crm", "R", "helpers_crm.R"))
source(here::here("bayesian_crm", "R", "bcrm_model_update.R"))
source(here::here("bayesian_crm", "R", "bcrm_example.R"))


# -----------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------
prior_spec <- "id70_p"
cap_method <- "capped"

# Paper example
example_scenario <- 4
example_sim <- 45

target_prob <- 0.70
N_cap <- 35
seed.no <- 250909

start.dose <- 10^4
group.size <- 5

target.doses <- 10^(3:7)


# -----------------------------------------------------------------------
# Prior settings
# -----------------------------------------------------------------------

id70_settings <- list(
  a_lower = 0.3,
  a_upper = 0.8,
  log10_id70_mean = 4,
  log10_id70_sd = 0.8,
  log10_id70_lower = 3,
  log10_id70_upper = 7
)

ab_settings <- list(
  a_lower = 0.1,
  a_upper = 0.7,
  b_lower = 20,
  b_upper = 600
)


# -----------------------------------------------------------------------
# MCMC settings
# -----------------------------------------------------------------------

chains <- 4
iter <- 3000
warmup <- 1000


# -----------------------------------------------------------------------
# Compile the selected Stan model
# -----------------------------------------------------------------------

stan_file <- switch(
  prior_spec,
  id70_p = "crm_dr_id70_prior.stan",
  ab_p = "crm_dr_ab_prior.stan"
)

compiled_model <- rstan::stan_model(
  here::here(
    "bayesian_crm",
    "stan",
    stan_file
  )
)


# -----------------------------------------------------------------------
# Output folder
# -----------------------------------------------------------------------

example_output_folder <- here::here(
  "bayesian_crm",
  "results",
  "examples",
  paste0(
    "scenario_",
    example_scenario,
    "_sim_",
    example_sim,
    "_",
    prior_spec,
    "_",
    cap_method
  )
)


# -----------------------------------------------------------------------
# Run one trial
# -----------------------------------------------------------------------

example_results <- run_bcrm_example(
  scenario_no = example_scenario,
  sim_no = example_sim,
  prior_spec = prior_spec,
  target_prob = target_prob,
  id70_settings = id70_settings,
  ab_settings = ab_settings,
  cap_method = cap_method,
  N_cap = N_cap,
  seed.no = seed.no,
  start.dose = start.dose,
  group.size = group.size,
  target.doses = target.doses,
  compiled_model = compiled_model,
  chains = chains,
  iter = iter,
  warmup = warmup,
  save.output.to = example_output_folder
)


# -----------------------------------------------------------------------
# Inspect the illustrative trial
# -----------------------------------------------------------------------

example_results$trial_data
example_results$true_scenario
example_results$final_administered_dose
example_results$final_recommended_dose
example_results$n_cohorts
example_results$total_N