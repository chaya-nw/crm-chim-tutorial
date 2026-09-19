# Run Bayesian CRM dose-finding simulations


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
source(here::here("bayesian_crm","R","bcrm_model_update.R"))
source(here::here("bayesian_crm","R","bcrm_simulation.R"))


# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

# Select: "capped" or "uncapped"
cap_method <- "capped"

# Maximum total sample size for the capped design
# This setting is ignored when cap_method = "uncapped"
N_cap <- 35

# Select: "random" or "scenarios"
generation_method <- "scenarios"

# Select prior: "id70_p" or "ab_p"; id70_p = Prior on a and log10(ID70); ab_p = Prior on a and b
prior_spec <- "id70_p"


# Target infection probability used to define ID70
target_prob <- 0.70


# Prior distribution settings for a and log10(ID70)
# id70_p
# log(a) ~ Uniform(log(a_lower), log(a_upper))
# log10(ID70) ~ Normal(log10_id70_mean, log10_id70_sd^2),
# truncated to (log10_id70_lower, log10_id70_upper)

id70_settings <- list(
  a_lower = 0.3,
  a_upper = 0.8,
  log10_id70_mean = 4,
  log10_id70_sd = 0.8,
  log10_id70_lower = 3,
  log10_id70_upper = 7
)


# Distribution settings for a and b
# ab_p
# log(a) ~ Uniform(log(a_lower), log(a_upper))
# log(b) ~ Uniform(log(b_lower), log(b_upper))

ab_settings <- list(
  a_lower = 0.1,
  a_upper = 0.7,
  b_lower = 20,
  b_upper = 600
)


# Number of random parameter sets
n_random <- 1000

# Number of simulations per selected scenario
n_sim_scenario <- 100

# Base simulation seed
seed.no <- 250909


# Starting dose and cohort size
designs <- tibble(
  start.dose = c(10^4),
  group.size = c(5)
)


# Doses that can be administered
target.doses <- c(10^3, 10^4, 10^5, 10^6, 10^7)


# MCMC settings
chains <- 4
iter   <- 3000
warmup <- 1000


# -------------------------------------------------------------------------
# Select and compile the Stan model
# -------------------------------------------------------------------------

stan_file <- switch(
  prior_spec,
  id70_p = "crm_dr_id70_prior.stan",
  ab_p = "crm_dr_ab_prior.stan"
)

compiled_model <- stan_model(here::here("bayesian_crm", "stan", stan_file))



# -------------------------------------------------------------------------
# Create the output folder
# -------------------------------------------------------------------------
out.folder <- here::here(
  "fig_tables",
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

dir.create(out.folder, recursive = TRUE, showWarnings = FALSE)
save.output.to <- out.folder


# -------------------------------------------------------------------------
# Run simulations
# -------------------------------------------------------------------------
bcrm_results <- run_bcrm(
  cap_method = cap_method,
  generation_method = generation_method,
  prior_spec = prior_spec,
  target_prob = target_prob,
  id70_settings = id70_settings,
  ab_settings = ab_settings,
  N_cap = N_cap,
  n_random = n_random,
  n_sim_scenario = n_sim_scenario,
  seed.no = seed.no,
  designs = designs,
  target.doses = target.doses,
  compiled_model = compiled_model,
  chains = chains,
  iter = iter,
  warmup = warmup,
  save.output.to = save.output.to
)


# -------------------------------------------------------------------------
# Completion message
# -------------------------------------------------------------------------
message(
  "Bayesian CRM simulation results saved to: ",
  save.output.to
)
