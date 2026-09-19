# Generate true dose-response parameter values used in the simulations

library(here)
library(readr)

# Load the parameter-generation functions
source(
  here::here(
    "data_gen",
    "param_helpers.R"
  )
)

# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------
# Select the parameter-generation specification: "id70_p" or "ab_p"
prior_spec <- "ab_p"

# Select the output type: "random" or "scenarios"
generation_method <- "scenarios"


# Target infection probability used to define ID70
target_prob <- 0.70

# Number of parameter sets for the random option
n_sim <- 1000

# Number of initial draws used to select the quantile scenarios
n_prior_draws <- 100000

# Percentiles used for the quantile-scenarios option
scenario_percent <- c(0.05, 0.25, 0.50, 0.75, 0.95)


# Distribution settings for a and log10(ID70)
#
# log(a) ~ Uniform(log(a_lower), log(a_upper))
# log10(ID70) ~ Normal(log10_id70_mean, log10_id70_sd^2), truncated to (log10_id70_lower, log10_id70_upper)

id70_settings <- list(
  a_lower = 0.3,
  a_upper = 0.8,
  log10_id70_mean = 4,
  log10_id70_sd = 0.8,
  log10_id70_lower = 3,
  log10_id70_upper = 7
)


# Distribution settings for a and b
#
# log(a) ~ Uniform(log(a_lower), log(a_upper))
# log(b) ~ Uniform(log(b_lower), log(b_upper))

ab_settings <- list(
  a_lower = 0.1,
  a_upper = 0.7,
  b_lower = 20,
  b_upper = 600
)

# -------------------------------------------------------------------------
# Generate parameter values
# -------------------------------------------------------------------------

if (generation_method == "random") {
  
  set.seed(1)
  
  # Generate independent parameter sets for the simulation replicates
  simulation_param <- generate_prior_draws(
    n_draws = n_sim,
    prior_spec = prior_spec,
    target_prob = target_prob,
    id70_settings = id70_settings,
    ab_settings = ab_settings
  )
  
} else if (generation_method == "scenarios") {
  
  set.seed(2025)
  
  # Generate a large sample from the selected parameter distributions
  param_draws <- generate_prior_draws(
    n_draws = n_prior_draws,
    prior_spec = prior_spec,
    target_prob = target_prob,
    id70_settings = id70_settings,
    ab_settings = ab_settings
  )
  
  # Select representative scenarios at the requested ID70 percentiles
  simulation_param <- select_quantile_scenarios(
    param_draws = param_draws,
    scenario_percent = scenario_percent
  )
  
} else {
  
  stop(
    "`generation_method` must be either ",
    "'random' or 'scenarios'."
  )
}

# -------------------------------------------------------------------------
# Save and display the generated parameter values
# -------------------------------------------------------------------------

output_file <- paste0(
  "pars_",
  generation_method,
  "_",
  prior_spec,
  ".csv"
)

readr::write_csv(
  simulation_param,
  here::here("data", output_file)
)

print(simulation_param)

message(
  "Parameter values saved to: ",
  here::here("data", output_file)
)