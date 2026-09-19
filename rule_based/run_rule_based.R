# Run rule-based dose-escalation simulations


# -------------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------------
library(here)


# -------------------------------------------------------------------------
# Load supporting functions
# -------------------------------------------------------------------------
source(here::here("rule_based","R","helpers_crm.R"))
source(here::here("rule_based","R","rule_based.R"))
source(here::here("rule_based","R","rb_simulation.R"))


# -------------------------------------------------------------------------
# User settings
# -------------------------------------------------------------------------

# Select: "capped" or "uncapped"
cap_method <- "capped"

# Maximum total sample size for the capped design
# This setting is ignored when cap_method = "uncapped"
n_cap <- 35

# Select: "random" or "scenarios"
generation_method <- "scenarios"

# Select: "id70_p" or "ab_p"
prior_spec <- "ab_p"

# Number of simulations per scenario
# This setting is used only when generation_method = "scenarios"
n_sim_scenario <- 1000

# Base seed used for the simulations
simulation_seed <- 250909

# -------------------------------------------------------------------------
# Run simulations
# -------------------------------------------------------------------------

rule_based_results <- run_rule_based(
  cap_method = cap_method,
  generation_method = generation_method,
  prior_spec = prior_spec,
  n_sim_scenario = n_sim_scenario,
  simulation_seed = simulation_seed,
  n_cap = n_cap
)

# -------------------------------------------------------------------------
# Save simulation results
# -------------------------------------------------------------------------

output_file <- paste0(
  "rb_",
  prior_spec,
  "_",
  generation_method,
  "_",
  cap_method,
  ".csv"
)

output_path <- here::here(
  "rule_based",
  "results",
  output_file
)

dir.create(
  dirname(output_path),
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  rule_based_results,
  output_path
)

message(
  "Rule-based simulation results saved to: ",
  output_path
)