# Functions for generating 'true' dose-response parameter values


# Generate draws from the parameter distributions
#
# Generates values of the dose-response parameters a and b using either:
# 1. prior distributions specified for a and log10(ID70); or
# 2. independent distributions specified for a and b.

generate_prior_draws <- function(
    n_draws,
    prior_spec,
    target_prob,
    id70_settings,
    ab_settings
) {
  
  prior_specification <- match.arg(
    prior_spec,
    choices = c("id70_p", "ab_p")
  )
  
  if (prior_spec == "id70_p") {
    
    # Generate a from a log-uniform distribution
    a_values <- exp(
      stats::runif(
        n_draws,
        min = log(id70_settings$a_lower),
        max = log(id70_settings$a_upper)
      )
    )
    
    # Generate log10(ID70) from a truncated normal distribution
    log10_id70_values <- truncnorm::rtruncnorm(
      n_draws,
      a = id70_settings$log10_id70_lower,
      b = id70_settings$log10_id70_upper,
      mean = id70_settings$log10_id70_mean,
      sd = id70_settings$log10_id70_sd
    )
    
    # Transform log10(ID70) to the original dose scale
    id70_values <- 10^log10_id70_values
    
    # Derive b using a, ID70, and the target probability
    b_values <- id70_values / (
      (1 - target_prob)^(-1 / a_values) - 1
    )
    
  } else if (prior_spec == "ab_p") {
    
    # Generate b first to reproduce the original simulation values
    b_values <- exp(
      stats::runif(
        n_draws,
        min = log(ab_settings$b_lower),
        max = log(ab_settings$b_upper)
      )
    )
    
    # Generate a second to reproduce the original simulation values
    a_values <- exp(
      stats::runif(
        n_draws,
        min = log(ab_settings$a_lower),
        max = log(ab_settings$a_upper)
      )
    )
    
    # Derive ID70 using a, b, and the target probability
    id70_values <- b_values * (
      (1 - target_prob)^(-1 / a_values) - 1
    )
    
    log10_id70_values <- log10(id70_values)
  }
  
  tibble::tibble(
    sim = seq_len(n_draws),
    a = a_values,
    b = b_values,
    ID70 = id70_values,
    ID70_log10 = log10_id70_values
  )
}


# Select scenarios at specified ID70 percentiles
#
# For each requested percentile, selects the prior draw with a value of
# log10(ID70) closest to the corresponding quantile.

select_quantile_scenarios <- function(
    param_draws,
    scenario_percent
) {
  
  if (any(scenario_percent <= 0 | scenario_percent >= 1)) {
    stop("All scenario percentiles must be strictly between 0 and 1.")
  }
  
  target_quantiles <- stats::quantile(
    param_draws$ID70_log10,
    probs = scenario_percent,
    names = FALSE
  )
  
  purrr::map_dfr(
    seq_along(target_quantiles),
    function(scenario_index) {
      
      target_quantile <- target_quantiles[scenario_index]
      
      selected_index <- which.min(
        abs(param_draws$ID70_log10 - target_quantile)
      )
      
      param_draws[selected_index, ] |>
        dplyr::transmute(
          scenario = scenario_index,
          percentile = scenario_percent[scenario_index],
          a = a,
          b = b,
          ID70 = ID70,
          log10_id70 = ID70_log10
        )
    }
  )
}