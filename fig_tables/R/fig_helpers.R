

# ======================================================================================================
# Fig 1: Conceptual illustration of a monotonic dose-response relationship used in dose-finding studies
# ======================================================================================================
create_fig_concept <- function(dose, TTL, target_inf) {

  # Logistic-type curve for illustration
  tox_prob <- 1 / (1 + exp(-(dose - 5)))
  inf_prob <- tox_prob

  # Approximate target doses
  MTD <- dose[which.min(abs(tox_prob - TTL))]
  ID70 <- dose[which.min(abs(inf_prob - target_inf))]

  df1 <- data.frame(dose = dose, prob = tox_prob)
  df2 <- data.frame(dose = dose, prob = inf_prob)

  fig <- ggplot(df2, aes(dose, prob)) +
    geom_line(linewidth = 1) +
    annotate("segment", x = 0,    xend = ID70, y = target_inf, yend = target_inf,
           linetype = "dotted", color = "blue") +
    annotate("segment", x = ID70, xend = ID70, y = 0, yend = target_inf,
           linetype = "dotted", color = "blue") +
    scale_x_continuous(
      breaks = ID70,
      labels = expression("MTD or " * italic(ID[p]))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = target_inf,
      labels = "Target DLT \n probability \n or \n Target \n infection \n probability (p)"
    ) +
    labs(
      x = "Dose level",
      y = "Probability of an outcome"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank()
    )


  fig

}




# =========================================================================
# Fig 2: Example dose-response curves under several commonly used models
# =========================================================================

create_fig_ds_curves <- function() {

  # -------------------------------------------------
  # CFU dose grid
  # -------------------------------------------------
  dose_cfu <- 10^seq(3, 7, length.out = 400)   # 10^3 to 10^7

  # Standardised dose scale for CRM working models
  z <- (log10(dose_cfu) - 3) / 4

  # -------------------------------------------------
  # Model definitions
  # -------------------------------------------------

  # Power model
  power_model <- function(z, beta) {
    z^(exp(beta))
  }

  # One-parameter logistic
  one_logistic <- function(z, beta) {
    eta <- -3 + exp(beta) * 6 * z
    exp(eta) / (1 + exp(eta))
  }

  # Two-parameter logistic
  two_logistic <- function(z, beta1, beta2) {
    eta <- beta1 + exp(beta2) * 6 * z
    exp(eta) / (1 + exp(eta))
  }

  # Beta-Poisson
  beta_poisson <- function(d, a, b) {
    1 - (1 + d / b)^(-a)
  }

  # -------------------------------------------------
  # One representative curve per model
  # -------------------------------------------------
  plot_df <- bind_rows(
    tibble(
      dose = dose_cfu,
      prob = power_model(z, beta = 0.5),
      model = "Power model"
    ),
    tibble(
      dose = dose_cfu,
      prob = one_logistic(z, beta = 0),
      model = "One-parameter logistic"
    ),
    tibble(
      dose = dose_cfu,
      prob = two_logistic(z, beta1 = -2, beta2 = 0),
      model = "Two-parameter logistic"
    ),
    tibble(
      dose = dose_cfu,
      prob = beta_poisson(dose_cfu, a = 0.6, b = 1e5),
      model = "Beta-Poisson"
    )
  )

  # -------------------------------------------------
  # Plot
  # -------------------------------------------------
  ggplot(plot_df, aes(x = dose, y = prob, colour = model)) +
    geom_line(linewidth = 1.2) +
    scale_x_log10(
      #labels = label_math(10^.x)
      breaks = NULL
    ) +
    scale_y_continuous(limits = c(0,1),
                     breaks = NULL) +
  
    # Nice color palette
    scale_colour_manual(values = c(
      "Power model" = "#1B9E77",
      "One-parameter logistic" = "#D95F02",
      "Two-parameter logistic" = "#7570B3",
      "Beta-Poisson" = "#E7298A"
    )) +
  
    labs(
      x = "Dose level",
      y = "Probability of an outcome",
      colour = "Dose-response model"
    ) +
  
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}




# ============================================================================================
# Fig 3: Illustration of posterior updating and adaptive dose allocation in the Bayesian CRM
# ============================================================================================





# ============================================================================================
# Fig 4: Distribution of trial sizes under the rule-based and Bayesian CRM designs (uncapped 
# and capped) across 1000 prior predictive simulations
# ============================================================================================

# ---------------------------------------------------------------------------------------
# Helper: Rule-based trial sizes
# ---------------------------------------------------------------------------------------
read_rulebased_sizes <- function(prior_spec, generation_method, cap_method) {
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
  
  sims <- readr::read_csv(sims_file, show_col_types = FALSE)
  
  sims %>%
    group_by(sim) %>%
    summarise(total_n = sum(total), .groups = "drop")
}


# -----------------------------------------------------------------------------------
# Helper: Bayesian CRM trial sizes
# -----------------------------------------------------------------------------------
read_crm_sizes <- function(prior_spec, generation_method, cap_method) {
  results_folder <- here::here(
    "bayesian_crm",
    "results",
    "sim_1000",
    paste0(
      "bcrm_",
      prior_spec,
      "_",
      generation_method,
      "_",
      cap_method
    )
  )
  
  sims_file = file.path(
    results_folder,
    "sims.csv"
  )
  
  sims <- readr::read_csv(
    sims_file, 
    col_names = c("dose", "infected", "total", "scenario", "design", "sim", "kk"), 
    show_col_types = FALSE
    )
  
  sims %>%
    group_by(sim) %>%
    summarise(total_n = sum(total), .groups = "drop")
}


# -----------------------------------------------------------------------------------
# Creating the figure
# -----------------------------------------------------------------------------------
create_fig_trial_size_random <- function(rule_cap, rule_uncap, crm_cap, crm_uncap) {

  #-----------------------------
  # Panel A: capped rule-based
  #-----------------------------
  pA <- ggplot(rule_cap, aes(x = total_n)) +
    geom_bar(aes(y = after_stat(prop))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "B. Rule-based (capped)",
      x = "Total trial size",
      y = "Proportion of simulations"
    ) +
    theme_bw(base_size = 13)

  #-----------------------------
  # Panel B: uncapped rule-based
  #-----------------------------
  pB <- ggplot(rule_uncap, aes(x = total_n)) +
    geom_bar(aes(y = after_stat(prop))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "A. Rule-based (uncapped)",
      x = "Total trial size",
      y = "Proportion of simulations"
    ) +
    theme_bw(base_size = 13)

  #-----------------------------
  # Panel C: capped Bayesian CRM
  #-----------------------------
  pC <- ggplot(crm_cap, aes(x = total_n)) +
    geom_bar(aes(y = after_stat(prop))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "D. Bayesian CRM (capped)",
      x = "Total trial size",
      y = "Proportion of simulations"
    ) +
    theme_bw(base_size = 13)

  #-----------------------------
  # Panel D: uncapped Bayesian CRM
  #-----------------------------
  pD <- ggplot(crm_uncap, aes(x = total_n)) +
    geom_bar(aes(y = after_stat(prop))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "C. Bayesian CRM (uncapped)",
      x = "Total trial size",
      y = "Proportion of simulations"
    ) +
    theme_bw(base_size = 13)


  # common x-range for B and D (uncapped)
  x_min_uncap <- min(rule_uncap$total_n, crm_uncap$total_n, na.rm = TRUE)
  x_max_uncap <- max(rule_uncap$total_n, crm_uncap$total_n, na.rm = TRUE)

  pB <- pB +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    coord_cartesian(xlim = c(x_min_uncap, x_max_uncap))

  pD <- pD +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    coord_cartesian(xlim = c(x_min_uncap, x_max_uncap))


  #-----------------------------
  # Combine into 2×2 grid
  #-----------------------------
  (pB | pA) /
    (pD | pC)

}




# ==========================================================================================================
# Fig 5: Prior-predictive dose-response relationship (Five scenarios: Prior on ID70 and a)
# ==========================================================================================================

create_fig_prior_pred_with5sc <- function(nsamp, id70_settings, ab_settings, prior_spec) {
  
  # -------------------------------------------------------------------------
  # Prior hyper parameters
  # -------------------------------------------------------------------------
  if (prior_spec == "id70_p") {

    a <- exp(runif(nsamp, log(id70_settings$a_lower), log(id70_settings$a_upper)))
    
    log10_id70 <- rtruncnorm(
      nsamp,
      a    = id70_settings$log10_id70_lower, 
      b    = id70_settings$log10_id70_upper,
      mean = id70_settings$log10_id70_mean, 
      sd   = id70_settings$log10_id70_sd
    )
    
    id70 <- 10^log10_id70
    
    b <- id70 / ((0.3)^(-1/a) - 1)

  } else{
    
    set.seed(2025)
    a <- exp(runif(nsamp, min = log(ab_settings$a_lower), max = log(ab_settings$a_upper)))
    b <- exp(runif(nsamp, min = log(ab_settings$b_lower),  max = log(ab_settings$b_upper)))
  }


  xgrid <- seq(1, 9, by = 0.01)
  dgrid <- 10^xgrid

  # compute prior curves
  prior_curves <- tibble(s = 1:nsamp, a = a, b = b) %>%
    group_by(s) %>% #for each draw from prior
    group_modify(~{ #generate the dose-response curve per draw
      tibble(
        x = xgrid,
        p = 1 - (1 + dgrid / .x$b)^(-.x$a)
      )
    }) %>%
    ungroup()

  # pointwise quantiles across prior draws
  prior_band <- prior_curves %>%
    group_by(x) %>% #for each dose,
    summarise(
      q05 = quantile(p, 0.025),
      q25 = quantile(p, 0.25),
      q50 = quantile(p, 0.50),
      q75 = quantile(p, 0.75),
      q95 = quantile(p, 0.975),
      .groups = "drop"
    )

  #Overlay 5 representative scenarios 
  pars_file <- here::here(
    "data",
    paste0(
      "pars_scenarios_",
      prior_spec,
      ".csv"
    )
  )
  
  scenarios <- readr::read_csv(pars_file, show_col_types = FALSE)
  
  scenarios5 <- scenarios %>%
    dplyr::select(scenario, a, b) %>%
    tidyr::expand_grid(x = xgrid) %>%
    dplyr::mutate(
      dose = 10^x,
      p = 1 - (1 + dose/b)^(-a)
    )

  # put small numeric labels "1..5" near each curve at its ID70
  sc_labels <- scenarios %>%
    dplyr::mutate(
      id70 = b * ((0.3)^(-1/a) - 1),
      x = log10(id70),
      y = 0.75
    ) %>%
    dplyr::select(scenario, x, y)

  ggplot(prior_band, aes(x = x)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "skyblue", alpha = 0.4) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.6) +
  geom_line(data = scenarios5,
          aes(y = p, group = factor(scenario)),
          linetype = "22", color = "black", linewidth = 0.9) +
  # optional numeric labels on the curves:
  geom_text(data = sc_labels,
            aes(x = x, y = y, label = scenario),
            inherit.aes = FALSE, size = 4) +
  geom_hline(yintercept = 0.7, color = "red") +
  geom_hline(yintercept = 0.6, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
  geom_rug(data = tibble(x = log10(c(1e3,1e4,1e5,1e6,1e7))),
           aes(x = x), sides = "b", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(1,9), breaks = 1:9) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = expression("Dose (log"[10]*" CFU)"),
       y = "Probability of infection") +
  theme_bw(base_size = 14)
  
}





# ==========================================================================================================
# Fig 6: Distribution of total study sizes across the five scenarios for the capped 5+5+10 rule-based 
# design and capped Bayesian CRM
# ==========================================================================================================
create_fig_trial_size_5scn <- function(rule_sims, crm_sims) {
  
  rule_sims <- rule_sims %>%
    mutate(
      total    = as.character(total),
      scenario = as.character(scenario)
    ) %>%
    filter(scenario != "scenario") %>%
    mutate(
      dose     = as.numeric(dose),
      infected = as.numeric(infected),
      total    = as.numeric(total),
      scenario = factor(scenario),
      sim      = as.integer(sim)
    )

  rule_study_size <- rule_sims %>%
    group_by(scenario, sim) %>%
    summarise(n_total = sum(total), .groups = "drop")



  crm_sims <- crm_sims %>%
    mutate(
      dose     = as.numeric(dose),
      infected = as.numeric(infected),
      total    = as.numeric(total),
      scenario = factor(scenario),
      sim      = as.integer(sim)
  )

  crm_study_size <- crm_sims %>%
    group_by(scenario, sim) %>%
    summarise(n_total = sum(total), .groups = "drop")

  ## Panel A: Rule-based
  p_rule <- ggplot(rule_study_size, aes(x = n_total)) +
    geom_histogram(
      binwidth = 5,
      fill  = "grey20",
      color = "white",
      width = 0.9
    ) +
    facet_wrap(~ scenario, nrow = 1) +
    labs(
      x = "Total study size",
      y = "Count",
      title = "A. Rule-based design (capped)"
    ) +
    scale_y_continuous(limits = c(0,1000)) +
    theme_bw(base_size = 14) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )

  ## Panel B: Bayesian CRM
  p_crm <- ggplot(crm_study_size, aes(x = n_total)) +
    geom_histogram(
      binwidth = 5,
      fill  = "grey20",
      color = "white",
      width = 0.9
    ) +
    facet_wrap(~ scenario, nrow = 1) +
    labs(
      x = "Total study size",
      y = "Count",
      title = "B. Bayesian CRM (capped)"
    ) +
    scale_y_continuous(limits = c(0,1000)) +
    theme_bw(base_size = 14) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )

  p_rule / p_crm 

}





# ==========================================================================================================
# Fig 7: Distribution of participants allocated to the target dose across the five scenarios for the capped 
# 5+5+10 rule-based design vs capped Bayesian CRM
# ==========================================================================================================

create_fig_target_5scn <- function(rule_sims, crm_sims) {
  # -------------------------------------
  # BAYESIAN CRM: read + process
  # -------------------------------------
  crm_sims <- crm_sims %>%
    mutate(
      dose     = as.numeric(dose),
      infected = as.numeric(infected),
      total    = as.numeric(total),
      scenario = factor(scenario),
      sim      = as.integer(sim)
    )

  # Final selected dose per simulation
  crm_final_dose <- crm_sims %>%
    group_by(scenario, sim) %>%
    summarise(final_dose = last(dose), .groups = "drop")

  # Number treated at that final (target) dose
  crm_n_target <- crm_sims %>%
    group_by(scenario, sim, dose) %>%
    summarise(n_at_dose = sum(total), .groups = "drop") %>%
    left_join(crm_final_dose, by = c("scenario", "sim")) %>%
    filter(dose == final_dose) %>%
    group_by(scenario, sim) %>%
    summarise(n_at_target = sum(n_at_dose), .groups = "drop")

  # Turn counts into % of simulations
  crm_n_target_pct <- crm_n_target %>%
    group_by(scenario, n_at_target) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ungroup()


  # -------------------------------------------------------------------------------------
  # RULE-BASED: read + process
  # How many participants were ever treated at the final dose during the trial?
  # -------------------------------------------------------------------------------------
  # Clean header row + coerce types
  rule_sims <- rule_sims %>%
    mutate(
      scenario = as.character(scenario),
      total    = as.character(total)
    ) %>%
    filter(scenario != "scenario") %>%      # drop fake header row
    mutate(
      dose     = as.numeric(dose),
      infected = as.numeric(infected),
      total    = as.numeric(total),
      scenario = factor(scenario),
      sim      = as.integer(sim)
    )

  # Final selected dose per simulation
  rule_final_dose <- rule_sims %>%
    group_by(scenario, sim) %>%
    summarise(final_dose = last(dose), .groups = "drop")

  # Number treated at that final (target) dose
  rule_n_target <- rule_sims %>%
    group_by(scenario, sim, dose) %>%
    summarise(n_at_dose = sum(total), .groups = "drop") %>%
    left_join(rule_final_dose, by = c("scenario", "sim")) %>%
    filter(dose == final_dose) %>%
    group_by(scenario, sim) %>%
    summarise(n_at_target = sum(n_at_dose), .groups = "drop")

  # Turn counts into % of simulations
  rule_n_target_pct <- rule_n_target %>%
    group_by(scenario, n_at_target) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ungroup()

  
  # -----------------------------
  # Make the two plots
  # -----------------------------

  # Panel A: rule-based
  p_rule <- ggplot(rule_n_target_pct, aes(x = n_at_target, y = pct)) +
    geom_col(fill = "grey20", width = 1.5) +
    geom_text(
      aes(
        label = ifelse(
          round(pct) == 0,
          "",
          paste0(round(pct), "%")
        )
      ),
      vjust = -0.3,
      size = 3.2,
      color = "black"
    ) +
    facet_wrap(~ scenario, nrow = 1) +
    labs(
      x = "Number of participants allocated to target dose",
      y = "Percent of simulations",
      title = "A. Rule-based design (capped)"
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )

  # Panel B: Bayesian CRM
  p_crm <- ggplot(crm_n_target_pct, aes(x = n_at_target, y = pct)) +
    geom_col(fill = "grey20", width = 1.5) +
    geom_text(
      aes(
        label = ifelse(
          round(pct) == 0,
          "",
          paste0(round(pct), "%")
        )
      ),
      vjust = -0.3,
      size = 3.2,
      color = "black"
    ) +
    facet_wrap(~ scenario, nrow = 1) +
    labs(
      x = "Number of participants allocated to target dose",
      y = "Percent of simulations",
      title = "B. Bayesian CRM design (capped)"
    ) +
    scale_x_continuous(
      breaks = c(5, 10, 15, 20),
      limits = c(4, 21)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )

  #------------------------------------------
  # Combine into one figure (Panel A & B)
  # ------------------------------------------
  combined_fig <- p_rule + p_crm +
    plot_layout(ncol = 1)

  combined_fig
}




# ==========================================================================================================
# Fig 2 in Appendix: 
# ==========================================================================================================
create_fig_prior_pred <- function(nsamp, id70_settings, ab_settings, prior_spec) {
  
  set.seed(2025)
  
  # -------------------------------------------------------------------------
  # Prior hyper parameters
  # -------------------------------------------------------------------------
  if (prior_spec == "id70_p") {
    
    a <- exp(runif(nsamp, log(id70_settings$a_lower), log(id70_settings$a_upper)))
    
    log10_id70 <- rtruncnorm(
      nsamp,
      a    = id70_settings$log10_id70_lower, 
      b    = id70_settings$log10_id70_upper,
      mean = id70_settings$log10_id70_mean, 
      sd   = id70_settings$log10_id70_sd
    )
    
    id70 <- 10^log10_id70
    
    b <- id70 / ((0.3)^(-1/a) - 1)
    
  } else{
    
    a <- exp(runif(nsamp, min = log(ab_settings$a_lower), max = log(ab_settings$a_upper)))
    b <- exp(runif(nsamp, min = log(ab_settings$b_lower),  max = log(ab_settings$b_upper)))
  }
  
  xgrid <- seq(1, 9, by = 0.01)
  dgrid <- 10^xgrid
  
  # compute prior curves
  prior_curves <- tibble(s = 1:nsamp, a = a, b = b) %>%
    group_by(s) %>% #for each draw from prior
    group_modify(~{ #generate the dose-response curve per draw
      tibble(
        x = xgrid,
        p = 1 - (1 + dgrid / .x$b)^(-.x$a)
      )
    }) %>%
    ungroup()
  
  # pointwise quantiles across prior draws
  prior_band <- prior_curves %>%
    group_by(x) %>% #for each dose,
    summarise(
      q05 = quantile(p, 0.025),
      q25 = quantile(p, 0.25),
      q50 = quantile(p, 0.50),
      q75 = quantile(p, 0.75),
      q95 = quantile(p, 0.975),
      .groups = "drop"
    )
  
  ggplot(prior_band, aes(x = x)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = "skyblue", alpha = 0.4) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.6) +
    #geom_line(aes(y = q50), color = "blue", linewidth = 1) +
    geom_hline(yintercept = 0.7, color = "red") +
    geom_hline(yintercept = 0.6, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0.8, color = "red", linetype = "dashed") +
    geom_rug(data = tibble(x = log10(c(1e3,1e4,1e5,1e6,1e7))),
             aes(x = x), sides = "b", inherit.aes = FALSE) +
    scale_x_continuous(limits = c(1,9), breaks = 1:9) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = expression("Dose (log"[10]*" CFU)"),
         y = "Probability of infection") + 
    theme_bw(base_size = 14)
  
}










# ==========================================================================================================
# Fig 3 in Appendix: 
# ==========================================================================================================

create_fig_estid70 <- function(prior_spec, generation_method, cap_method, id70, pars_file, rule_sims) {
  
# Keep only the FINAL snapshot for each simulated trial
id70_final <- id70 %>%
  group_by(scenario, sim) %>%
  filter(n.doses == max(n.doses, na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(
    scenario,
    sim = as.integer(sim),
    # 90% CI and 50% CI bounds on log10(ID70)
    lo90 = q05, hi90 = q95,
    lo50 = q25, hi50 = q75,
    med  = q50
  )

# True ID70 (log10) per scenario for the blue horizontal line

scenarios <- readr::read_csv(pars_file, show_col_types = FALSE)%>%
  mutate(scenario = as.integer(scenario)) %>%
  transmute(scenario, true_id70 = log10_id70)

id70_plotdf <- id70_final %>%
  left_join(scenarios, by = "scenario") %>%
  mutate(panel_lab = paste0("Scenario: ", scenario))


# rule-based final-dose overlay data
rule_final <- rule_sims %>%
  group_by(scenario, sim) %>%
  summarise(
    rb_final_dose = dplyr::last(dose),
    rb_stop_reason = dplyr::last(stop_reason),  # same within trial
    .groups = "drop"
  ) %>%
  mutate(
    rb_final_log10 = log10(rb_final_dose),
    rb_status = dplyr::case_when(
      rb_stop_reason == "identified" ~ "Dose identified",
      rb_stop_reason == "cap" ~ "Stopped due to cap"))

id70_plotdf2 <- id70_plotdf %>%
  left_join(rule_final, by = c("scenario", "sim"))



# # Select 100 simulations per scenario for visualisation
# set.seed(250909)
# 
# plot_simulations <- id70_plotdf %>%
#   distinct(
#     scenario,
#     sim
#   ) %>%
#   group_by(scenario) %>%
#   slice_sample(n = 100) %>%
#   ungroup()
# 
# 
# # Restrict both plotting datasets to the same selected simulations
# id70_plotdf <- id70_plotdf %>%
#   semi_join(
#     plot_simulations,
#     by = c(
#       "scenario",
#       "sim"
#     )
#   )
# 
# 
# id70_plotdf2 <- id70_plotdf2 %>%
#   semi_join(
#     plot_simulations,
#     by = c(
#       "scenario",
#       "sim"
#     )
#   )


#plot
ggplot(id70_plotdf, aes(x = sim)) +
  geom_linerange(aes(ymin = lo90, ymax = hi90), size = 0.3, colour = "grey40", alpha = 0.7) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), size = 1.0, colour = "grey40") +
  geom_point(aes(y = med), size = 0.5, colour = "grey40") +
  geom_hline(aes(yintercept = true_id70), colour = "#2C7FB8", linewidth = 0.8) +
  
  # Rule-based final dose overlay
  geom_point(
    data = id70_plotdf2,
    aes(x = sim, y = rb_final_log10, colour = rb_status),
    size = 1.5, alpha = 0.9, shape = 16, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Dose identified" = "#8E6BBE",
      "Stopped due to cap"        = "#EB6123"
    )
  ) +
  
  facet_wrap(~ panel_lab, ncol = 1) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000), limits = c(0,1000)) +
  labs(
    x = "Simulation Number",
    y = expression(log[10]*" CFU"),
    colour = "Rule-based outcome"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


}






# Create posterior-updating plot for one illustrative CRM trial

create_fig_posterior_updating <- function(
    trial_dat,
    post_curves,
    true_a,
    true_b,
    id70_settings,
    target_prob = 0.70,
    n_cohorts_to_plot = 5,
    prior_seed = 2025,
    n_prior_draws = 4000
) {
  
  # -----------------------------------------------------------------------
  # Retain cohorts displayed in the figure
  # -----------------------------------------------------------------------
  
  trial_dat <- trial_dat %>%
    dplyr::filter(
      cohort <= n_cohorts_to_plot
    )
  
  post_curves <- post_curves %>%
    dplyr::filter(
      cohort <= n_cohorts_to_plot
    )
  
  
  # -----------------------------------------------------------------------
  # Construct prior dose-response band
  # -----------------------------------------------------------------------
  
  set.seed(prior_seed)
  
  a_prior <- exp(
    stats::runif(
      n_prior_draws,
      min = log(id70_settings$a_lower),
      max = log(id70_settings$a_upper)
    )
  )
  
  log10_id70_prior <- truncnorm::rtruncnorm(
    n_prior_draws,
    a = id70_settings$log10_id70_lower,
    b = id70_settings$log10_id70_upper,
    mean = id70_settings$log10_id70_mean,
    sd = id70_settings$log10_id70_sd
  )
  
  id70_prior <- 10^log10_id70_prior
  
  b_prior <- id70_prior / (
    (1 - target_prob)^(-1 / a_prior) - 1
  )
  
  
  xgrid <- seq(
    1,
    9,
    by = 0.01
  )
  
  dgrid <- 10^xgrid
  
  
  prior_band <- purrr::map_dfr(
    seq_along(dgrid),
    function(index) {
      
      posterior_probability <- dr.model(
        dgrid[index],
        b = b_prior,
        a = a_prior
      )
      
      tibble::tibble(
        dose = dgrid[index],
        log10_dose = xgrid[index],
        q05 = stats::quantile(
          posterior_probability,
          0.05
        ),
        q25 = stats::quantile(
          posterior_probability,
          0.25
        ),
        q50 = stats::quantile(
          posterior_probability,
          0.50
        ),
        q75 = stats::quantile(
          posterior_probability,
          0.75
        ),
        q95 = stats::quantile(
          posterior_probability,
          0.95
        ),
        cohort = 0
      )
    }
  )
  
  
  # -----------------------------------------------------------------------
  # Construct true dose-response curve
  # -----------------------------------------------------------------------
  
  true_curve <- tibble::tibble(
    dose = dgrid,
    log10_dose = xgrid,
    p_true = dr.model(
      dgrid,
      b = true_b,
      a = true_a
    )
  )
  
  
  # -----------------------------------------------------------------------
  # Combine prior and posterior panels
  # -----------------------------------------------------------------------
  
  panel_levels <- 0:n_cohorts_to_plot
  
  panel_labels <- c(
    "A. Prior",
    paste0(
      LETTERS[
        2:(n_cohorts_to_plot + 1)
      ],
      ". After cohort ",
      seq_len(n_cohorts_to_plot)
    )
  )
  
  plot_df <- dplyr::bind_rows(
    prior_band,
    post_curves %>%
      dplyr::select(
        dose,
        log10_dose,
        q05,
        q25,
        q50,
        q75,
        q95,
        cohort
      )
  ) %>%
    dplyr::mutate(
      panel = factor(
        cohort,
        levels = panel_levels,
        labels = panel_labels
      )
    )
  
  
  # -----------------------------------------------------------------------
  # Construct observed cohort points
  # -----------------------------------------------------------------------
  
  observed_points <- trial_dat %>%
    dplyr::mutate(
      log10_dose = log10(ds),
      p_obs = obs / nn
    )
  
  observed_points_panel <- purrr::map_dfr(
    seq_len(n_cohorts_to_plot),
    function(panel_cohort) {
      
      observed_points %>%
        dplyr::filter(
          cohort <= panel_cohort
        ) %>%
        dplyr::mutate(
          panel = panel_labels[
            panel_cohort + 1
          ]
        )
    }
  ) %>%
    dplyr::mutate(
      panel = factor(
        panel,
        levels = panel_labels
      )
    )
  
  
  # Adjust label locations
  observed_points_panel <- observed_points_panel %>%
    dplyr::mutate(
      x_lab = dplyr::case_when(
        cohort == 1 ~ log10_dose - 0.60,
        cohort == 3 ~ log10_dose - 0.70,
        cohort == 4 ~ log10_dose + 0.60,
        cohort == 5 ~ log10_dose + 0.60,
        TRUE ~ log10_dose
      ),
      
      y_lab = dplyr::case_when(
        cohort == 1 ~ p_obs,
        cohort == 4 ~ p_obs + 0.01,
        cohort == 5 ~ p_obs + 0.01,
        TRUE ~ pmin(
          p_obs + 0.12,
          0.99
        )
      )
    )
  
  
  # -----------------------------------------------------------------------
  # Create figure
  # -----------------------------------------------------------------------
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = log10_dose
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = q05,
        ymax = q95
      ),
      fill = "skyblue",
      alpha = 0.35
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = q25,
        ymax = q75
      ),
      fill = "steelblue",
      alpha = 0.55
    ) +
    ggplot2::geom_line(
      data = prior_band,
      ggplot2::aes(
        x = log10_dose,
        y = q05
      ),
      inherit.aes = FALSE,
      colour = "grey30",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::geom_line(
      data = prior_band,
      ggplot2::aes(
        x = log10_dose,
        y = q95
      ),
      inherit.aes = FALSE,
      colour = "grey30",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::geom_line(
      data = true_curve,
      ggplot2::aes(
        x = log10_dose,
        y = p_true
      ),
      inherit.aes = FALSE,
      colour = "darkgreen",
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = target_prob,
      colour = "red",
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = c(0.6, 0.8),
      colour = "red",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_rug(
      data = tibble::tibble(
        log10_dose = log10(
          c(
            1e3,
            1e4,
            1e5,
            1e6,
            1e7
          )
        )
      ),
      ggplot2::aes(
        x = log10_dose
      ),
      inherit.aes = FALSE,
      sides = "b"
    ) +
    ggplot2::geom_point(
      data = observed_points_panel,
      ggplot2::aes(
        x = log10_dose,
        y = p_obs
      ),
      inherit.aes = FALSE,
      size = 1.6,
      colour = "black"
    ) +
    ggplot2::geom_label(
      data = observed_points_panel,
      ggplot2::aes(
        x = x_lab,
        y = y_lab,
        label = cohort
      ),
      inherit.aes = FALSE,
      size = 4.2,
      fill = "white",
      colour = "black",
      label.size = 0.2,
      label.r = grid::unit(
        0.12,
        "lines"
      )
    ) +
    ggplot2::facet_wrap(
      ~ panel,
      ncol = 3
    ) +
    ggplot2::scale_x_continuous(
      limits = c(1, 9),
      breaks = 1:9
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = expression(
        "Dose (" * log[10] * " CFU)"
      ),
      y = "Probability of infection"
    ) +
    ggplot2::theme_bw(
      base_size = 18
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 20
      ),
      axis.title.y = ggplot2::element_text(
        size = 20
      ),
      axis.text.x = ggplot2::element_text(
        size = 15
      ),
      axis.text.y = ggplot2::element_text(
        size = 15
      ),
      strip.text = ggplot2::element_text(
        size = 17,
        face = "bold",
        margin = ggplot2::margin(
          t = 3,
          b = 3
        )
      ),
      panel.grid.minor =
        ggplot2::element_blank(),
      panel.spacing.x = grid::unit(
        0.35,
        "lines"
      ),
      panel.spacing.y = grid::unit(
        0.30,
        "lines"
      ),
      plot.margin = ggplot2::margin(
        t = 1.9,
        r = 1.9,
        b = -3.8,
        l = 1.9
      )
    )
}






