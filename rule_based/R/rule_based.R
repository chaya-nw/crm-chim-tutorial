# Simulate one trial using the rule-based dose-escalation design

rulebased_path <- function(
    b,
    a,
    scenario,
    sim,
    cap_method,
    N_cap
) {
  
  cap_method <- match.arg(
    cap_method,
    choices = c("capped", "uncapped")
  )
  
  # -----------------------------------------------------------------------
  # Design settings
  # -----------------------------------------------------------------------
  
  # Starting values
  dose <- 4 #log10 dose level; start at 10^4
  n1   <- 5
  n2   <- 10
  
  # Restrict doses to the prespecified log10 dose range
  clamp_dose <- function(dose) {
    pmin(pmax(dose, 3), 7)
  }
  
  
  # -----------------------------------------------------------------------
  # Initialise trial information
  # -----------------------------------------------------------------------
  
  total_N     <- 0
  identified  <- FALSE
  stopped     <- FALSE
  stop_reason <- NA_character_
  
  # Store cohort-level trial information
  dose_idx_vec <- numeric(0)  # log10 dose administered
  n_vec        <- integer(0)  # cohort sizes
  inf_vec      <- integer(0)  # infected counts
  
  
  # -----------------------------------------------------------------------
  # Run the trial
  # -----------------------------------------------------------------------
  
  while (!identified && !stopped) {
    
    # Reset infection counts when starting at a dose
    xA <- 0
    xB <- 0
    xC <- 0
    
    # -----------------------------------------------------------------------------------------------------
    # Step 1: challenge five participants at the current dose; ; if <2/5, escalate and repeat 
    # -----------------------------------------------------------------------------------------------------
    
    repeat {
      
      # Stop if the next cohort would exceed the sample-size cap
      if (
        cap_method == "capped" &&
        total_N + n1 > N_cap
      ) {
        stopped <- TRUE
        stop_reason <- "cap"
        break
      }
      
      p.inf <- dr.model(10^dose,b,a)
      xA    <- rbinom(1, size = n1, prob = p.inf)
      
      dose_idx_vec <- c(dose_idx_vec, dose)
      n_vec        <- c(n_vec, n1)
      inf_vec      <- c(inf_vec, xA)
      total_N      <- total_N + n1
      
      # Continue at this dose if at least 2/5 are infected
      if (xA >= 2) {
        break
      }
      
      # Otherwise, escalate by one dose level
      dose <- clamp_dose(dose + 1)     # NOTE: constrained escalation
    }
    
    if (stopped) {
      break
    }
    
    # Recalculate the infection probability at the current dose
    p.inf <- dr.model(10^dose,b,a)
    
    
    
    # ---------------------------------------------------------------------
    # Step 2: challenge another five participants at the same dose
    # ---------------------------------------------------------------------
    if (
      cap_method == "capped" &&
      total_N + n1 > N_cap
    ) {
      stopped <- TRUE
      stop_reason <- "cap" #sample size cap
      break
    }
    
    xB <- rbinom(1, size = n1, prob = p.inf)
  
    dose_idx_vec <- c(dose_idx_vec, dose)
    n_vec        <- c(n_vec, n1)
    inf_vec      <- c(inf_vec, xB)
    total_N      <- total_N + n1
    
    
    
    # ---------------------------------------------------------------------
    # Make the decision after 10 participants at the dose
    # ---------------------------------------------------------------------
    
    if (xA + xB < 6) {
      dose <- clamp_dose(dose + 1)
      next
    }
    
    
    # ---------------------------------------------------------------------
    # Step 3: challenge 10 additional participants at the same dose
    # ---------------------------------------------------------------------
    
    if (
      cap_method == "capped" &&
      total_N + n2 > N_cap
    ) {
      stopped <- TRUE
      stop_reason <- "cap"
      break
    }
    
    xC <- rbinom(1, size = n2, prob = p.inf)
    
    dose_idx_vec <- c(dose_idx_vec, dose)
    n_vec        <- c(n_vec, n2)
    inf_vec      <- c(inf_vec, xC)
    total_N      <- total_N + n2
    
    
    
    # ---------------------------------------------------------------------
    # Make the decision after 20 participants at the dose
    # ---------------------------------------------------------------------
    total_infected_20 <- xA + xB + xC
    
    if (
      total_infected_20 >= 12 &&
      total_infected_20 <= 16
    ) {
      
      identified  <- TRUE
      stop_reason <- "identified" #dose identified
      
    } else if (total_infected_20 < 12) {
      dose <- clamp_dose(dose + 1)
      next
      
    } else {
      dose <- clamp_dose(dose - 1)
      next
    }
  }
  
  # -----------------------------------------------------------------------
  # Return the simulated trial path
  # -----------------------------------------------------------------------
  
  tibble::tibble(
    dose        = 10^dose_idx_vec, #CFU
    infected    = inf_vec,
    total       = n_vec,
    scenario    = scenario,
    design      = "rule_based",
    sim         = sim,
    kk          = seq_along(dose_idx_vec),
    stop_reason = stop_reason,
    cap_method  = cap_method
  )
}