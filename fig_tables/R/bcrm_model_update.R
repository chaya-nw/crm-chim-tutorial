# Fit the Bayesian CRM model at one interim analysis
#
# This function updates the CRM using all accumulated cohort data, obtains posterior draws and summary statistics, 
# identifies the model-estimated ID70 dose category, and returns the posterior infection-probability summary 
# for the current dose.

# Arguments:
# dat:            accumulated cohort-level trial data containing dose,
#                 infected count, and cohort size
# scenario.no:    index of the true parameter set or scenario
# des.no:         CRM design index
# sim.no:         simulation replicate number
# target.doses:   vector of candidate doses that may be administered
# save.output.to: folder in which posterior summaries are saved
# row.numb:       overall simulation-job index
# compiled_model: compiled Stan model corresponding to the selected prior
# chains:         number of MCMC chains
# iter:           total MCMC iterations per chain, including warm-up
# warmup:         number of warm-up iterations per chain
# prior_spec:     selected prior specification, either "id70_p" or "ab_p"
# target_prob:    target infection probability used to define ID70
# id70_settings:  prior settings for a and log10(ID70)
# ab_settings:    prior settings for a and b

update_crm_model <- function(
    dat,
    scenario.no,
    des.no,
    sim.no,
    target.doses,
    save.output.to,
    row.numb,
    compiled_model,
    chains,
    iter,
    warmup,
    prior_spec,
    target_prob,
    id70_settings,
    ab_settings
  ) {
  
  # -----------------------------------------------------------------------
  # Simulation indices
  # -----------------------------------------------------------------------
  
  dd <- des.no
  i <- sim.no
  ss <- scenario.no

  # Doses administered up to the current interim analysis
  doses <- dat$ds
  
  # Number of cohorts observed so far
  j <- length(doses)

  
  
  # -----------------------------------------------------------------------
  # Prepare accumulated trial data for Stan
  # -----------------------------------------------------------------------
  
  dat <- dat %>%
    dplyr::select(ds, obs, nn) %>%
    tidyr::unnest(ds:nn)
  
  ds_vec <- as.numeric(dat$ds)
  nn_vec <- as.integer(dat$nn)
  obs_vec <- as.integer(dat$obs)
  
  N <- length(ds_vec)
  
  # Data common to both prior specifications
  data_list <- list(
    N = as.integer(N),
    ds = array(ds_vec, dim = N),
    nn = array(nn_vec, dim = N),
    obs = array(obs_vec, dim = N)
  )
  
  
  # -----------------------------------------------------------------------
  # Add settings for the selected prior specification
  # -----------------------------------------------------------------------
  
  if (prior_spec == "id70_p") {
    
    data_list$a_lower          <- as.numeric(id70_settings$a_lower)
    data_list$a_upper          <- as.numeric(id70_settings$a_upper)
    
    data_list$log10_id70_mean  <- as.numeric(id70_settings$log10_id70_mean)
    data_list$log10_id70_sd    <- as.numeric(id70_settings$log10_id70_sd)
    
    data_list$log10_id70_lower <- as.numeric(id70_settings$log10_id70_lower)
    data_list$log10_id70_upper <- as.numeric(id70_settings$log10_id70_upper)
    
    data_list$target_prob      <- as.numeric(target_prob)
    
  } else if (prior_spec == "ab_p") {
    
    data_list$a_lower          <- as.numeric(ab_settings$a_lower)
    data_list$a_upper          <- as.numeric(ab_settings$a_upper)
    
    data_list$b_lower          <- as.numeric(ab_settings$b_lower)
    data_list$b_upper          <- as.numeric(ab_settings$b_upper)
  }
  
  
  # Ensure that the storage modes match the Stan data declarations
  storage.mode(data_list$ds) <- "double"
  storage.mode(data_list$nn) <- "integer"
  storage.mode(data_list$obs) <- "integer"
  

  # -----------------------------------------------------------------------
  # Perform MCMC sampling in Stan
  # -----------------------------------------------------------------------
  fit <- sampling(
    object = compiled_model,
    data = data_list,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = 0              # Suppresses iteration updates
  )
  

  # Extract posterior samples
  post <- extract(fit)
  
  
  # -----------------------------------------------------------------------
  # Construct posterior draws of a, b, and ID70
  # -----------------------------------------------------------------------
  
  if (prior_spec == "id70_p") {
    
    # The Stan model directly samples log(a) and log10(ID70)
    long.draws <- tibble(
      a = exp(post$log_a),
      log10_id70 = post$log10_id70,
      id70 = 10^log10_id70
    ) %>%
      mutate(
        b = id70 / (
          (1 - target_prob)^(-1 / a) - 1
        )
      )
    
  } else if (prior_spec == "ab_p") {
    
    # The Stan model directly samples log(a) and log(b)
    long.draws <- tibble(
      a = exp(post$log_a),
      b = exp(post$log_b)
    ) %>%
      mutate(
        id70 = b * (
          (1 - target_prob)^(-1 / a) - 1
        ),
        log10_id70 = log10(id70)
      )
  }
  
  
  long.draws <- long.draws %>%
    dplyr::select(
      a,
      b,
      id70,
      log10_id70
    )
    
    
  # Posterior draws of ID70 on the original and log10 dose scales
  id70_vec <- long.draws$id70
  id70_log10 <- long.draws$log10_id70
    
    
    
  # -----------------------------------------------------------------------
  # Determine the model-based ID70 dose category
  # -----------------------------------------------------------------------
  # Get model-based estimate of the ID70. Bin posterior prob on log10 dose scale
  # specify target dose as the "most likely" dose to correspond to ID70
  
  # Assign each posterior draw of log10(ID70) to a dose category
  id70.cat <- long.draws %>% 
    mutate(pinf = id70_log10) %>% 
    mutate(pinf.cat = cut(pinf, breaks = c(0,2.5, 3.5, 4.5, 5.5, 6.5, Inf))) %>% 
    group_by(pinf.cat) %>% 
    summarise(n = n()) %>% 
    filter(n == max(n)) %>% 
    pull(pinf.cat)
  
  
  # Identify the log10 dose category containing the largest number of posterior ID70 draws
  levels_vec<-levels(cut(as.numeric(NA), breaks = c(0,2.5, 3.5, 4.5, 5.5, 6.5, Inf)))    # integer index saying which log-dose decade the posterior mode of log(ID70) falls into:(i.e. bin index representing which 10-fold dose range the posterior mode lies in)
  model.est.id70 <- match(as.character(id70.cat), levels_vec) + 1
  
  # If two or more dose categories are equally likely, choose the smaller dose (conservative/focus on safety)
  # handles the rare case of a tie (two bins have the same count). It picks the smaller index --> the lower dose, to be conservative (safety-oriented de-escalation bias)
  if(length(model.est.id70)>1){
    model.est.id70 <- min(model.est.id70)
  }
  
  
  
  # -----------------------------------------------------------------------------------
  # Summarise posterior infection probabilities at the candidate doses (using quantiles)
  # -----------------------------------------------------------------------------------
  est.prob.summary <- long.draws %>% 
    rowwise() %>%                                                                     #treat each row independently
    
    #Infection probabilities for each posterior draw
    mutate(est = list(enframe(dr.model(target.doses, as.vector(b),as.vector(a)), value = "p",  name = "dose"))) %>%   ## changed this line to "target.doses" to still get est probs for proposed schedule at each step
    unnest(cols = est) %>% #expand little tables into longer format
    group_by(dose) %>% 
    summarise(qs = list(enframe(quantile(p, c(0.05, 0.25, 0.5, 0.75, 0.95))))) %>%    #posterior infection-probability quantile summaries for each dose
    unnest(cols = qs) %>% 
    pivot_wider(names_from = "name", values_from = "value") %>% 
    mutate(sim = i, n.doses = j, design = dd, scenario = ss, rn = row.numb)
  
  
  # Save posterior infection-probability summaries
  # (appends the newly computed posterior infection-probability quantile summaries for this simulation step to a growing file)
  write_csv(x = est.prob.summary,
            file = paste0(save.output.to, "/est.probs.csv"),
            append = TRUE, col_names = FALSE)
  
  
  
  # -----------------------------------------------------------------------
  # Summarise the posterior distribution of log10(ID70)  (using quantiles)
  # -----------------------------------------------------------------------
  # Posterior quantiles of the estimated ID70 (dose achieving 70% infection).
  id70.summary <- long.draws %>% 
    as_tibble() %>% 
    summarise(id70.ests = list(enframe(quantile(log10_id70, c(0.05, 0.25, 0.50, 0.75, 0.95))))) %>% 
    unnest(cols = id70.ests) %>% 
    pivot_wider(names_from = "name", values_from = "value") %>% 
    mutate(sim = i, n.doses = j, design = dd, scenario = ss, rn = row.numb)
  
  # Save posterior log10(ID70) summaries
  write_csv(x = id70.summary,
            file = paste0(save.output.to, "/est.id70s.csv"),
            append = TRUE, col_names = FALSE)
  
  
  # -----------------------------------------------------------------------
  # Return the current-dose summary
  # -----------------------------------------------------------------------
  current.dose <- doses[j]
  summary.stats <- est.prob.summary %>% 
    filter(dose == current.dose) %>% 
    select(`5%`:`95%` )
  
  
  return(cbind("ds" = current.dose, summary.stats, "model.est.id70" = model.est.id70))
 
}

