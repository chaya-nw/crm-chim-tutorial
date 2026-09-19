# Function for running Bayesian CRM dose-finding simulations

run_bcrm <- function(
    cap_method,
    generation_method,
    prior_spec,
    target_prob,
    id70_settings,
    ab_settings,
    N_cap,
    n_random,
    n_sim_scenario,
    seed.no,
    designs,
    target.doses,
    compiled_model,
    chains,
    iter,
    warmup,
    save.output.to
) {
  
  # -----------------------------------------------------------------------
  # Check selected options
  # -----------------------------------------------------------------------
  cap_method <- match.arg(
    cap_method,
    choices = c("capped", "uncapped")
  )
  
  generation_method <- match.arg(
    generation_method,
    choices = c("random", "scenarios")
  )
  
  prior_spec <- match.arg(
    prior_spec,
    choices = c("id70_p", "ab_p")
  )
  
  # -----------------------------------------------------------------------
  # Read true parameter values and construct simulation index
  # -----------------------------------------------------------------------
  if (generation_method == "random") {
    
    nsims <- n_random
    
    parameter_file <- paste0(
      "pars_random_",
      prior_spec,
      ".csv"
    )
    
    pars_random <- readr::read_csv(
      here::here(
        "data",
        parameter_file
      ),
      show_col_types = FALSE
    )
    
    scenarios <- pars_random %>%
      select(a, b)
    
    if (nrow(scenarios) < nsims) {
      stop(
        "The random-parameter file contains fewer than ",
        nsims,
        " parameter sets."
      )
    }
    
    my.inputs <- tibble(ss = 1:nsims, dd = 1, i = 1:nsims)
    
    #Index tibble telling the outer loop what to run:
    #ss: scenario index (1…nsims),
    #dd: design index (only 1 design),
    #i: simulation replicate ID (here it equals ss).
    
  } else {
    
    nsims <- n_sim_scenario
    
    parameter_file <- paste0(
      "pars_scenarios_",
      prior_spec,
      ".csv"
    )
    
    scenarios <- readr::read_csv(
      here::here(
        "data",
        parameter_file
      ),
      show_col_types = FALSE
    ) %>%
      select(a, b)
    
    # Run nsims replicate trials for each scenario
    my.inputs <- expand.grid(
      ss = 1:nrow(scenarios),
      dd = 1,
      i = 1:nsims
    ) %>%
      as_tibble()
  
  }
  
  # Total number of simulation jobs
  nnn <- nrow(my.inputs)
  
  # -----------------------------------------------------------------------
  # Define output files
  # -----------------------------------------------------------------------
  
  sims_file <- file.path(
    save.output.to,
    "sims.csv"
  )
  
  summary_stats_file <- file.path(
    save.output.to,
    "summary_stats.csv"
  )
  
  # Prevent results from being appended to existing output files
  if (
    file.exists(sims_file) ||
    file.exists(summary_stats_file)
  ) {
    stop(
      "Simulation output files already exist in: ",
      save.output.to,
      ". Remove or rename the existing files before rerunning."
    )
  }
  
  # -----------------------------------------------------------------------
  # Run simulation jobs
  # -----------------------------------------------------------------------
  
  tic <- Sys.time()
  
  for (j in seq_len(nrow(my.inputs))) {
    
    print(paste0(j, "/", nnn))
    
    # Obtain scenario, design, and simulation indices for the current run
    ss <- my.inputs$ss[j]
    dd <- my.inputs$dd[j]
    i  <- my.inputs$i[j]
    
    # Initialise the simulated trial
    my.dat <- NULL                              # accumulating trial data as it unfolds (rows: dose, infected, total)
    my.summary.stats <- NULL                    # accumulating per-cohort summaries from run_crm() (posterior/ID estimates etc)
    new.dose <- designs$start.dose[dd]          # current dose to administer (starts at 10^4)
    n.at.target.dose <- 0                       # count of participants at the current recommended dose (run this at each cohort)
    k <- 1                                      # cohort counter
    
    # Seed settings
    if (generation_method == "random") {
      
      set.seed(seed.no + j)
      
    } else {
      
      rep_seed <- seed.no + 100 * ss + i        # a function of scenario (1,2,3,4,5) and replication (1,2,...,100)
      set.seed(rep_seed)
    }
    
    
    # Stopping criteria for the capped and uncapped designs
    #
    # Both designs continue until 20 participants have been treated at the
    # current target dose.
    #
    # For the capped design, the trial also stops if enrolling another cohort
    # would exceed the maximum total sample size (N_cap).
    #
    # For the uncapped design, no maximum total sample size is applied.
    
    continue_trial <- function() {
      
      if (cap_method == "capped") {
        
        n.at.target.dose < 20 &&
          k * designs$group.size[dd] <= N_cap
        
      } else {
        
        n.at.target.dose < 20
      }
    }
    
    
    # ---------------------------------------------------------------------
    # Sequential CRM trial
    # ---------------------------------------------------------------------
    while (continue_trial()) {
      
      # Simulate one cohort at current dose
      my.dat <- my.sim.fn(
        dat = my.dat,                                    # accumulating data (rows: dose, infected, total)
        scenario = scenarios[ss, ],                      # true parameter values (a, b)
        new.dose = new.dose,                             # current dose
        new.group.size = designs$group.size[dd]          # cohort size (fixed at 5 for all cohorts)
      )
      
      # Perform the Bayesian CRM update
      summary.stats <- update_crm_model(
        dat = my.dat,
        scenario.no = ss,
        des.no = dd,
        sim.no = i,
        target.doses = target.doses,
        save.output.to = save.output.to,
        row.numb = j,
        compiled_model = compiled_model,
        chains = chains,
        iter = iter,
        warmup = warmup,
        prior_spec = prior_spec,
        target_prob = target_prob,
        id70_settings = id70_settings,
        ab_settings = ab_settings
      )
      
      my.summary.stats <- rbind(
        my.summary.stats,
        summary.stats
      )
      
      
      # Apply the dose-allocation criteria
      # The decision function also restricts the proposed dose to the allowable
      # range from 10^3 to 10^7 CFU.
      dec <- decision.criteria(
        dat = my.dat,
        summary.stats = summary.stats
      )
      
      
      # Update the cumulative number of participants treated at the current
      # target dose; the trial stops once this number reaches 20
      n.at.target.dose <- as.numeric(
        dec["target.n"]
      )
      
      # Obtain the 'proposed next dose'
      prop.dose <- as.numeric(
        dec["new.dose"]
      )
      
      # Convert the 'current dose' from the original CFU scale to the log10 scale
      new.dose <- log10(new.dose)
      
      
      # If the proposed dose is above the current dose and the model-estimated
      # ID70 is above the proposed dose, move further towards the estimated ID70,
      # subject to a maximum increase of two log10 dose levels from the current dose.
      if (
        # i.e. model-estimated ID70 > proposed dose > current dose
        summary.stats$model.est.id70 > prop.dose &&
        prop.dose > new.dose
      ) {
        
        new.dose <- dplyr::case_when(
          
          # If the estimated ID70 is more than two log10 levels above the
          # current dose, escalate by the maximum permitted amount
          summary.stats$model.est.id70 - new.dose > 2 ~
            new.dose + 2,
          
          # Otherwise, move directly to the model-estimated ID70
          summary.stats$model.est.id70 - new.dose <= 2 ~
            summary.stats$model.est.id70
        )
        
      } else {
        
        # Otherwise, use the dose proposed by the decision criteria
        new.dose <- prop.dose
      }
      
      # Convert the selected next dose back to the original CFU scale
      new.dose <- 10^new.dose
      
      # Move to the next cohort
      k <- k + 1
    }
    
    
    # ---------------------------------------------------------------------
    # Save the completed simulated trial
    # ---------------------------------------------------------------------
    
    # Append the cohort-level trial data to the simulation-results file.
    # The added variables identify the true scenario, design, simulation
    # replicate, and simulation-job index.
    readr::write_csv(
      x = my.dat %>%
        dplyr::mutate(
          scenario = ss,  # index of the true parameter scenario: for "random", it indexes one of the randomly generated true parameter sets; for "scenarios", it indexes one of the five fixed true parameter scenarios.
          design = dd,    # index of the CRM design
          sim = i,        # simulation replicate number within the scenario
          kk = j          # overall simulation-job index (row number in my.inputs)
        ),
      file = sims_file,
      append = TRUE,
      col_names = FALSE
    )
    
    
    # Append the posterior summary statistics from all interim analyses in
    # the current trial to the summary-statistics file.
    readr::write_csv(
      x = my.summary.stats %>%
        dplyr::mutate(
          scenario = ss,
          design = dd,
          sim = i,
          kk = j
        ),
      file = summary_stats_file,
      append = TRUE,
      col_names = FALSE
    )
  }
  
  # -----------------------------------------------------------------------
  # Save computation time
  # -----------------------------------------------------------------------
  
  # Calculate the elapsed time since the simulation process started
  toc <- Sys.time() - tic
  
  # Convert the elapsed time to seconds
  elapsed_time <- as.numeric(
    toc,
    units = "secs"
  )
  
  # Define the path for the computation-time file
  time_file <- file.path(
    save.output.to,
    "time_exp.csv"
  )
  
  # Save the computation time in numeric and readable formats
  readr::write_csv(
    data.frame(
      elapsed_seconds = elapsed_time,
      pretty = paste(
        round(elapsed_time, 2),
        "secs"
      )
    ),
    file = time_file
  )
  
  
  # Return the output-file locations and elapsed time without automatically
  # printing them in the Console
  invisible(
    list(
      sims_file = sims_file,
      summary_stats_file = summary_stats_file,
      time_file = time_file,
      elapsed_time = elapsed_time
    )
  )
}