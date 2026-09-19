# Helper functions for the Bayesian CRM dose-finding simulations



# Evaluate the probability of infection under the dose-response model
#
# Arguments:
# DD: administered dose, which may be a single value or a vector
# b:  dose-response scale parameter
# a:  dose-response shape parameter
#
# Returns the infection probability:
# p(DD) = 1 - (1 + DD / b)^(-a)
dr.model <- function(DD, b, a){
  pp <- 1 - (1 + DD/b)^-a  
  names(pp) <- DD
  return(pp)
}



# Simulate one cohort at the current dose
#
# The function:
# 1. calculates the true probability of infection at the current dose using the dose-response parameters a and b;
# 2. simulates the number of infected participants in the cohort from a binomial distribution; and
# 3. appends the new cohort results to the accumulated trial data.
#
# Arguments:
# dat:            accumulated trial data from previously simulated cohorts
# scenario:       true dose-response parameter values a and b
# new.dose:       dose administered to the new cohort
# new.group.size: number of participants in the new cohort
#
# Returns a tibble with:
# ds:  administered dose
# obs: number of infected participants
# nn:  total number of participants in the cohort
my.sim.fn <- function(dat, scenario, new.dose, new.group.size){
  
  p.inf <- dr.model(new.dose, scenario$b, scenario$a)                            # calculate the true probability of infection at the current dose
  new.dat <- rbinom(1, size = new.group.size, prob = p.inf)                      # Simulate the number of infected participants in the new cohort
  out <- rbind(dat, c(new.dose, new.dat, new.group.size))                        # Append the new cohort results to the accumulated trial data
  colnames(out) <- c("ds","obs","nn")                                            # Assign variable names to the accumulated trial data
  return(as_tibble(out))
}




# Apply the CRM decision rules for dose allocation
#
# The function uses posterior infection-probability quantiles at the most
# recently administered dose to determine whether the next cohort should
# receive a higher, lower, or unchanged dose.
#
# The proposed dose is restricted to the allowable range from 10^3 to
# 10^7 CFU.
#
# Arguments:
# dat:           accumulated cohort-level trial data
# summary.stats: posterior infection-probability quantiles at the current dose
#
# Returns:
# target.n: total number of participants treated at the current dose
# new.dose: proposed next dose on the log10 scale

decision.criteria <- function(dat, summary.stats){
  
  target.dose <- dat %>% dplyr::slice(nrow(dat)) %>% pull(ds)                    # Identify the most recently administered dose
  target.n <- dat %>% filter(ds == target.dose) %>% summarise(nn = sum(nn))      # Calculate the total number of participants treated at this dose
  
  #Reads posterior quantiles of infection probability at current dose, and determine the proposed dose for the next cohort
  new.dose <- case_when(
    #If even the 25th percentile < 0.6 --> likely under-dosing --> escalate +1 log10 (×10 dose)
    summary.stats$`25%` < 0.6 ~ log10(target.dose) + 1,
    
    #Else if even the 5th percentile > 0.8 --> very likely over-dosing --> de-escalate −1 log10 (÷10 dose)
    summary.stats$`5%` > 0.8 ~ log10(target.dose) - 1,
    
    #Else → stay at the same log10 dose.
    TRUE ~ log10(target.dose)
  )
  
  # Restrict the proposed dose to the allowable log10 dose range
  if (new.dose < 3){
    new.dose <- 3
  } else if (new.dose > 7){
    new.dose <- 7
  } 
  
  return(c("target.n" = as.numeric(target.n), "new.dose" = new.dose))
}




