// Dose–response: p(d) = 1 - (1 + d/b)^(-a)
// Priors: log10(ID70) ~ TruncNormal(mean=4, sd=0.8) on [3,7]
//         log(a) ~ Uniform(log(0.3), log(0.8))
// b is derived from ID70 and a via p(ID70)=0.7

data {
  int<lower=1> N;
  array[N] real<lower=0> ds;           // doses
  array[N] int<lower=0>  nn;           // cohort sizes
  array[N] int<lower=0>  obs;          // infected counts
}

parameters {
  real<lower=log(0.3), upper=log(0.8)> log_a;
  
  // log10(ID70) truncated to [3,7]
  real<lower=3, upper=7> log10_id70;
}

transformed parameters {
  real a = exp(log_a);
  
  // ID70 on natural scale
  real ID70 = pow(10.0, log10_id70);

  // derive b so that p(ID70)=0.7
  real b = ID70 / (pow(0.3, -1.0 / a) - 1.0);
  
  
  array[N] real<lower=0, upper=1> p;
  for (n in 1:N)
    p[n] = 1 - pow(1 + ds[n] / b, -a);
}

model {
  
  // Truncated normal prior for log10_id70:
  log10_id70 ~ normal(4, 0.8);
  
  for (n in 1:N)
    obs[n] ~ binomial(nn[n], p[n]);
}
