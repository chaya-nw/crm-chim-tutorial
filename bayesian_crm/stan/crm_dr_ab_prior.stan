// Dose–response: p(d) = 1 - (1 + d/b)^(-a)
// Priors: log(a) ~ Uniform(log(a_lower), log(a_upper))
//         log(b) ~ Uniform(log(b_lower), log(b_upper))

data {
  int<lower=1> N;
  array[N] real<lower=0> ds;           // doses
  array[N] int<lower=0>  nn;           // cohort sizes
  array[N] int<lower=0>  obs;          // infected counts
  
  // User-specified prior settings
  real<lower=0> a_lower;
  real<lower=a_lower> a_upper;

  real<lower=0> b_lower;
  real<lower=b_lower> b_upper;
}

parameters {
  real<lower=log(a_lower), upper=log(a_upper)> log_a;
  real<lower=log(b_lower),  upper=log(b_upper)> log_b;
}

transformed parameters {
  real a = exp(log_a);
  real b = exp(log_b);
  array[N] real<lower=0, upper=1> p;
  for (n in 1:N)
    p[n] = 1 - pow(1 + ds[n] / b, -a);
}

model {
  for (n in 1:N)
    obs[n] ~ binomial(nn[n], p[n]);
}
