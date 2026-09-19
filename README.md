# Bayesian CRM for dose-finding in controlled human infection model studies

This repository contains the R and Stan code used in the accompanying manuscript:

> Weerasinghe, C., Osowicki, J., Simpson, J. A., Crocker-Buque, T., McCarthy, J., Williams, E., and Price, D. J. (2026). *Rationale and guidance for implementing the continual reassessment method for dose-finding in controlled human infection model studies.*

The repository provides code to:

- generate true dose-response parameter values;
- implement the Bayesian continual reassessment method (CRM);
- simulate prior-predictive and fixed-scenario evaluations;
- implement the 5+5+10 rule-based comparator;
- reproduce the manuscript and supplementary tables and figures; and
- run an illustrative Bayesian CRM trial showing posterior updating after successive cohorts.


## Repository structure

```text
crm-chim-tutorial/
├── data/
├── data_gen/
├── bayesian_crm/
│   ├── R/
│   ├── stan/
│   └── results/
├── rule_based/
│   ├── R/
│   └── results/
├── fig_tables/
│   ├── R/
│   └── save_plots/
└── implementation_guide.qmd
```


## Software requirements

The analysis requires R, Stan and the following R packages:

```r
install.packages(
  c(
    "here",
    "tidyverse",
    "rstan",
    "truncnorm",
    "kableExtra",
    "patchwork"
  )
)
```

A working C++ toolchain is required to compile the Stan models.


## Implementation guide

Detailed instructions are provided in
[`implementation_guide.qmd`](implementation_guide.qmd).

The guide describes:

- parameter generation;
- Bayesian CRM implementation;
- the illustrative simulated trial;
- the rule-based comparator;
- reproduction of selected tables and figures; and
- reproduction of the complete simulation study.

To render the guide from the project root, use:

```bash
quarto render implementation_guide.qmd
```


## Running the simulations

Generate the true parameter values using:

```text
data_gen/generate_parameters.R
```

Run the Bayesian CRM simulations using:

```text
bayesian_crm/run_bcrm.R
```

Run the rule-based simulations using:

```text
rule_based/run_rule_based.R
```

The user-defined settings within these scripts determine:

- the prior specification;
- whether prior-predictive or fixed-scenario simulations are conducted;
- whether the capped or uncapped design is used; and
- the number of simulation replicates.


## Reproducing tables and figures

The manuscript and supplementary results can be reproduced using:

```text
fig_tables/tables.Rmd
fig_tables/figures.Rmd
```

The supporting summary, table and figure functions are provided in:

```text
fig_tables/R/
```


## Computational considerations

The full Bayesian CRM simulation study requires repeated Stan fitting and may take substantial computation time. Saved simulation outputs are included so that the manuscript tables and figures can be reproduced without rerunning the complete simulation study.


## Citation

If you use or adapt this code, please cite the accompanying manuscript:

> Weerasinghe, C., Osowicki, J., Simpson, J. A., Crocker-Buque, T., McCarthy, J., Williams, E., and Price, D. J. (2026). *Rationale and guidance for implementing the continual reassessment method for dose-finding in controlled human infection model studies.*
