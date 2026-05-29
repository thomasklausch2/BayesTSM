# BayesTSM

`BayesTSM` fits Bayesian progressive three-state semi-Markov models for screening panel data with censoring after intervention. It is designed for settings where individuals are repeatedly screened for a progressive disease process, for example:

1. state 1: disease-free,
2. state 2: pre-state disease / early disease,
3. state 3: disease / late disease.

The package models the latent transition time from baseline to state 2, denoted `X`, and the transition time from state 2 to state 3, denoted `S`. The total time from baseline to state 3 is `Y = X + S`.

`BayesTSM` uses accelerated failure time models for both transition times,

\[
\log X_i = z_{Xi}^\top \beta_X + \sigma_X \epsilon_i,
\qquad
\log S_i = z_{Si}^\top \beta_S + \sigma_S \xi_i,
\]

and estimates the model using Bayesian data augmentation and MCMC. Weibull, lognormal, and log-logistic transition time distributions are currently available.

## Installation

Install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github(
  "thomasklausch2/BayesTSM",
  build_vignettes = TRUE
)
```

## Basic usage

```r
library(BayesTSM)

set.seed(1)

dat <- gendat(
  n = 1000,
  p = 2,
  sigma.X = 0.3,
  mu.X    = 2,
  beta.X  = c(0.5, 0.5),
  sigma.S = 0.5,
  mu.S    = 1,
  beta.S  = c(0.5, 0.5),
  dist.X  = "weibull",
  dist.S  = "weibull",
  v.min   = 1,
  v.max   = 5,
  Tmax    = 200,
  mean.rc = 10
)

head(dat)
table(dat$d)
```

The variables passed to `bayestsm()` are:

- `d`: event type, where `1` denotes right censoring, `2` a pre-state event, and `3` a terminal-state event;
- `L`: left bound of the final screening interval;
- `R`: right bound of the final screening interval, or `Inf` for right-censored observations;
- `Z.X`: covariate matrix for the `X` transition model;
- `Z.S`: covariate matrix for the `S` transition model.

```r
d <- dat$d
L <- dat$L
R <- dat$R
Z <- dat[, c("Z.1", "Z.2")]

mod <- bayestsm(
  d = d,
  L = L,
  R = R,
  Z.X = Z,
  Z.S = Z,
  mc = 1e4,
  warmup = 5e2,
  thinning = 10,
  chains = 4,
  update_till_convergence = TRUE,
  mc_update = 1e4,
  MH = FALSE,
  dist.X = "weibull",
  dist.S = "weibull",
  seed_chains = 1:4
)
```

By default, `BayesTSM` uses slice sampling for the model parameters. Metropolis sampling can be used instead by setting `MH = TRUE`.

## Posterior summaries

```r
summary(mod, warmup = 500)
plot(mod, warmup = 500)
```

The summary method reports posterior medians, 95% credible intervals, R-hat values, and effective sample sizes.

## Posterior predictive transition probabilities

Posterior predictive cumulative transition probabilities can be obtained with `ppCIF()`.

```r
ppCIF(
  mod,
  type = "quantiles",
  warmup = 500,
  quant = c(5, 10)
)
```

For plotting cumulative transition probability curves:

```r
pp_grid <- ppCIF(
  mod,
  type = "quantiles",
  warmup = 500
)

plot(pp_grid, xlim = c(0, 50))
```

Conditional predictions can be obtained by fixing covariate values through `fix_Z.X` and `fix_Z.S`.

```r
ppCIF(
  mod,
  type = "quantiles",
  warmup = 500,
  quant = c(5, 10),
  fix_Z.X = c(1, NA),
  fix_Z.S = c(1, NA)
)
```

Here, the first covariate is fixed at 1, while the second covariate is marginalized over its observed distribution.

## Model comparison

Different transition time distributions can be compared using information criteria.

```r
get_IC(mod, warmup = 500, cores = NULL)
```

`get_IC()` returns the deviance information criterion (DIC) and two versions of the widely applicable information criterion (WAIC-1 and WAIC-2).

## References

Akwiwu, E. U., Coupé, V. M. H., Berkhof, J., & Klausch, T. (2026). A comparison of methods for modeling multistate cancer progression using screening data with censoring after intervention. *Medical Decision Making*. https://doi.org/10.1177/0272989X261422681

Klausch, T., Akwiwu, E. U., van de Wiel, M. A., Coupé, V. M. H., & Berkhof, J. (2023). A Bayesian accelerated failure time model for interval censored three-state screening outcomes. *The Annals of Applied Statistics, 17*(2), 1285–1306. https://doi.org/10.1214/22-AOAS1669

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. *Bayesian Analysis, 16*(2), 667–718. https://doi.org/10.1214/20-BA1221

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research, 11*, 3571–3594.
