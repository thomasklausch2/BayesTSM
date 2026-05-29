#' Generate initial parameter values for the X model from survreg estimates
#'
#' Fits a parametric survival regression model for the latent `X` transition
#' using interval-censored data and draws chain-specific initial values from a
#' multivariate normal approximation centered at the maximum likelihood
#' estimates.
#'
#' @param chains Integer scalar. Number of MCMC chains.
#' @param L Numeric vector of left interval bounds.
#' @param R Numeric vector of right interval bounds. May contain `Inf` for
#'   right-censored observations.
#' @param Z Data frame or matrix of covariates for the `X` model. If no
#'   covariates are used, pass a data frame with zero columns.
#' @param dist Character scalar. Distribution passed to
#'   [survival::survreg()]. Supported values are those accepted by
#'   `survreg()`, for example `"weibull"`, `"lognormal"`, and
#'   `"loglogistic"`.
#'
#' @return A numeric matrix with `chains` rows containing sampled initial
#'   parameter values. Columns correspond to the regression coefficients from
#'   `survreg()` followed by the scale parameter on the natural scale.
#'
#' @details
#' The function fits
#' \deqn{\texttt{Surv(L, R, type = "interval2")} \sim Z}
#' using [survival::survreg()].
#' Initial values are then sampled as
#' \deqn{\theta^{(c)} \sim N(\hat\theta, \widehat{\mathrm{Var}}(\hat\theta))}
#' for chains \eqn{c = 1, \dots, \texttt{chains}}, where the final component is
#' transformed back from `log(scale)` to `scale`.
#'
#' This initializer assumes that the parameterization used by `survreg()` is
#' compatible with the parameterization used in the MCMC sampler.
#'
#' @noRd
ini.par.x <- function(L, R, Z, dist) {
  Z_df <- as.data.frame(Z)

  L[L==0] <- -Inf
  if (ncol(Z_df) == 0L) {
    formula.x <- survival::Surv(time = L, time2 = R, type = "interval2") ~ 1
  } else {
    formula.x <- survival::Surv(time = L, time2 = R, type = "interval2") ~ .
  }

  dist_temp = dist
  if(dist == 'loglog') dist_temp = "loglogistic"
  mod.x <- survival::survreg(
    formula = formula.x,
    data = Z_df,
    dist = dist_temp
  )

  mu <- c(stats::coef(mod.x), log(mod.x$scale))
  Sigma <- stats::vcov(mod.x)

  par_ini <- MASS::mvrnorm(
    n = 1,
    mu = mu,
    Sigma = Sigma
  )

  par_ini <- as.numeric(par_ini)
  par_ini[length(par_ini)] <- exp(par_ini[length(par_ini)])

  par_ini
}
