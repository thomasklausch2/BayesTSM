ini.par.s <- function(time, event, Z, dist) {
  Z_df <- as.data.frame(Z)

  if (ncol(Z_df) == 0L) {
    formula.s <- survival::Surv(time = time, event = event) ~ 1
  } else {
    formula.s <- survival::Surv(time = time, event = event) ~ .
  }

  mod.s <- survival::survreg(
    formula = formula.s,
    data = Z_df,
    dist = dist
  )

  mu <- c(stats::coef(mod.s), log(mod.s$scale))
  Sigma <- stats::vcov(mod.s)

  par_ini <- MASS::mvrnorm(
    n = 1,
    mu = mu,
    Sigma = Sigma
  )

  par_ini <- as.numeric(par_ini)
  par_ini[length(par_ini)] <- exp(par_ini[length(par_ini)])

  par_ini
}
