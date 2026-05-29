#' Find valid starting values for one MCMC chain
#'
#' Repeatedly generates candidate starting values for the latent variables and
#' transition-model parameters of the progressive three-state model until both
#' posterior kernels are finite or the maximum number of attempts is reached.
#'
#' Starting values for the `X` model are obtained from `ini.par.x()`, followed
#' by latent `X` draws from the truncated transition-time distribution. Given
#' these `X` values, provisional `S` values are constructed for observations
#' with `d = 2` or `d = 3`, after which starting values for the `S` model are
#' obtained from `ini.par.s()`. The resulting parameter vectors are converted to
#' the internal parameterization used by the sampler, where the scale parameter
#' is stored on the log scale.
#'
#' @param L Numeric vector of left interval bounds.
#' @param R Numeric vector of right interval bounds.
#' @param d Integer vector of outcome categories taking values in `1`, `2`, or `3`.
#' @param Z.X Covariate matrix for the `X` model, without intercept.
#' @param Z.S Covariate matrix for the `S` model, without intercept.
#' @param Z1.X Design matrix for the `X` model including intercept.
#' @param Z1.S Design matrix for the `S` model including intercept.
#' @param p.X Integer. Number of covariates in the `X` model, excluding intercept.
#' @param p.S Integer. Number of covariates in the `S` model, excluding intercept.
#' @param dist.X Character string specifying the distribution for `X`.
#' @param dist.S Character string specifying the distribution for `S`.
#' @param beta.prior.X Prior hyperparameter for the `X` regression coefficients.
#' @param beta.prior.S Prior hyperparameter for the `S` regression coefficients.
#' @param sig.prior.X Prior hyperparameter or fixed value for the `X` scale.
#' @param sig.prior.S Prior hyperparameter or fixed value for the `S` scale.
#' @param fix.sigma.X Logical; whether `sigma_X` is fixed.
#' @param fix.sigma.S Logical; whether `sigma_S` is fixed.
#' @param beta.prior Character string specifying the prior family.
#' @param max_tries Integer. Maximum number of attempts to find finite starting
#'   values. Default is `100`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`x.ini`}{Numeric vector of initial latent `X` values.}
#'   \item{`s.ini`}{Numeric vector of initial latent `S` values, with `NA` for `d = 1`.}
#'   \item{`beta.x.ini`}{Initial `X` parameter vector excluding log-scale.}
#'   \item{`sigma.x.ini`}{Initial log-scale parameter for the `X` model.}
#'   \item{`beta.s.ini`}{Initial `S` parameter vector excluding log-scale.}
#'   \item{`sigma.s.ini`}{Initial log-scale parameter for the `S` model.}
#'   \item{`par.x`}{Raw starting parameter vector returned by `ini.par.x()`.}
#'   \item{`par.s`}{Raw starting parameter vector returned by `ini.par.s()`.}
#'   \item{`log.pst.X`}{Posterior kernel at the selected `X` starting values.}
#'   \item{`log.pst.S`}{Posterior kernel at the selected `S` starting values.}
#'   \item{`counter`}{Number of failed attempts before success.}
#' }
#'
#' @noRd
find_starting_values <- function(L, R, d, Z.X, Z.S, Z1.X, Z1.S, p.X, p.S,
                                 dist.X, dist.S,
                                 beta.prior.X, beta.prior.S,
                                 sig.prior.X, sig.prior.S,
                                 fix.sigma.X, fix.sigma.S,
                                 beta.prior,
                                 max_tries = 100,
                                 log_prior_fun) {
  n <- length(d)
  d23 <- d == 2 | d == 3

  go <- FALSE
  counter <- 0L

  while (!go) {
    par.x <- ini.par.x(L = L, R = R, Z = Z.X, dist = dist.X)
    if(dist.X != 'lognormal') trans.par.x <- trans.par(Z1 = Z1.X, par = par.x)
    if(dist.X == 'lognormal') trans.par.x <- trans.par.ind.norm(Z1 = Z1.X,
                                                                p = par.x[1:(length(par.x)-1)],
                                                                v= par.x[length(par.x)])
    x.ini <- r.trdist(
      par = trans.par.x,
      a = L,
      b = R,
      dist = dist.X
    )

    s.ini <- rep(NA_real_, n)
    s.ini[d == 3] <- stats::runif(
      n = sum(d == 3),
      min = 0,
      max = R[d == 3] - x.ini[d == 3]
    )
    s.ini[d == 2] <- R[d == 2] - x.ini[d == 2]


    par.s <- ini.par.s(
      time = s.ini[d23],
      event = d[d23] == 3,
      Z = Z.S[d23, , drop = FALSE],
      dist = dist.S
    )

    beta.x.ini <- par.x[1:(p.X + 1)]
    sigma.x.ini <- log(par.x[p.X + 2])
    beta.s.ini <- par.s[1:(p.S + 1)]
    sigma.s.ini <- log(par.s[p.S + 2])

    if (fix.sigma.X) {
      sigma.x.ini <- log(sig.prior.X)
    }
    if (fix.sigma.S) {
      sigma.s.ini <- log(sig.prior.S)
    }

    log.pst.X <- log_aft_compdatpst(
      eta = c(beta.x.ini, sigma.x.ini),
      y = x.ini,
      Z = Z1.X,
      tau = beta.prior.X,
      sig.prior = sig.prior.X,
      dist = dist.X,
      beta.prior = beta.prior,
      log_prior_fun = log_prior_fun
    )

    log.pst.S <- log_aft_compdatpst(
      eta = c(beta.s.ini, sigma.s.ini),
      y = s.ini[d23],
      Z = Z1.S[d23, , drop = FALSE],
      tau = beta.prior.S,
      sig.prior = sig.prior.S,
      dist = dist.S,
      beta.prior = beta.prior,
      log_prior_fun = log_prior_fun
    )

    if (is.finite(log.pst.X) && is.finite(log.pst.S)) {
      go <- TRUE
    } else {
      counter <- counter + 1L
      if (!is.finite(log.pst.X)) {
        message("Non-finite posterior at starting values for X. Repeating search.")
      }
      if (!is.finite(log.pst.S)) {
        message("Non-finite posterior at starting values for S. Repeating search.")
      }
    }

    if (counter > max_tries) {
      stop("Could not find starting values.", call. = FALSE)
    }
  }

  list(
    x.ini = x.ini,
    s.ini = s.ini,
    beta.x.ini = beta.x.ini,
    sigma.x.ini = sigma.x.ini,
    beta.s.ini = beta.s.ini,
    sigma.s.ini = sigma.s.ini,
    par.x = par.x,
    par.s = par.s,
    log.pst.X = log.pst.X,
    log.pst.S = log.pst.S,
    counter = counter
  )
}
