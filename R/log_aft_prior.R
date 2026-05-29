#' Log-prior for accelerated failure time models
#'
#' Evaluates the log-prior density for the parameter vector of an accelerated
#' failure time (AFT) model. The parameter vector `eta` contains the AFT
#' regression parameters and the log-residual standard deviation. The intercept
#' is stored in the first position, followed by the slope coefficients, and the
#' final element is the log of the residual standard deviation.
#'
#' That is,
#'
#' \deqn{
#'   \eta = (\beta_0, \beta_1, \ldots, \beta_p, \log(\sigma)).
#' }
#'
#' Custom AFT prior functions can be passed through the `bayestsm` interface.
#' Such functions must take `eta`, `tau`, `sig.prior` and `beta.prior` as arguments and must
#' return a single numeric value: the log-prior density. Only `tau` and
#' `sig.prior` are passed from the `bayestsm` interface. If a custom prior uses
#' additional prior parameters, these currently have to be hardcoded inside the
#' custom prior function.
#'
#' The AFT priors for the `x` and `s` models currently have to use the same
#' prior function. However, the prior parameters can differ between the two
#' models, since `bayestsm` passes `tau.x` and `sig.prior.x` to the `x` model
#' prior, and `tau.s` and `sig.prior.s` to the `s` model prior.
#'
#' The default prior places either independent Student-\eqn{t} priors or
#' independent normal priors on the regression coefficients. A normal prior is
#' placed on \eqn{\sigma}, with the corresponding Jacobian adjustment included
#' because the model is parameterized in terms of \eqn{\log(\sigma)}.
#'
#' @param eta Numeric parameter vector. The first element is the intercept,
#'   followed by the slope coefficients, with `log(sigma)` in the final
#'   position.
#' @param tau Numeric prior parameter. For `beta.prior = "t"`, this is the
#'   degrees of freedom of the Student-\eqn{t} prior. For
#'   `beta.prior = "norm"`, this is the standard deviation of the normal prior.
#' @param sig.prior Numeric standard deviation of the normal prior on
#'   \eqn{\sigma}.
#' @param beta.prior Character string specifying the prior for the regression
#'   coefficients. Currently `"t"` and `"norm"` are supported.
#'
#' @return A single numeric value giving the log-prior density.
#'
#' @examples
#' ## Example custom prior:
#' ## Same prior as the default Student-t prior, except that the prior on the
#' ## intercept eta[1] is relaxed by a factor 5.
#' log_aft_prior_relaxed_intercept <- function(eta, tau = 4, sig.prior = 1, beta.prior = 't') {
#'   p        <- length(eta)
#'   beta     <- eta[1:(p - 1)]
#'   logsigma <- eta[p]
#'   sigma    <- exp(logsigma)
#'
#'   dt(beta[1] / 5, df = tau, log = TRUE) - log(5) +
#'     sum(dt(beta[-1], df = tau, log = TRUE)) +
#'     dnorm(sigma, sd = sig.prior, log = TRUE) +
#'     logsigma
#' }
#'
#' @keywords internal
#' @export
log_aft_prior <- function(eta, tau = 4, sig.prior = 1, beta.prior = "t") {
  p        <- length(eta)
  beta     <- eta[1:(p - 1)]
  logsigma <- eta[p]
  sigma    <- exp(logsigma)

  if (beta.prior == "t") {
    return(
      sum(dt(beta, df = tau, log = TRUE)) +
        dnorm(sigma, sd = sig.prior, log = TRUE) +
        logsigma
    )
  }

  if (beta.prior == "norm") {
    return(
      sum(dnorm(beta, sd = tau, log = TRUE)) +
        dnorm(sigma, sd = sig.prior, log = TRUE) +
        logsigma
    )
  }
}
