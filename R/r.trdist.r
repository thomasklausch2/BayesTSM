#' Random generation from a truncated distribution
#'
#' Generates random draws from a truncated distribution for any distribution
#' implemented through the package wrappers `pdist()` and `qdist()`.
#'
#' @param par Numeric matrix of distribution parameters, with one row per
#'   observation.
#' @param a Numeric vector of lower truncation bounds.
#' @param b Numeric vector of upper truncation bounds. Defaults to `Inf`.
#' @param dist Character string specifying the distribution.
#'
#' @return A numeric vector of random draws from the truncated distribution.
#'
#' @details
#' Sampling is performed by inverse transform sampling on the truncated support.
#' For each observation, the function evaluates the distribution function at the
#' lower and upper truncation bounds, draws a uniform random value on the
#' corresponding truncated CDF interval, and maps it back through the quantile
#' function.
#'
#' If the truncated probability mass is numerically negligible, the function
#' returns the lower bound `a` for that observation.
#' @keywords internal
#' @noRd
r.trdist = function(par, a = 0, b=Inf, dist){
  n = length(a)
  cdf.a = pdist(a, par, dist)
  cdf.b = pdist(b, par, dist)
  dif = cdf.b-cdf.a
  dif = dif < 1e-08
  u   = ifelse(dif, a, qdist( cdf.a + runif(n) * (cdf.b-cdf.a), par, dist ))
  return(u)
}
