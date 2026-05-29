#' Extreme value distribution utilities
#'
#' Quantile evaluation and random number generation for the standard extreme
#' value distribution.
#'
#' @param p Numeric vector of probabilities.
#' @param n Integer; number of observations to generate.
#'
#' @return
#' `q.ev()` returns a numeric vector of quantiles.
#' `r.ev()` returns a numeric vector of random draws.
#'
#' @name extremevalue
#' @noRd
q.ev = function(p) log(-log(1-p))

#' @rdname extremevalue
#' @noRd
r.ev = function(n){
  u = runif(n)
  q.ev(u)
}
