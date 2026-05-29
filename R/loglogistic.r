#' Logistic and log-logistic distribution utilities
#'
#' Internal helper functions for the standard logistic distribution and the
#' parameterized log-logistic distribution.
#'
#' @param x Numeric vector of evaluation points.
#' @param p Numeric vector of probabilities.
#' @param n Integer; number of observations to generate.
#' @param lambda Scale parameter of the log-logistic distribution.
#' @param gamma Shape parameter of the log-logistic distribution.
#'
#' @return
#' `dlog()` returns logistic density values.
#' `qlog()` returns logistic quantiles.
#' `rlog()` returns random draws from the logistic distribution.
#' `dloglog()` returns log-logistic density values.
#' `ploglog()` returns log-logistic cumulative probabilities.
#' `qloglog()` returns log-logistic quantiles.
#' `rloglog()` returns random draws from the log-logistic distribution.
#'
#' @details
#' These functions are internal convenience wrappers used throughout the package.
#' The functions `dloglog()`, `ploglog()`, `qloglog()`, and `rloglog()` wrap the
#' corresponding log-logistic functions from \pkg{actuar} using the argument
#' names `lambda` for scale and `gamma` for shape.
#'
#' @name loglogistic
#' @keywords internal
#' @noRd
dlog = function(x) exp(x) * (1 + exp(x))^(-2)

#' @rdname loglogistic
#' @keywords internal
#' @noRd
qlog = function(p) log(p) - log(1 - p)

#' @rdname loglogistic
#' @keywords internal
#' @noRd
rlog = function(n){
  u = runif(n)
  qlog(u)
}

#' @rdname loglogistic
#' @keywords internal
#' @noRd
dloglog = function(x, lambda, gamma){ dllogis(x, shape = gamma, scale = lambda) }

#' @rdname loglogistic
#' @keywords internal
#' @noRd
ploglog = function(x, lambda, gamma){ pllogis(x, shape = gamma, scale = lambda) }

#' @rdname loglogistic
#' @keywords internal
#' @noRd
qloglog = function(x, lambda, gamma){ qllogis(x, shape = gamma, scale = lambda) }

#' @rdname loglogistic
#' @keywords internal
#' @noRd
rloglog = function(x, lambda, gamma){ rllogis(x, shape = gamma, scale = lambda) }
