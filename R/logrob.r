#' Numerically robust logarithm
#'
#' Internal helper that computes a logarithm after replacing exact zeros by a
#' small tolerance value.
#'
#' @param x Numeric vector.
#' @param tol Positive numeric tolerance added to elements of `x` that are
#'   exactly zero.
#'
#' @return A numeric vector containing `log(x)`, with exact zeros replaced by
#'   `tol` before taking logs.
#'
#' @details
#' This function is used to avoid taking the logarithm of zero in situations
#' where small numerical underflow may occur.
#' @keywords internal
#' @noRd
logrob = function(x,tol){
  i = x==0
  x[i] = x[i]+tol
  log(x)
  }
