#' Transform a correlation matrix to a covariance matrix
#'
#' Converts a correlation matrix to a variance-covariance matrix using a vector
#' of standard deviations.
#'
#' @param R Correlation matrix.
#' @param S Numeric vector of standard deviations.
#'
#' @return A variance-covariance matrix with standard deviations given by `S`
#'   and correlation structure given by `R`.
#' @keywords internal
#' @noRd
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}
