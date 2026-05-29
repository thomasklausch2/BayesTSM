#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib BayesTSM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dexp dgamma dlnorm dnorm dt dweibull ecdf median optim pexp pgamma plnorm pnorm pweibull qexp qgamma qlnorm qnorm quantile qweibull rbeta rbinom rchisq rexp rgamma rlnorm rmultinom rnorm runif rweibull var complete.cases integrate time na.omit
#' @importFrom coda gelman.diag effectiveSize mcmc
#' @importFrom actuar dllogis pllogis qllogis rllogis
#' @importFrom MASS mvrnorm
#' @importFrom utils capture.output
#' @importFrom rlang .data
#' @importFrom graphics plot
#' @import   doParallel coda mvtnorm MCMCpack foreach parallel
## usethis namespace: end
NULL
