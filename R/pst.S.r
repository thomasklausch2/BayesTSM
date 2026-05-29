#' Data augmentation for the latent `S` times
#'
#' Draws augmented values of the latent transition time `S` conditional on the
#' observed data and current parameter values.
#'
#' @param par Numeric matrix of distribution parameters, with one row per
#'   observation.
#' @param dist Character string specifying the distribution used for `S`.
#' @param X Numeric vector of latent `X` times.
#' @param d Numeric vector of observed event categories.
#' @param L Numeric vector of left interval bounds.
#' @param R Numeric vector of right interval bounds.
#'
#' @return A numeric vector of augmented latent `S` values.
#'
#' @details
#' The sampling scheme depends on the observed event category `d`. For
#' right-censored observations (`d = 1`), `S` is sampled from its full
#' conditional distribution. For observations with first-event detection
#' (`d = 2`) or second-event detection (`d = 3`), `S` is sampled from the
#' corresponding truncated conditional distribution implied by the interval
#' constraints.
#' @keywords internal
#' @noRd
pst.S   = function(par, dist, X, d, L, R){
  n = length(X)
  q = rep(NA, n)
  q[d==1] = rdist(sum(d==1), par[d==1,], dist)
  q[d==2] = r.trdist(par[d==2,], a = R[d==2]-X[d==2], b= rep(Inf,sum(d==2)), dist)
  q[d==3] = r.trdist(par[d==3,], a =rep(0,sum(d==3)) , b= apply( cbind(R[d==3]-X[d==3], R[d==3]-L[d==3]), 1, min), dist) #min!
  return(q)
}
