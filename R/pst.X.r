#' Data augmentation for the latent `X` times
#'
#' Draws augmented values of the latent transition time `X` conditional on the
#' observed data, the latent `S` times, and current parameter values.
#'
#' @param par Numeric matrix of distribution parameters, with one row per
#'   observation.
#' @param dist Character string specifying the distribution used for `X`.
#' @param S Numeric vector of latent `S` times.
#' @param d Numeric vector of observed event categories.
#' @param L Numeric vector of left interval bounds.
#' @param R Numeric vector of right interval bounds.
#'
#' @return A numeric vector of augmented latent `X` values.
#'
#' @details
#' The sampling scheme depends on the observed event category `d`. For
#' right-censored observations (`d = 1`), `X` is sampled from its truncated
#' conditional distribution over the interval `[L, R]`. For observations with
#' `d = 2`, the truncation interval depends on both the observation interval and
#' the latent `S` value. For observations with `d = 3`, `X` is sampled from the
#' truncated interval implied by the requirement that both transitions occur
#' before the observed right endpoint.
#' @keywords internal
#' @noRd
# pst.X  = function(par, dist, S, d, L, R){
#   q = rep(NA, length(d))
#   q[d==1] = r.trdist(par[d==1,], a = L[d==1], b = R[d==1], dist)
#   q[d==2] = r.trdist(par[d==2,], a = apply ( cbind(L[d==2], R[d==2]-S[d==2]), 1, max), b=R[d==2], dist)
#   q[d==3] = r.trdist(par[d==3,], a = L[d==3] , b=R[d==3]-S[d==3], dist)
#   return(q)
# }
pst.X  = function(par, dist, S, d, L, R){
  q = rep(NA, length(d))
  q[d==1] = r.trdist(par[d==1,], a = L[d==1], b = R[d==1], dist)
  q[d==2] = r.trdist(par[d==2,], a = pmax(L[d==2], R[d==2]-S[d==2]), b=R[d==2], dist)
  q[d==3] = r.trdist(par[d==3,], a = L[d==3] , b=R[d==3]-S[d==3], dist)
  return(q)
}
