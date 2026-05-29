#' Log-likelihood for accelerated failure time models
#'
#' @keywords internal
#' @noRd
LL_aft = function(eta, y, Z, dist){
  p = length(eta)
  beta=eta[1:(p-1)]
  sigma=exp(eta[p])
  log.sigma = eta[p]
  n = length(y)
  V = (log(y)- Z %*% as.matrix(beta))/sigma

  if(dist=='weibull') LL = -n*log.sigma + sum( V - exp(V) )

  if(dist=='loglog'){
    logexpV = log(  1 + exp(-V) )
    logexpV[is.infinite(logexpV)]= -V[is.infinite(logexpV)]
    LL = -n*log.sigma - sum( V ) - 2 * sum( logexpV )
  }

  if(dist=='lognormal') LL = -n*log.sigma - 1/2 * sum( V^2 )

  LL
}
