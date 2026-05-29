#' Joint Metropolis-Hastings update for `X` and `S` AFT parameters
#'
#' @keywords internal
#' @noRd
mhstep_aft_XS = function(x, X, S, d23, Z.X, Z.S, tau.X, sig.prior.X,
                         tau.S, sig.prior.S, dist.X, dist.S,
                         prop.var=diag(1,length(x)), fix.sigma.S=F, fix.sigma.X=F, beta.prior = 't',
                         log_prior_fun){

  x_ = mvrnorm(1, x, prop.var)
  p1.X = ncol(Z.X)
  if(fix.sigma.S){ x_[length(x_)] = log(sig.prior.S) }
  if(fix.sigma.X){ x_[p1.X+1] = log(sig.prior.X) }

  p.Y = log_aft_compdatpst_XS(eta = x_, X = X, S = S,
                              d23 = d23, Z.X = Z.X, Z.S = Z.S,
                              tau.X = tau.X, sig.prior.X = sig.prior.X,
                              tau.S = tau.S, sig.prior.S = sig.prior.S,
                              dist.X = dist.X, dist.S = dist.S,
                              beta.prior = beta.prior,
                              log_prior_fun = log_prior_fun)

  p.X = log_aft_compdatpst_XS(eta=x, X = X, S = S,
                              d23 = d23, Z.X = Z.X, Z.S = Z.S,
                              tau.X = tau.X, sig.prior.X = sig.prior.X,
                              tau.S = tau.S, sig.prior.S = sig.prior.S,
                              dist.X = dist.X, dist.S = dist.S,
                              beta.prior = beta.prior,
                              log_prior_fun = log_prior_fun)
  r   = p.Y - p.X
  u  = runif(1, 0, 1)
  ret = list()
  if (is.nan(r))
    r <- -Inf
  if( log(u) <= r ){
    ret$s = x_
    ret$p.X = p.Y
    ret$ac = 1
  }
  else{
    ret$s = x
    ret$p.X = p.X
    ret$ac = 0
  }
  return(ret)
}

