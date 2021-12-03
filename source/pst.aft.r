# Complete data posteriors for latent times variables
pst.aft = function(par, t, Z, tau, sig.prior, w, dist, beta.prior = 't'){
  p = length(par)
  beta = par[1:(p-1)]
  si   = exp(par[p])
  if(beta.prior == 't'){
  LL = LL.aft(par, t=t, Z=Z, w=w, dist=dist) + 
    log( dnorm(si,0,sig.prior) ) +
    sum( log( dt(beta, df = tau ) ) ) 
  }
  if(beta.prior == 'norm'){
  LL =  LL.aft(par, t=t, Z=Z, w=w, dist=dist) + 
      log( dnorm(si,0,sig.prior) ) +
      sum( log( dnorm(beta, 0, tau ) ) ) 
  }
  LL
}