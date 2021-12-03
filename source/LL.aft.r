# Log likelihood AFT for different time distributions
LL.aft = function(par, t, Z , w, dist){
  p = length(par)
  beta=par[1:(p-1)]
  sigma=exp(par[p])
  log.sigma = par[p]
  n = length(t)
  V = (log(t)- Z %*% as.matrix(beta))/sigma
  
  if(dist=='weibull'){
    if(length(w)==1){ 
      if(w==1) LL= -n*log.sigma + sum( ( V - exp(V) ) ) 
      }
    else{ LL= -sum(w)*log.sigma + sum( ( V - exp(V) ) *w )}
  }
  
  if(dist=='loglog'){
    if(length(w)==1){ 
      if(w==1) {
        logexpV = log(  1 + exp(-V) )
        logexpV[is.infinite(logexpV)]= -V[is.infinite(logexpV)]
        LL = -n*log.sigma - sum( V ) - 2 * sum( logexpV ) 
      }
      }
    else{ LL = -sum(w)*log.sigma - sum( w*V ) - 2 * sum( w * log(  1 + exp(-V) ) ) }
  }
  
  if(dist=='lognormal'){
    if(length(w)==1){ 
      if(w==1) LL = -n*log.sigma - 1/2 * sum( V^2 ) 
    }
    else{ stop('LL for lognormal with weights not available yet.') }
  }
  LL
}