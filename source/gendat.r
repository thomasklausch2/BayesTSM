# Simulates interval censored 3 state screening data assuming MV-normal covariates
# See demo.r for details
gendat = function(n, p = 3, p.discrete = 0, r=.1, s = 1, 
                  sigma.X=1/2, mu.X = 1, beta.X = NULL, 
                  sigma.S=1/2, mu.S = 1, beta.S = NULL, 
                  cor.X.S = 0,
                  Tmax = 20, v.min = 1, v.max = 6, mean.rc = 40,
                  dist.X='weibull', dist.S='weibull',
                  do.b.XS = F, b.XS = NULL){
  
  if(p ==0){
    beta.X = as.matrix(mu.X)
    beta.S = as.matrix(mu.S)
  }   else{
  beta.X = as.matrix(c(mu.X, beta.X))
  beta.S = as.matrix(c(mu.S, beta.S))
  }
  
  # Sim Z
  R = matrix(r, p, p)
  diag(R) =  1
  S = rep(s, p)
  Sigma = cor2cov(R, S)
  if(p>0) {Z = mvrnorm(n, mu = rep(0,p), Sigma)} else Z = NULL
  # Sim discrete Z
  if(p.discrete == 1){
    Z.discrete = rbinom(n, 1, 0.5)
    Z = cbind( Z, Z.discrete)
    colnames(Z) = paste(1:ncol(Z))
  }  
  Z1 = cbind( as.matrix(rep(1,n)),Z)
  
  # Sim.X
  if( dist.X=='weibull' | dist.X=='loglog' | dist.X=='lognormal' ){
    if(dist.X=='weibull')    e.X = r.ev(n)
    if(dist.X=='loglog')     e.X = rlog(n)
    if(dist.X=='lognormal')  e.X = rnorm(n)
    X = Z1 %*% beta.X + sigma.X * e.X # log surv times
    X = exp(X) # surv times
    X  = as.numeric(X)
  }
  # Sim.S
  if( dist.S=='weibull' | dist.S=='loglog' | dist.S=='lognormal' ){
    if(dist.S=='weibull')   e.S = r.ev(n)
    if(dist.S=='loglog')    e.S = rlog(n)
    if(dist.S=='lognormal') e.S = rnorm(n)
    S = Z1 %*% beta.S + sigma.S * e.S # log surv times
    if(do.b.XS){
      #logstdX = scale(X)
      S = Z1 %*% beta.S + b.XS * log(X) + sigma.S * e.S # log surv times
    } 
    S = exp(S) # surv times
    S = as.numeric(S)
  }
  
  if(dist.X=='bv-lognormal'){
    Sigma.e  = matrix( c(1, cor.X.S, cor.X.S, 1), 2, 2)
    e = mvrnorm(n, c(0,0), Sigma.e)
    X = Z1 %*% rbind(mu.X, beta.X) + sigma.X * e[,1] # log surv times
    S = Z1 %*% rbind(mu.S, beta.S) + sigma.S * e[,2] # log surv times
    X = as.numeric(exp(X))
    S = as.numeric(exp(S))
  }
  
  # Create total time
  XS = as.numeric(X+S)  # Total time from baseline to event 2
  
  # Generate screening sequences
  t.rc = rexp(n, 1/mean.rc ) #runif(n, start.rc, vmax*(visitdist+leniancy))
  V = as.matrix(runif(n, v.min, v.max ))
  for(i in 2:Tmax){
    V = cbind(V, runif(n, V[,(i-1)]+v.min, V[,(i-1)]+v.max ))
  }
  t.rc = V[,1] + t.rc
  i = apply(V < t.rc, 1, sum)
  t.rc = apply( cbind(i, V), 1, function(x) x[x[1]+1] )
  t.rc = matrix( rep( t.rc, Tmax), nrow = n, ncol = Tmax)
  V    = (V < t.rc) * V + (V >= t.rc) * t.rc
  V    = cbind(0, V)
  
  # Find indices of events
  ind.h  = rowSums(V < X) # Last healthy time index
  ind.x  = ind.h+1        # First stage 1 
  ind.xs = rowSums(V < XS) + 1 # First stage 2
  ind.x[ind.x>ncol(V)] = ncol(V)
  
  # Define delta (event)
  d = rep(-99,n)
  d[ind.x < ind.xs ]  = 2
  d[ind.x == ind.xs  & ind.xs != (ncol(V)+1)] = 3
  d[ind.h == ind.x]  = 1
  
  # Define left and right interval bound
  L = V[ cbind(1:nrow(V), ind.h) ]
  R = V[ cbind(1:nrow(V), ind.x) ]
  R[ ind.h == ind.x ] = Inf
  
  if(p>0) dat = data.frame(L=L, R=R, d=d, X=X, S=S, Z=Z)
  if(p==0) dat = data.frame(L=L, R=R, d=d, X=X, S=S)
  
  return(dat)
}