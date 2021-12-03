# Data augmentation for the S variable
pst.S   = function(par, dist, X, d, L, R){
  n = length(X)
  q = rep(NA, n)
  q[d==1] = rdist(sum(d==1), par[d==1,], dist)
  q[d==2] = r.trdist(par[d==2,], a = R[d==2]-X[d==2], b= rep(Inf,sum(d==2)), dist)
  q[d==3] = r.trdist(par[d==3,], a =rep(0,sum(d==3)) , b= apply( cbind(R[d==3]-X[d==3], R[d==3]-L[d==3]), 1, min), dist) #min!
  return(q)
}