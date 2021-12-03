# Data augmentation for the X variable
pst.X  = function(par, dist, S, d, L, R){
  q = rep(NA, length(d))
  q[d==1] = r.trdist(par[d==1,], a = L[d==1], b = R[d==1], dist)
  q[d==2] = r.trdist(par[d==2,], a = apply ( cbind(L[d==2], R[d==2]-S[d==2]), 1, max), b=R[d==2], dist)
  q[d==3] = r.trdist(par[d==3,], a = L[d==3] , b=R[d==3]-S[d==3], dist)
  return(q)
}