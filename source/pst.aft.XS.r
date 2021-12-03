# Complete data posterior of AFT model, joint specification for X and S
pst.aft.XS = function(par, X=X, S=S, Z.X, Z.S, tau.X, sig.prior.X, tau.S, sig.prior.S, w=1, 
                      dist.X, dist.S, beta.prior = 't'){
  p1.X = ncol(Z.X)
  p1.S = ncol(Z.S)
  par.X = (par[1:(p1.X+1)])
  par.S = (par[(p1.X+2):(p1.X+p1.S+2)])
  
  pst.X = pst.aft(par=par.X, t=X, Z=Z.X, tau=tau.X, sig.prior=sig.prior.X, w=w, dist = dist.X, beta.prior = beta.prior)
  pst.S = pst.aft(par=par.S, t=S, Z=Z.S, tau=tau.S, sig.prior=sig.prior.S, w=w, dist = dist.S, beta.prior = beta.prior)
  pst.X + pst.S
}