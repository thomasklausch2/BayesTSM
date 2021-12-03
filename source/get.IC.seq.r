# This function computes information criteria DIC, WAIC-1 and WAIC-2 for a bts_survreg model
get.IC.seq = function(mod, samples = nrow(mod$par.X.bi), seed = NULL){
  set.seed(seed)
  m.X = trim.mcmc(mod$par.X.all, burnin = round(nrow(as.matrix(mod$par.X.all[1]))/2)+1)
  m.S = trim.mcmc(mod$par.S.all, burnin = round(nrow(as.matrix(mod$par.S.all[1]))/2)+1)
  m.X = as.matrix(m.X)
  m.S = as.matrix(m.S)
  if(samples >  nrow(m.X)) {stop ('more samples than mcmc draws selected')}
  
  ncol.X= ncol(m.X)
  ncol.S= ncol(m.S)
  m.X[,ncol.X] = log(m.X[,ncol.X])
  m.S[,ncol.S] = log(m.S[,ncol.S])
  m = cbind(m.X,m.S)
  m.s = m[sample(1:nrow(m), samples, replace=F),]
  
  pst.mean = apply(m, 2, mean) 
  run = apply(m.s,1, function(x) obsLL(x, dat=mod$dat, Z.X = mod$Z.X, Z.S = mod$Z.S, dist.X=mod$dist.X, 
                                            dist.S=mod$dist.S, log.scale = F, sum=F) ) 
  
  lppd = sum (log(apply(run,1, mean)))
  q1   = sum (apply( log(run), 1, mean))
  q2   = sum( apply( log(run), 1, var) )
  q3   = mean( apply( log(run),2, sum) )
  q4   = sum(obsLL(pst.mean, dat=mod$dat, Z.X = mod$Z.X, Z.S = mod$Z.S, dist.X=mod$dist.X, 
                   dist.S=mod$dist.S, log.scale = T, sum=F))
  DIC  = -2* ( q4 - 2*(q4-q3) )
  WAIC1    = -2*(-lppd + 2*q1)
  WAIC2    = -2*(lppd - q2)
  mat=matrix(nrow=1, ncol=3)
  mat[1,]=c(WAIC1, WAIC2, DIC)
  colnames(mat) = c('WAIC1', 'WAIC2', 'DIC')
  mat
}

