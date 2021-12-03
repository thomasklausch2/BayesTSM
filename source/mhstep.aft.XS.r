# Metropolis steps for updating afT parameters of X and S jointly
mhstep.aft.XS = function(x, X, S, Z.X, Z.S, tau.X, sig.prior.X, 
                         tau.S, sig.prior.S, w=1, dist.X, dist.S,
                         prop.var=diag(1,length(x)),fix.sigma.S=F,fix.sigma.X=F, beta.prior = 't'){
  x_ = mvrnorm(1, x, prop.var)
  p1.X = ncol(Z.X)
  if(fix.sigma.S){ x_[length(x_)] = log(sig.prior.S) }
  if(fix.sigma.X){ x_[p1.X+1] = log(sig.prior.X) }
  p.Y = pst.aft.XS(par=x_, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S, beta.prior = beta.prior)
  p.X = pst.aft.XS(par=x, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S, beta.prior = beta.prior)
  r   = p.Y - p.X
  u  = runif(1, 0, 1)
  ret = list()
  if (is.nan(r)) 
    r <- -Inf
  if( log(u) <= r ){ 
    ret$s = x_
    ret$p.X = p.Y
  }
  else{
    ret$s = x
    ret$p.X = p.X
  }
  return(ret)
}

# Metropolis steps for updating afT parameters jointly for conditional independence sensitivity analysis
mhstep.aft.XS.sens = function(x, X, S, Z.X, Z.S, tau.X, sig.prior.X, 
                             tau.S, sig.prior.S, w=1, dist.X, dist.S,
                             prop.var=diag(1,length(x)),fix.sigma.S=F,fix.sigma.X=F, fix.b.XS = T){
  x_ = mvrnorm(1, x, prop.var)
  p1.X = ncol(Z.X)
  p1.S = ncol(Z.S)
  if(fix.sigma.X){ x_[p1.X+1] = log(sig.prior.X) }
  if(fix.sigma.S){ x_[length(x_)] = log(sig.prior.S) }
  if(fix.b.XS){ x_[p1.X+p1.S+1] = x[p1.X+p1.S+1] }
  p.Y = pst.aft.XS(par=x_, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S)
  p.X = pst.aft.XS(par=x, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S)
  r   = p.Y - p.X
  u  = runif(1, 0, 1)
  ret = list()
  if (is.nan(r)) 
    r <- -Inf
  if( log(u) <= r ){ 
    ret$s = x_
    ret$p.X = p.Y
  }
  else{
    ret$s = x
    ret$p.X = p.X
  }
  return(ret)
}

mhstep.aft.XS.sens2 = function(x, X, S, Z.X, Z.S, tau.X, sig.prior.X, 
                              tau.S, sig.prior.S, w=1, dist.X, dist.S,
                              prop.var=diag(1,length(x)),fix.sigma.S=F,fix.sigma.X=F, fix.b.XS = T){
  x_ = mvrnorm(1, x, prop.var)
  p1.X = ncol(Z.X)
  p1.S = ncol(Z.S)
  if(fix.sigma.X){ x_[p1.X+1] = log(sig.prior.X) }
  if(fix.sigma.S){ x_[length(x_)] = log(sig.prior.S) }
  if(fix.b.XS){ 
    x_[p1.X]        = x[p1.X] 
    x_[p1.X+p1.S+1] = x[p1.X+p1.S+1] 
  }
  p.Y = pst.aft.XS(par=x_, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S)
  p.X = pst.aft.XS(par=x, X=X, S=S, Z.X=Z.X, Z.S=Z.S, tau.X=tau.X, sig.prior.X=sig.prior.X, 
                   tau.S = tau.S, sig.prior.S = sig.prior.S, w=w, dist.X=dist.X, dist.S=dist.S)
  r   = p.Y - p.X
  u  = runif(1, 0, 1)
  ret = list()
  if (is.nan(r)) 
    r <- -Inf
  if( log(u) <= r ){ 
    ret$s = x_
    ret$logp_old = p.X
    ret$logp_new = p.Y
    ret$r = r
    ret$logu = log(u)
    ret$acc = log(u) <= r
  }
  else{
    ret$s = x
    ret$logp_old = p.X
    ret$logp_new = p.Y
    ret$r = r
    ret$logu = log(u)
    ret$acc = log(u) <= r
    }
  return(ret)
}