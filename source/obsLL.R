# Impelemntion of the observed likelihood obsLL for use in get.IC.r

f.c2 = function(x, par.X, par.S, dist.X, dist.S, R ){
  ddist(x, par.X, dist.X) * (1 - pdist(R-x, par.S, dist.S) )
}
f.c3 = function(x, par.X, par.S, dist.X, dist.S, R ){
  ddist(x, par.X, dist.X) * pdist(R-x, par.S, dist.S)
}

obsLL.c = function( dat, par.X, par.S , dist.X, dist.S, tol=1e-9 ){
  L = dat$L
  R = dat$R
  d = dat$d
  p1.X = ncol(par.X)
  p1.S = ncol(par.S)

  p = rep(NA, length(d))
  M = cbind(L,R, par.X, par.S)
  r1 = 1- pdist(L[d==1], par.X[d==1,], dist.X)
  r2 = apply( M[d==2,], 1, function(x){ integrate( f=f.c2 , lower=x[1], upper=x[2], 
                                                   par.X = t(as.matrix(x[3:4])), 
                                                   par.S = t(as.matrix(x[5:6])), 
                                                   dist.X, dist.S, R=x[2],  stop.on.error = F )$value })
  r3 = apply( M[d==3,], 1, function(x){ integrate( f=f.c3 , lower=x[1], upper=x[2], 
                                                   par.X = t(as.matrix(x[3:4])), 
                                                   par.S = t(as.matrix(x[5:6])), 
                                                   dist.X, dist.S, R=x[2],  stop.on.error = F )$value })
  p[d==1] = r1
  p[d==2] = r2
  p[d==3] = r3
  p <- (1-2*tol) * p + tol 
  ll=log(p)
  ll
}

obsLL = function(p, dat, Z.X, Z.S, dist.X, dist.S, ridge = F, lambda = NULL, sum.ll=T, log.scale=T){ 
  Z1.X = as.matrix(cbind(rep(1,nrow(dat)), Z.X))
  Z1.S = as.matrix(cbind(rep(1,nrow(dat)), Z.S))
  p1.X = ncol(Z1.X)
  p1.S = ncol(Z1.S)
  if( length(p) != p1.X+p1.S+2 ) stop('check p length')
  par.X = p[1:(p1.X+1)]
  par.S = p[(p1.X+2):(p1.X+p1.S+2)]
  if(dist.X!='lognormal'){
    tpar.X = trans.par(Z1.X, par.X)
  }
  if(dist.S!='lognormal'){
    tpar.S = trans.par(Z1.S, par.S)
  }
  if(dist.X=='lognormal'){
    m1 = Z1.X %*% as.matrix(par.X[1:p1.X])
    sd1 = exp(par.X[p1.X+1])
    tpar.X = cbind(m1, sd1)
  }
  if(dist.S=='lognormal'){
    m2 = Z1.S %*% as.matrix(par.S[1:p1.S])
    sd2 = exp(par.S[p1.S+1])
    tpar.S = cbind(m2, sd2)
  }
  ll = obsLL.c( dat= dat, par.X = tpar.X, par.S = tpar.S, 
               dist.X = dist.X, dist.S = dist.S)
  if(ridge) ll = ll - sum(p^2) * lambda /length(ll)
  if(sum.ll & log.scale) return( -( sum(ll) ) )
  if(!sum.ll & log.scale) return( ll )
  if(sum.ll & !log.scale) return( -( sum(exp(ll)) ) )
  if(!sum.ll & !log.scale) return( exp(ll) )
}

