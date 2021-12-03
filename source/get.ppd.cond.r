# Calculates conditional predictive CIF from model mod fit by bts_survreg
get.ppd.cond = function(mod, Z = NULL, pst.samples=10^3, grid = seq(0,20,0.01)) {
  par.list.X = (mod$par.X.bi)
  par.list.S = (mod$par.S.bi)
  if( !is.null(mod$Z.X) ) Z.X = as.matrix(mod$Z.X) else Z.X = mod$Z.X
  if( !is.null(mod$Z.X) ) Z.S = as.matrix(mod$Z.S) else Z.S = mod$Z.S
  dist.X = mod$dist.X
  dist.S = mod$dist.S
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  
  ret = list()
  Q = unlist.par.cond(par.list.X, par.list.S, Z.X = Z, Z.S = Z, dist.X, dist.S, pst.samples, s, grid = grid, mc.xs = 1e4)
  ret$med.cdf.x     = apply(Q$p.X, 1, median)
  ret$med.cdf.s     = apply(Q$p.S, 1, median)
  ret$med.cdf.xs    = apply(Q$p.XS, 1, median)
  ret$med.cdf.x.ci  = apply(Q$p.X, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.s.ci  = apply(Q$p.S, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.xs.ci = apply(Q$p.XS, 1, quantile, c(0.025, 0.975))
  ret$grid = grid
  ret
}

unlist.par.cond = function(par.list.X, par.list.S, Z.X, Z.S, dist.X, dist.S, pst.samples, s, grid, mc.xs = 1e5) {
  par.X = as.matrix(par.list.X[1])
  par.S = as.matrix(par.list.S[1])
  p.X = ncol(par.X)
  p.S = ncol(par.S)
  beta.X = par.X[,1:p.X-1, drop=F]
  beta.S = par.S[,1:p.S-1, drop=F]
  sigma.X = par.X[,p.X]
  sigma.S = par.S[,p.S]
  for(i in 1:length(par.list.X)){
    par.X = as.matrix(par.list.X[i])
    beta.X = rbind(beta.X,par.X[,1:(p.X-1), drop=F])
    sigma.X = c(sigma.X, par.X[,p.X])
  }
  for(i in 1:length(par.list.S)){
    par.S = as.matrix(par.list.S[i])
    beta.S = rbind(beta.S,par.S[,1:(p.S-1), drop=F])
    sigma.S = c(sigma.S, par.S[,p.S])
  }
  linterm.X = beta.X[s,] %*% as.matrix(c(1,Z.X)) 
  linterm.S = beta.S[s,] %*% as.matrix(c(1,Z.S)) 
  if(dist.X != 'lognormal'){
    a.X = sigma.X[s]^-1
    b.X = exp(linterm.X)
  }
  if(dist.X == 'lognormal'){
    a.X = sigma.X[s] 
    b.X = linterm.X 
  }
  if(dist.S != 'lognormal'){
    a.S = sigma.S[s]^-1
    b.S = exp(linterm.S)
  }
  if(dist.S == 'lognormal'){
    a.S = sigma.S[s] 
    b.S = linterm.S 
  }
  ret=list()
  ret$p.X = apply(rbind(a.X,c(b.X)),2, function(x) pdist(q = grid, par = cbind(x[2:length(x)],x[1]), dist = dist.X ) )
  ret$p.S = apply(rbind(a.S,c(b.S)),2, function(x) pdist(q = grid, par = cbind(x[2:length(x)],x[1]), dist = dist.S ) )
  # 
  x = apply(rbind(a.X,c(b.X)),2, function(x) rdist(n = mc.xs, par = cbind(x[2:length(x)],x[1]), dist = dist.X ) )
  s = apply(rbind(a.S,c(b.S)),2, function(x) rdist(n = mc.xs, par = cbind(x[2:length(x)],x[1]), dist = dist.S ) )

  ret$p.XS = apply(x+s, 2, function(y){ ecdf.x = ecdf(y); 
                                        ecdf.x(grid) } )
  ret
  }

