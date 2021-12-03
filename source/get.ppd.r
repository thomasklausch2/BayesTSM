# Calculates marginal predictive CIF from model mod fit by bts_survreg

unlist.par = function(par.list, Z, dist, pst.samples, s) {
  par = as.matrix(par.list[1])
  p = ncol(par)
  beta = par[,1:p-1, drop=F]
  sigma = par[,p]
  for(i in 1:length(par.list)){
    par = as.matrix(par.list[i])
    beta = rbind(beta,par[,1:(p-1), drop=F])
    sigma = c(sigma, par[,p])
  }
  linterm = cbind(1,Z) %*% t(beta[s,])
  if(dist != 'lognormal'){
    a = sigma[s]^-1
    b = exp(linterm)
  }
  if(dist == 'lognormal'){
    a = sigma[s] 
    b = linterm 
  }
  x = apply(rbind(a,b),2, function(x) rdist(n = length(x)-1, par = cbind(x[2:length(x)],x[1]), dist = dist ) )
  as.matrix(x)
}

get.ppd = function(mod, pst.samples=10^3, perc = seq(0, 1, 0.01)) {
  par.list.X = (mod$par.X.bi)
  par.list.S = (mod$par.S.bi)
  if( !is.null(mod$Z.X) ) Z.X = as.matrix(mod$Z.X) else Z.X = mod$Z.X # mod$Z.X[sample(1:nrow(mod$Z.X),nrow(mod$Z.X), replace=T),] 
  if( !is.null(mod$Z.X) ) Z.S = as.matrix(mod$Z.S) else Z.S = mod$Z.S # Z.X
  dist.X = mod$dist.X
  dist.S = mod$dist.S
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ret = list()
  x = unlist.par(par.list=par.list.X, Z.X, dist.X, pst.samples, s)
  s = unlist.par(par.list.S, Z.S, dist.S, pst.samples, s)
  q.x  = apply( x, 2, quantile, perc )
  q.s  = apply( s, 2, quantile, perc )
  q.xs = apply( x+s, 2, quantile, perc )
  ret$med.cdf.x     = apply(q.x, 1, median)
  ret$med.cdf.s     = apply(q.s, 1, median)
  ret$med.cdf.xs    = apply(q.xs, 1, median)
  ret$med.cdf.x.ci  = apply(q.x, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.s.ci  = apply(q.s, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.xs.ci = apply(q.xs, 1, quantile, c(0.025, 0.975))
  ret$perc = perc
  # ret$x.sample = as.numeric(ret$x)[sample(1:length(as.numeric(ret$x)), pst.samples, replace=F)]
  # ret$s.sample = as.numeric(ret$s)[sample(1:length(as.numeric(ret$s)), pst.samples, replace=F)]
  ret
}

get.ppd.grid = function(mod, pst.samples=10^3, q = seq(0, 20, 0.01)) {
  par.list.X = (mod$par.X.bi)
  par.list.S = (mod$par.S.bi)
  if( !is.null(mod$Z.X) ) Z.X = as.matrix(mod$Z.X) else Z.X = mod$Z.X # mod$Z.X[sample(1:nrow(mod$Z.X),nrow(mod$Z.X), replace=T),] 
  if( !is.null(mod$Z.X) ) Z.S = as.matrix(mod$Z.S) else Z.S = mod$Z.S # Z.X
  dist.X = mod$dist.X
  dist.S = mod$dist.S
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ret = list()
  x = unlist.par(par.list=par.list.X, Z.X, dist.X, pst.samples, s)
  s = unlist.par(par.list.S, Z.S, dist.S, pst.samples, s)
  p.x  = apply(x, 2, function(y){ ecdf.x = ecdf(y); ecdf.x(q) })
  p.s  = apply(s, 2, function(y){ ecdf.s = ecdf(y); ecdf.s(q) })
  p.xs = apply(x + s, 2, function(y){ ecdf.s = ecdf(y); ecdf.s(q) })
  ret$med.cdf.x     = apply(p.x, 1, median)
  ret$med.cdf.s     = apply(p.s, 1, median)
  ret$med.cdf.xs    = apply(p.xs, 1, median)
  ret$med.cdf.x.ci  = apply(p.x, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.s.ci  = apply(p.s, 1, quantile, c(0.025, 0.975))
  ret$med.cdf.xs.ci = apply(p.xs, 1, quantile, c(0.025, 0.975))
  ret$q = q
  ret
}

get.pCIF.q = get.ppd
get.pCIF.p = get.ppd.grid


# Function for use with bayes.2S only
get.ppd.2S = function(mod, pst.samples=10^3) {
  par.list.X = (mod$par.X.bi)
  Z.X = as.matrix(mod$Z.X)
  dist.X = mod$dist.X
  s = sample(1:(length(par.list.X)*nrow(as.matrix(par.list.X[1]))), pst.samples, replace=F)
  ppd = unlist.par(par.list.X, Z.X, dist.X, pst.samples, s)
  ret = list()
  perc = seq(0, 1, 0.01)
  q = apply( ppd, 2, quantile, perc )
  ret$med.cdf = apply(q, 1, median)
  ret$med.cdf.ci = apply(q, 1, quantile, c(0.025, 0.975))
  ret$perc = perc
  ret$x = ppd
  ret$x.sample = as.numeric(ppd)[sample(1:length(as.numeric(ppd)), pst.samples, replace=F)]
  ret
}