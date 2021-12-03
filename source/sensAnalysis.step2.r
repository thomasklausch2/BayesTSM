# Implements sensitvity analysis by modeling a latent common cause of X ond S 
sensAnalysis.step2 = function(mod, b.UX = 1, b.US = 1, parallel = T, 
                       mc=1000, chains=1, burnin=round(mc/2), thining=1, prop.sd.X=NULL, prop.sd.S=prop.sd.X,
                       update.burnin = T, fix.b.XS = T, prev.run = NULL){
  
    if(is.null(prev.run)){
    dims = dim(as.matrix(mod$par.X.all[1]))
    start.val.X = matrix(ncol=dims[2], nrow= chains )
    start.val.S = matrix(ncol=dims[2], nrow= chains )
    for(i in 1:chains){
      start.val.X[i,] = as.matrix(mod$par.X.all[i])[dims[1],]
      start.val.S[i,] = as.matrix(mod$par.S.all[i])[dims[1],]
    }
    d = mod$dat$d
    L = mod$dat$L
    R = mod$dat$R
    Z.X = mod$Z.X
    Z.S = mod$Z.S
    if(is.null(prop.sd.X)){ prop.sd.X =  mod$prop.sd.X }
    if(is.null(prop.sd.S)){ prop.sd.S =  mod$prop.sd.S }
    beta.prior.X = mod$priors$beta.prior.X
    beta.prior.S = mod$priors$beta.prior.S
    sig.prior.X = mod$priors$sig.prior.X
    sig.prior.S = mod$priors$sig.prior.S
    w = mod$w
    X.prev = mod$X
    S.prev = mod$S
    dist.X = mod$dist.X
    dist.S = mod$dist.S
    fix.sigma.X = mod$fix.sigma.X
    fix.sigma.S = mod$fix.sigma.S
    }
   if(!is.null(prev.run)){
     b.UX = prev.run$b.UX
     b.US = prev.run$b.US
     chains = length(prev.run$par.X.all)
     dims = dim(as.matrix(prev.run$par.X.all[1]))
     start.val.X = matrix(ncol=dims[2], nrow= chains )
     start.val.S = matrix(ncol=dims[2], nrow= chains )
     for(i in 1:chains){
       start.val.X[i,] = as.matrix(prev.run$par.X.all[i])[dims[1],]
       start.val.S[i,] = as.matrix(prev.run$par.S.all[i])[dims[1],]
     }
     d = prev.run$dat$d
     L = prev.run$dat$L
     R = prev.run$dat$R
     Z.X = prev.run$Z.X
     Z.S = prev.run$Z.S
     if(is.null(prop.sd.X)){ prop.sd.X =  prev.run$prop.sd.X }
     if(is.null(prop.sd.S)){ prop.sd.S =  prev.run$prop.sd.S }
     beta.prior.X = prev.run$priors$beta.prior.X
     beta.prior.S = prev.run$priors$beta.prior.S
     sig.prior.X = prev.run$priors$sig.prior.X
     sig.prior.S = prev.run$priors$sig.prior.S
     thining = prev.run$thining
     print('Updating previous MCMC run.') #Overriding data, chains, thining, burn-in, and priors
     w = prev.run$w
     X.prev = prev.run$X
     S.prev = prev.run$S
     dist.X = prev.run$dist.X
     dist.S = prev.run$dist.S
     fix.sigma.X = prev.run$fix.sigma.X
     fix.sigma.S = prev.run$fix.sigma.S
     beta.prior = prev.run$beta.prior
   }
  
  if(is.null(prop.sd.X)  ){
    stop("Provide proposal SD")
  } 
    if(!is.null(Z.X)) {
      Z.X = as.matrix(Z.X)
      if( is.null(colnames(Z.X)) ) colnames(Z.X) = paste('ZX.',1:ncol(Z.X))
      p.X = ncol(Z.X)
      Z1.X = cbind(1,Z.X)
    }
    if(!is.null(Z.S)) {
      Z.S = as.matrix(Z.S)
      if( is.null(colnames(Z.S)) ) colnames(Z.S) = paste('ZS.',1:ncol(Z.S))
      p.S = ncol(Z.S)
      Z1.S = cbind(1,Z.S)
    }
    if(is.null(Z.X)){  
      p.X = 0
      Z1.X  = matrix(1, ncol=1, nrow=length(d))
    }
    if(is.null(Z.S)){  
      p.S = 0
      Z1.S  = matrix(1, ncol=1, nrow=length(d))
    }  
  p1.X = ncol(Z1.X)
  p1.S = ncol(Z1.S)
  
  if(is.null(prev.run)){
  U = rnorm( length(d) )
  Z1.XU = cbind(Z1.X, U)
  Z1.SU = cbind(Z1.S, U)
  } else{
    Z1.XU = Z1.X
    Z1.SU = Z1.S
  }
  
  s = sample(1:10^5,chains, replace=F)
  ## Start run 
  print(paste("Starting Gibbs sampler with", chains, "chains and", mc, "iterations."))
  `%swichdo%` <- ifelse(parallel, `%dopar%`, `%do%`)
  if(parallel){
    cl    = makePSOCKcluster(chains)
    clusterSetRNGStream(cl)
    registerDoParallel(cl)
  }
  run = foreach(j = 1:chains,
                .export = c('pst.X', 'r.trdist', 'pdist', 'qdist', 'pst.S', 'rdist',
                              'trans.par', 'mhstep.aft', 'pst.aft', 'LL.aft', 'ploglog', 'rloglog',
                              'qloglog', 'pllogis', 'rllogis', 'dllogis', 'qllogis', 
                              'trans.par.norm', 'cor2cov', 'logrob', 'trans.par.ind.norm',
                              'mhstep.aft.XS', 'pst.aft.XS','mhstep.aft.XS.sens','mhstep.aft.XS.sens2')
  ) %swichdo% {
    library(MCMCpack)
    library(mvtnorm)
    set.seed(s[j])
    ndraws = mc
    
    # Get starting values from previous run
    if(p.X > 0) beta.X.ini = start.val.X[j,2:(p.X+1)] else beta.X.ini = NULL
    if(p.S > 0) beta.S.ini = start.val.S[j,2:(p.S+1)] else beta.S.ini = NULL
    sigma.X.ini = log(start.val.X[j,(p.X+2)])
    sigma.S.ini = log(start.val.S[j,(p.S+2)])
    mu.X.ini = start.val.X[j,1]
    mu.S.ini = start.val.S[j,1]
    n.prev = length(d)
    X = X.prev[ ((n.prev*(j-1))+1): (n.prev*j) ]
    S = S.prev[ ((n.prev*(j-1))+1): (n.prev*j) ]
    
    if(is.null(prev.run)){
      p.S  = p.S + 1
      p1.S = p1.S + 1
      p.X  = p.X + 1
      p1.X = p1.X + 1
      cur.par.Xreg = matrix(ncol=p.X+2, nrow=ndraws+1)
      cur.par.Sreg = matrix(ncol=p.S+2, nrow=ndraws+1)
      cur.par.Xreg[1,] = c(mu.X.ini,beta.X.ini, b.UX, sigma.X.ini) 
      cur.par.Sreg[1,] = c(mu.S.ini,beta.S.ini, b.US, sigma.S.ini) 
    }
    if(!is.null(prev.run)){
      cur.par.Xreg = matrix(ncol=p.X+2, nrow=ndraws+1)
      cur.par.Sreg = matrix(ncol=p.S+2, nrow=ndraws+1)
      cur.par.Xreg[1,] = c(mu.X.ini,beta.X.ini, sigma.X.ini) 
      cur.par.Sreg[1,] = c(mu.S.ini,beta.S.ini, sigma.S.ini) 
    }
    
    if( length(prop.sd.X) == 1 ) prop.sd.X.mat = diag(prop.sd.X^2,p1.X+1)
    if( length(prop.sd.S) == 1 ) prop.sd.S.mat = diag(prop.sd.S^2,p1.S+1)
    if( length(prop.sd.X) == 1 ) prop.sd.X.mat = diag(prop.sd.X^2,p1.X+p1.S+2) 
    
    #logXstd = scale( log(X) ) #
    #Z1.S = cbind(Z1.S, logXstd)
    ac.X = ac.S = 1
    
    # Begin Gibbs      

    for(i in 1:ndraws){
    # Update  X parameters
      par = t(as.matrix( c(cur.par.Xreg[i,, drop=F], cur.par.Sreg[i,, drop=F]) ) )
      mh  = mhstep.aft.XS.sens2( x = par,
                           X = X, S = S, Z.X = Z1.XU, Z.S = Z1.SU, 
                           tau.X = beta.prior.X, sig.prior.X = sig.prior.X, 
                           tau.S = beta.prior.S, sig.prior.S = sig.prior.S,
                           w = w, dist.X = dist.X, dist.S = dist.S,
                           prop.var=prop.sd.X.mat, fix.sigma.S=fix.sigma.S, fix.sigma.X=fix.sigma.X,
                           fix.b.XS = fix.b.XS)
      cur.par.Xreg[i+1,] = mh$s[1:(p1.X+1)]
      cur.par.Sreg[i+1,] = mh$s[(p1.X+2):(p1.X+p1.S+2)]
    
      # Assess acceptance for MH steps
      ac.X[i+1] = cur.par.Xreg[i+1,1] != cur.par.Xreg[i,1]
      ac.S[i+1] = cur.par.Sreg[i+1,1] != cur.par.Sreg[i,1]
      
      # augment X, S
      if(dist.X != 'lognormal') cur.par.X = trans.par(Z1.XU, par = cur.par.Xreg[(i+1),]) 
      if(dist.S != 'lognormal') cur.par.S = trans.par(Z1.SU, par = cur.par.Sreg[(i+1),])
      if(dist.X == 'lognormal') cur.par.X = trans.par.ind.norm(Z1 = Z1.XU, p = cur.par.Xreg[(i+1),1:p1.X], v= cur.par.Xreg[(i+1),(p1.X+1)])
      if(dist.S == 'lognormal') cur.par.S = trans.par.ind.norm(Z1 = Z1.SU, p = cur.par.Sreg[(i+1),1:p1.S], v= cur.par.Sreg[(i+1),(p1.S+1)])
      X = pst.X(par=cur.par.X, S=S, d = d, L= L, R= R, dist = dist.X) 
      S = pst.S(par=cur.par.S, X=X, d = d, L= L, R= R, dist = dist.S) 
      X[X==0] = 10^-300
      S[S==0] = 10^-300
    }
    cur.par.Xreg = cur.par.Xreg[2:(ndraws+1),]
    cur.par.Sreg = cur.par.Sreg[2:(ndraws+1),]
    out=list()
    out$X = X
    out$S = S
    out$par.X = cur.par.Xreg
    out$par.S = cur.par.Sreg
    out$ac.X = ac.X
    out$ac.S = ac.S
    out
  }
  if(parallel) stopCluster(cl)
  
  # Unwrap chains wihout trimming
  mcmc.par.X = mcmc.par.S = list()
  i=1
  par.X = run[[i]]$par.X
  par.S = run[[i]]$par.S
  # if(fix.b.XS) {
    # par.S = par.S[,-(p1.S+1)]
    # par.X = par.X[,-(p1.X+1)]
    # colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
    # colnames(par.S) = c("Intercept", colnames(Z.S), "sigma")
  # } else {
  if(is.null(prev.run)){
    colnames(par.X) = c("Intercept", colnames(Z.X), 'beta.UX', "sigma")
    colnames(par.S) = c("Intercept", colnames(Z.S), 'beta.US', "sigma")
  }  else{
    colnames(par.X) = c("Intercept", colnames(Z.X)[1:length(colnames(Z.X))-1], 'beta.UX', "sigma")
    colnames(par.S) = c("Intercept", colnames(Z.S)[1:length(colnames(Z.S))-1], 'beta.US', "sigma")
  }
  # }
  par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
  par.S[,ncol(par.S)] = exp(par.S[,ncol(par.S)])
  
  mcmc.par.X[[i]] = mcmc(par.X)
  mcmc.par.S[[i]] = mcmc(par.S)
  ac.X = run[[i]]$ac.X
  ac.S = run[[i]]$ac.S
  X.draw = run[[i]]$X
  S.draw = run[[i]]$S
  
  if(length(run)>1){
    for(i in 2:length(run)){
      par.X = run[[i]]$par.X
      par.S = run[[i]]$par.S
      # if(fix.b.XS) {
      # par.S = par.S[,-(p1.S+1)]
      # par.X = par.X[,-(p1.X+1)]
      # colnames(par.X) = c("Intercept", colnames(Z.X), "sigma")
      # colnames(par.S) = c("Intercept", colnames(Z.S), "sigma")
      # } else {
      if(is.null(prev.run)){
        colnames(par.X) = c("Intercept", colnames(Z.X), 'beta.UX', "sigma")
        colnames(par.S) = c("Intercept", colnames(Z.S), 'beta.US', "sigma")
      }  else{
        colnames(par.X) = c("Intercept", colnames(Z.X)[1:length(colnames(Z.X))-1], 'beta.UX', "sigma")
        colnames(par.S) = c("Intercept", colnames(Z.S)[1:length(colnames(Z.S))-1], 'beta.US', "sigma")
      }
      # }
      par.X[,ncol(par.X)] = exp(par.X[,ncol(par.X)])
      par.S[,ncol(par.S)] = exp(par.S[,ncol(par.S)])
      mcmc.par.X[[i]] = mcmc(par.X)
      mcmc.par.S[[i]] = mcmc(par.S)
      ac.X = cbind(ac.X, run[[i]]$ac.X)
      ac.S = cbind(ac.S, run[[i]]$ac.S)
      X.draw = c(X.draw, run[[i]]$X)
      S.draw = c(S.draw, run[[i]]$S)
    }
    }
  
  par.X.all = mcmc.list(mcmc.par.X)
  par.S.all = mcmc.list(mcmc.par.S)
 
  if(!is.null(prev.run)){
    par.X.all = bind.mcmclists(prev.run$par.X.all, par.X.all)
    par.S.all = bind.mcmclists(prev.run$par.S.all, par.S.all)
    ac.X.cur = ac.X
    ac.S.cur = ac.S
    ac.X = rbind(prev.run$ac.X, ac.X)
    ac.S = rbind(prev.run$ac.X, ac.S)
    mc.prev  = nrow(as.matrix(prev.run$par.X.all[1]))
    if(update.burnin) burnin   = round((mc.prev + mc)/2)
  }
  nr       = nrow(as.matrix(par.X.all[1]))
  par.X.bi = mcmc.list(lapply( par.X.all, function(x) mcmc(x[seq(burnin,nr,thining),], start = burnin, end=nr, thin = thining) ))
  par.S.bi = mcmc.list(lapply( par.S.all, function(x) mcmc(x[seq(burnin,nr,thining),], start = burnin, end=nr, thin = thining) ))
  
  ## save built info
  dat = data.frame(d=d, L=L, R=R)
  priors = list()
  priors$beta.prior.X = beta.prior.X
  priors$beta.prior.S = beta.prior.S
  priors$sig.prior.X = sig.prior.X
  priors$sig.prior.S = sig.prior.S
  
  
  ret = list()
  ret$par.X.all = par.X.all
  ret$par.S.all = par.S.all
  ret$par.X.bi = par.X.bi
  ret$par.S.bi = par.S.bi
  ret$X = X.draw
  ret$S = S.draw
  ret$ac.X = ac.X
  ret$ac.S = ac.S
  # if(!is.null(mod)){
  #   ret$ac.X.cur = ac.X.cur
  #   ret$ac.S.cur = ac.S.cur}
  ret$dat = dat
  ret$priors = priors
  ret$thining = thining
  ret$Z.S = Z1.SU[,-1]
  ret$Z.X = Z1.XU[,-1]
  ret$w = w
  ret$prop.sd.X = prop.sd.X
  ret$prop.sd.S = prop.sd.S
  ret$dist.X = dist.X
  ret$dist.S = dist.S
  ret$fix.sigma.X = fix.sigma.X
  ret$fix.sigma.S = fix.sigma.S
  ret$burnin = burnin
  ret$b.UX = b.UX
  ret$b.US = b.US
  return(ret)
}
