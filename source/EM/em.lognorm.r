# EM algorithm for the lognormal-lognormal model
em.lognorm = function(d, L, R, Z.X, Z.S, 
                      sig.x = NULL, sig.s = NULL, beta.x = NULL, beta.s = NULL, tau = Inf,
                      n.rej = 1e3, tol = 1e-4 , max.it = 1e3, silent = F){
  dat = data.frame(d = d, L = L, R = R)
  n   = nrow(dat)
  Z1.X      = as.matrix(cbind( rep(1,n), Z.X))
  Z1.S      = as.matrix(cbind( rep(1,n), Z.S))
  
  if(is.null(sig.x)) sig.x = runif(1,exp(-2), exp(2))
  if(is.null(sig.s)) sig.s = runif(1,exp(-2), exp(2))
  if(is.null(beta.x)) beta.x = (as.matrix(runif(ncol(Z1.X),-1,1)))
  if(is.null(beta.s)) beta.s = (as.matrix(runif(ncol(Z1.S),-1,1)))
  par.X = cbind( Z1.X %*% (beta.x), sig.x)
  par.S = cbind( Z1.S %*% (beta.s), sig.s)
  if( is.infinite(tau) ){
    ll    =  -obsLL(p = c(beta.x,log(sig.x),beta.s,log(sig.s)) ,
                        dat = dat,
                        Z.X = Z1.X[,-1],
                        Z.S = Z1.S[,-1],
                        dist.X = 'lognormal',
                        dist.S = 'lognormal') 
  }
  if( !is.infinite(tau) ){
    ll    =  -obsLL(p = c(beta.x,log(sig.x),beta.s,log(sig.s)) ,
                        dat = dat,
                        Z.X = Z1.X[,-1],
                        Z.S = Z1.S[,-1],
                        dist.X = 'lognormal',
                        dist.S = 'lognormal')  +
                      (-(sig.x^2)^2/(2*tau^2) ) + 
                      (-(sig.s^2)^2/(2*tau^2) )
  }
  beta.x.seq = t(as.matrix(beta.x))
  beta.s.seq = t(as.matrix(beta.s))
  sig.x.seq = sig.x
  sig.s.seq = sig.s
  ll.seq = ll
  
  i = 1
  er = 1
  while(er > tol & i <= max.it){
    if(!silent){
    cat('Iteration', i-1,'\n')
    cat('logL: ', ll,'\n')
    cat('beta_x: ', beta.x,'\n')
    cat('beta_s: ', beta.s,'\n')
    cat('sigma_x: ',sig.x,'\n')
    cat('sigma_s: ',sig.s,'\n')
    }
    # E-step
    int.dat   = cbind(dat[,c('L','R','d')],par.X, par.S)
    int.dat$L = log(int.dat$L )
    int.dat$R = log(int.dat$R )
    int.dat1  = int.dat[int.dat$d == 1,]
    int.dat23 = int.dat[int.dat$d == 2 | int.dat$d == 3,]
    
    Ex = t( apply( cbind(int.dat23, 1),1, 
                   function(y) {Ex_mc_imp(left = y[1], right = y[2], d = y[3],
                                          cur.par.X = y[4:5], 
                                          cur.par.S = y[6:7],
                                          n = n.rej)} ) )

    Es = t( apply( cbind(int.dat23, 1),1, 
                   function(y) {Es_mc_imp(left = y[1], right = y[2], d = y[3],
                                          cur.par.X = y[4:5], 
                                          cur.par.S = y[6:7],
                                          n = n.rej)} ) )
    
    Ex_d1 = EX.truncnorm.rightcens( 
                mu  = int.dat1[,4], 
                sig = int.dat1[,5], 
                a   = int.dat1[,1]
              )
    
    Es_d1 = EX.norm(mu = int.dat1[,6], 
                    sig = int.dat1[,7])
    
    Ex.j = rep(NA, nrow(dat))
    Ex.j[dat$d == 1] = Ex_d1[,1]
    Ex.j[dat$d == 2 | dat$d == 3] = Ex[,1]
    
    Ex2.j = rep(NA, nrow(dat))
    Ex2.j[dat$d == 1] = Ex_d1[,2]
    Ex2.j[dat$d == 2 | dat$d == 3] = Ex[,2]
    
    Es.j = rep(NA, nrow(dat))
    Es.j[dat$d == 1] = Es_d1[,1]
    Es.j[dat$d == 2 | dat$d == 3] = Es[,1]
    
    Es2.j = rep(NA, nrow(dat))
    Es2.j[dat$d == 1] = Es_d1[,2]
    Es2.j[dat$d == 2 | dat$d == 3] = Es[,2]
    
    
    # M-step
    I = diag(ncol(Z1.X))
    
    D.X = t(Z1.X) %*% Z1.X  
    D.S = t(Z1.S) %*% Z1.S  
    beta.x = solve( D.X ) %*% t(Z1.X) %*% Ex.j
    beta.s = solve( D.S ) %*% t(Z1.S) %*% Es.j
    lt.x = Z1.X %*% beta.x
    lt.s = Z1.S %*% beta.s
    if( !is.infinite(tau) ){
      sig.x = pen.sig.em(E = Ex.j, E2 = Ex2.j, X=Z1.X, b=beta.x, tau = tau^2)
      sig.s = pen.sig.em(E = Es.j, E2 = Es2.j, X=Z1.S, b=beta.s, tau = tau^2)
      ll    =  -obsLL(p = c(beta.x,log(sig.x),beta.s,log(sig.s)) ,
                      dat = dat,
                      Z.X = Z1.X[,-1],
                      Z.S = Z1.S[,-1],
                      dist.X = 'lognormal',
                      dist.S = 'lognormal')  +
        (-(sig.x^2)^2/(2*tau^2) ) + 
        (-(sig.s^2)^2/(2*tau^2) )
      }
    if( is.infinite(tau) ){
      sig.x = sqrt( mean( Ex2.j - 2 * Ex.j  *lt.x + lt.x^2 ) )
      sig.s = sqrt( mean( Es2.j - 2 * Es.j  *lt.s + lt.s^2 ) )
      ll    =  -obsLL(p = c(beta.x,log(sig.x),beta.s,log(sig.s)) ,
                      dat = dat,
                      Z.X = Z1.X[,-1],
                      Z.S = Z1.S[,-1],
                      dist.X = 'lognormal',
                      dist.S = 'lognormal') 
      }
    
    par.X = cbind(lt.x, sig.x)
    par.S = cbind(lt.s, sig.s)
    ll.seq[i+1] = ll
    beta.x.seq = rbind(beta.x.seq, t(beta.x))
    beta.s.seq = rbind(beta.s.seq, t(beta.s))
    sig.x.seq[i+1] = sig.x
    sig.s.seq[i+1] = sig.s
    er = abs(ll.seq[i+1] - ll.seq[i])
    i = i + 1
  }
  if(!silent){
  cat('Iteration', i-1,'\n')
  cat('logL: ', ll,'\n')
  cat('beta_x: ', beta.x,'\n')
  cat('beta_s: ', beta.s,'\n')
  cat('sigma_x: ',sig.x,'\n')
  cat('sigma_s: ',sig.s,'\n')
  }
  ret = list()
  ret$beta.x = beta.x
  ret$beta.s = beta.s
  ret$sig.x  = sig.x
  ret$sig.s  = sig.s
  ret$beta.x.seq = beta.x.seq
  ret$beta.s.seq = beta.s.seq
  ret$sig.x.seq = sig.x.seq
  ret$sig.s.seq = sig.s.seq
  ret$ll = ll
  ret$ll.seq = ll.seq
  ret$convergence = er < tol
  ret$er = er
  ret$tol = tol
  ret
}

