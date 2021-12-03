# Helper functions for the E step

ps = function(s, left, right, d, cur.par.S, cur.par.X){
  ab = truncbounds.x(s, l = left, r = right, d = d)
  Fx.a = pnorm(ab[,1], mean = cur.par.X[1], sd = cur.par.X[2])
  Fx.b = pnorm(ab[,2], mean = cur.par.X[1], sd = cur.par.X[2])
  dnorm(s, mean = cur.par.S[1], sd = cur.par.S[2]) * (Fx.b - Fx.a)
}


Es_mc_imp = function(left, right, d, cur.par.S, cur.par.X, n = 1e3, max.n = 1e4){
  conv = F
  while(!conv){
  if(d ==2){
    xk = rnorm( n, mean = cur.par.S[1], sd = cur.par.S[2])
    q  = dnorm(xk, mean = cur.par.S[1], sd = cur.par.S[2])
   }
  if(d == 3){
    xk = rtruncnorm( n, a=-Inf, b=log(exp(right)-exp(left)), mean = cur.par.S[1], sd = cur.par.S[2])
    q  = dtruncnorm(xk, a=-Inf, b=log(exp(right)-exp(left)), mean = cur.par.S[1], sd = cur.par.S[2])
   }
  p  = ps(xk, left, right, d, cur.par.S, cur.par.X) 
  r  = p/q
  nc = mean( r )
  Ex = mean( xk * r ) /nc
  Ex2 = mean( xk^2 * r ) / nc
  if( !is.nan(Ex) | max.n < n ) conv = T else n = 10*n
  }
  c(Ex,Ex2)
}


