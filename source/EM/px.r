# Helper functions for the E step

px = function(x, left, right, d, cur.par.S, cur.par.X){
  cd = truncbounds.s(x, l = left, r = right, d = d)
  Fs.c = pnorm(cd[,1], mean = cur.par.S[1], sd = cur.par.S[2])
  Fs.d = pnorm(cd[,2], mean = cur.par.S[1], sd = cur.par.S[2])
  ind_x = as.numeric( x > left & x <= right )
  ind_x * dnorm(x, mean = cur.par.X[1], sd = cur.par.X[2]) * (Fs.d - Fs.c)
}

library(truncnorm)
Ex_mc_imp = function(left, right, d, cur.par.S, cur.par.X, n = 1e3, max.n = 1e4){
  conv = F
  while(!conv){
  xk = rtruncnorm( n, a=left, b = right, mean = cur.par.X[1], sd = cur.par.X[2])
  q  = dtruncnorm(xk, a=left, b = right, mean = cur.par.X[1], sd = cur.par.X[2])
  p  = px(xk, left, right, d, cur.par.S, cur.par.X )
  r  = p / q
  nc = mean( r )
  Ex = mean( xk * r ) / nc
  Ex2 = mean( xk^2 * r ) / nc
  if( !is.nan(Ex) | max.n < n ) conv = T else n = 10*n
  }
  c(Ex,Ex2)
}

ind = function(x, s, l, r, d){
  if(d == 1) return( as.numeric( log(l) < x ) )
  if(d == 2) return( as.numeric( log(l) < x & x < log(r) & log(r) < log(exp(x) + exp(s)) )) 
  if(d == 3) return( as.numeric( log(l) < x & x < log(r) & log(exp(x) + exp(s)) < log(r) ))
}

px2.int = function(s, x, l, r, d, cur.par.S, cur.par.X){
  dnorm(x, mean = cur.par.X[1], sd = cur.par.X[2]) *
  dnorm(s, mean = cur.par.S[1], sd = cur.par.S[2]) * 
  ind(x, s , l, r, d)
}

px2 = function(x, l, r, d, cur.par.S, cur.par.X){
  apply(as.matrix(x), 1, function(y){
integral(px2.int, xmin = -Inf, xmax = Inf,
         method = c("Kronrod"),
         no_intervals = 8, random = FALSE,
         reltol = 1e-8, abstol = 0,
         l = l, r = r, d = d,
         cur.par.X = cur.par.X,
         cur.par.S = cur.par.S,
         x = y)})
}
