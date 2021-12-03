# Helper functions for the E step
truncbounds.s = function(x, l, r, d){
  q = matrix(NA,ncol=2, nrow=length(x))
  q[d==1,] = cbind( rep( -Inf, sum(d==1)), rep(Inf, sum(d==1)))
  q[d==2,] = cbind( log(apply(as.matrix( (exp(r[d==2])-exp(x))), 1, function(y) max(0,y))), Inf)
  q[d==3,] = cbind( -Inf, log(apply( cbind(0, exp(r[d==3])-exp(x)), 1, function(y) max(y[1],y[2]))))
  q
}
