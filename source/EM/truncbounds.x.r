# Helper functions for the E step
truncbounds.x = function(s, l, r, d){
  q = matrix(NA, nrow = length(s), ncol = 2)
  q[ d == 1, ] = c( l[d==1], Inf)
  q[ d == 2, ] = cbind( log( apply( cbind(exp(r[d==2])-exp(s), exp(l[d==2])), 1, function(y){ max(y[1],y[2])} )), r )
  q[ d == 3, ] = cbind( l, log( apply( cbind(exp(r[d==3])-exp(s), exp(l[d==3])), 1, function(y){ max(y[1],y[2])} ) ) )
  return(q)
}
