# The logistic and log-logistic distributions, densities, and random number generators
dlog = function(x) exp(x)*(1+exp(x))^(-2)
qlog = function(p) log(p)-log(1-p)
rlog = function(n){
  u = runif(n)
  qlog(u)
}

dloglog = function(x,lambda, gamma){ dllogis(x, shape = gamma, scale= lambda ) }
ploglog = function(x,lambda, gamma){ pllogis(x, shape = gamma, scale= lambda ) }
qloglog = function(x,lambda, gamma){ qllogis(x, shape = gamma, scale= lambda ) }
rloglog = function(x,lambda, gamma){ rllogis(x, shape = gamma, scale= lambda ) }