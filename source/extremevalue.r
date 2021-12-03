# Quantile and random number generation for the extreme value distribution
q.ev = function(p) log(-log(1-p))
r.ev = function(n){
  u = runif(n)
  q.ev(u)
}