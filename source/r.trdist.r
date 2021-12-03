# Truncated distribution for any distribution implemented in pdist and qdist
# Needs known densities and quantile functions to find samples from truncated distribution using standard sampling by uniform-CDF inversion
r.trdist = function(par, a = 0, b=Inf, dist){ 
  n = length(a)
  cdf.a = pdist(a, par, dist)
  cdf.b = pdist(b, par, dist)
  dif = cdf.b-cdf.a
  dif = dif < 1e-08
  u   = ifelse(dif, a, qdist( cdf.a + runif(n) * (cdf.b-cdf.a), par, dist ))
  return(u)
}
