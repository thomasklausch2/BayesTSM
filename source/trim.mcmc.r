# Helper function useful for trimming and thinning mcmc lists
trim.mcmc = function(obj, burnin = 1, end = nrow (as.matrix(obj[1])), thining = 1){
  mcmc.list(lapply( obj, function(x) mcmc(x[seq(burnin,end,thining),], start = burnin, end=end, thin = thining) ))
}
