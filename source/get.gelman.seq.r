# Helper file for the simulation to extract a sequence of convergence statistics from a MCMC list
get.gelman.seq = function(mclist, start, end, steps){
  d.mclist = dim( mclist[[1]])
  j = seq(start,end,steps)
  pest = upb = neff = matrix(ncol=d.mclist[2], nrow=length(j))
  mpsrf = numeric()
  for(i in 1:length(j)){
    s = round(j[i]/2)
    l = trim.mcmc(mclist, burnin = s, end = j[i] ) # mcmc.list( lapply(mclist[s:j[i],],as.mcmc ))
    g = gelman.diag( l )
    if(d.mclist[2] > 1) mpsrf[i] = g$mpsrf else mpsrf[i] = NA
    pest[i,] = g$psrf[,1]
    upb[i,]  = g$psrf[,2]
    neff[i,] = effectiveSize( l )
  }
  ret=list()
  ret$s = j
  ret$mpsrf = mpsrf
  ret$pest = pest
  ret$upb = upb
  ret$neff = neff
  ret
}
