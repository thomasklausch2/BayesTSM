# Wrapper function to switch between distributions, quantile functions, and random draws
pdist = function (q, par, dist = "exp"){
  if(dist == "exp"){
    return(pexp(q, rate = par[,1]))
  }
  if(dist == "weibull2"){
    return(pweibull2(q, beta = par[,2], theta = par[,1]))
  }
  if(dist == "gamma"){
    return(pgamma(q, shape = par[,1], rate = par[,2]))
  }
  if(dist == "weibull"){
    return(pweibull(q, shape = par[,2], scale = par[,1]))
  }
  if(dist == "loglog"){
    return(ploglog(q, lambda = par[,1], gamma = par[,2]))
  }
  if(dist == "lognormal"){
    return(plnorm(q, meanlog = par[,1], sdlog = par[,2]))
  }
}

qdist = function (p, par, dist = "exp"){
  if(dist == "exp"){
    return(qexp(p, rate = par[,1]))
  }
  if(dist == "weibull2"){
    return(q.weibull2(p, beta = par[,2], theta = par[,1]))
  }
  if(dist == "gamma"){
    return(qgamma(p, shape = par[,1], rate = par[,2]))
  }
  if(dist == "weibull"){
    return(qweibull(p, shape = par[,2], scale = par[,1]))
  }
  if(dist == "loglog"){
    return(qloglog(p, lambda = par[,1], gamma = par[,2]))
  }
  if(dist == "lognormal"){
    return(qlnorm(p, meanlog = par[,1], sdlog = par[,2]))
  }
}

rdist = function (n, par, dist = "exp"){
  if(dist == "exp"){
    return(rexp(n, rate = par[,1]))
  }
  if(dist == "weibull2"){
    return(rweibull2(n, beta = par[,2], theta = par[,1]))
  }
  if(dist == "gamma"){
    return(rgamma(n, shape = par[,1], rate = par[,2]))
  }
  if(dist == "weibull"){
    return(rweibull(n, shape = par[,2], scale = par[,1]))
  }
  if(dist == "loglog"){
    return(rloglog(n, lambda = par[,1], gamma = par[,2]))
  }
  if(dist == "lognormal"){
    return(rlnorm(n, meanlog = par[,1], sdlog = par[,2]))
  }
}

ddist = function (x, par, dist = "exp"){
  if(dist == "exp"){
    return(dexp(x, rate = par[,1]))
  }
  if(dist == "weibull2"){
    return(dweibull2(x, beta = par[,2], theta = par[,1]))
  }
  if(dist == "gamma"){
    return(dgamma(x, shape = par[,1], rate = par[,2]))
  }
  if(dist == "weibull"){
    return(dweibull(x, shape = par[,2], scale = par[,1]))
  }
  if(dist == "loglog"){
    return(dloglog(x, lambda = par[,1], gamma = par[,2]))
  }
  if(dist == "lognormal"){
    return(dlnorm(x, meanlog = par[,1], sdlog = par[,2]))
  }
}
