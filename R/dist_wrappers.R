#' Distribution wrappers
#'
#' Unified wrappers for density, distribution, quantile, and random generation
#' functions across the supported parametric families.
#'
#' @param x Numeric vector of evaluation points.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer; number of observations to generate.
#' @param par Numeric matrix of distribution parameters. The first column
#'   contains the first distribution parameter and the second column, when
#'   required, contains the second distribution parameter.
#' @param dist Character string specifying the distribution. Supported values
#'   are `"exp"`, `"weibull2"`, `"gamma"`, `"weibull"`, `"loglog"`, and
#'   `"lognormal"`.
#'
#' @return
#' `pdist()` returns cumulative probabilities.
#' `qdist()` returns quantiles.
#' `rdist()` returns random draws.
#' `ddist()` returns density values.
#'
#' @details
#' These functions provide a common interface to distribution-specific density,
#' distribution, quantile, and random generation functions. The interpretation
#' of `par` depends on `dist`:
#' \describe{
#'   \item{`"exp"`}{`par[, 1]` is the rate.}
#'   \item{`"weibull2"`}{`par[, 1]` is `theta` and `par[, 2]` is `beta`.}
#'   \item{`"gamma"`}{`par[, 1]` is the shape and `par[, 2]` is the rate.}
#'   \item{`"weibull"`}{`par[, 1]` is the scale and `par[, 2]` is the shape.}
#'   \item{`"loglog"`}{`par[, 1]` is `lambda` and `par[, 2]` is `gamma`.}
#'   \item{`"lognormal"`}{`par[, 1]` is `meanlog` and `par[, 2]` is `sdlog`.}
#' }
#'
#' @name dist_wrappers
#' @noRd
pdist = function (q, par, dist = "exp"){
  if(dist == "exp"){
    return(pexp(q, rate = par[,1]))
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

#' @rdname dist_wrappers
#' @noRd
qdist = function (p, par, dist = "exp"){
  if(dist == "exp"){
    return(qexp(p, rate = par[,1]))
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

#' @rdname dist_wrappers
#' @noRd
rdist = function (n, par, dist = "exp"){
  if(dist == "exp"){
    return(rexp(n, rate = par[,1]))
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


#' @rdname dist_wrappers
#' @noRd
ddist = function (x, par, dist = "exp"){
  if(dist == "exp"){
    return(dexp(x, rate = par[,1]))
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
