#' Compute Posterior Predictive CIF Building Blocks
#'
#' Internal helper used by \code{ppCIF()} to compute posterior predictive
#' cumulative incidence function (CIF) quantities for \eqn{X}, \eqn{S}, and
#' \eqn{Y=X+S} from posterior draws. Depending on \code{method} and
#' \code{type}, the function returns posterior predictive CDF values on a grid
#' of time points or posterior predictive quantiles on a grid of probabilities.
#'
#' @param beta.X Matrix of posterior draws of regression coefficients for the
#'   \eqn{X}-submodel.
#' @param beta.S Matrix of posterior draws of regression coefficients for the
#'   \eqn{S}-submodel.
#' @param sigma.X Numeric vector of posterior draws of scale parameters for the
#'   \eqn{X}-submodel.
#' @param sigma.S Numeric vector of posterior draws of scale parameters for the
#'   \eqn{S}-submodel.
#' @param Z1.X Design matrix for the \eqn{X}-submodel, including the intercept.
#' @param Z1.S Design matrix for the \eqn{S}-submodel, including the intercept.
#' @param type Character string indicating whether posterior predictive
#'   inference is returned as CDF values on a time grid
#'   (\code{"quantiles"}) or as quantiles on a probability grid
#'   (\code{"percentiles"}).
#' @param method Character string specifying whether posterior predictive
#'   inference is obtained by closed-form CDF evaluation where available
#'   (\code{"analytic"}) or by Monte Carlo simulation (\code{"simulation"}).
#' @param dist.X Character string naming the distribution used for \eqn{X}.
#' @param dist.S Character string naming the distribution used for \eqn{S}.
#' @param grid Numeric vector giving the evaluation grid, interpreted according
#'   to \code{type}.
#'
#' @return A list of posterior predictive quantities used by \code{ppCIF()}.
#' @noRd
get_pp_stats = function(beta.X, beta.S, sigma.X, sigma.S, Z1.X, Z1.S, type, method,
                        dist.X, dist.S, grid) {

  linterm.X = Z1.X %*% t(beta.X)
  linterm.S = Z1.S %*% t(beta.S)
  if(dist.X != 'lognormal'){
    a.X = sigma.X^-1
    b.X = exp(linterm.X)
  }
  if(dist.X == 'lognormal'){
    a.X = sigma.X
    b.X = linterm.X
  }
  if(dist.S != 'lognormal'){
    a.S = sigma.S^-1
    b.S = exp(linterm.S)
  }
  if(dist.S == 'lognormal'){
    a.S = sigma.S
    b.S = linterm.S
  }
  n = nrow(Z1.X)
  L.X = list()
  for(i in seq_along(sigma.X)){
    L.X[[i]] <- cbind(b.X[,i], a.X[i])
  }
  L.S = list()
  for(i in seq_along(sigma.S)){
    L.S[[i]] <- cbind(b.S[,i], a.S[i])
  }

  # Generate random variates for simulation based inference
  X.gen = sapply(L.X, function(x) rdist(n = n, par = x, dist = dist.X ) )
  S.gen = sapply(L.S, function(x) rdist(n = n, par = x, dist = dist.S ) )
  Y.gen = X.gen + S.gen

  # Do analytic method
  ret=list()
  if(method == 'analytic' && type == 'percentiles') stop('Inference via percentiles is not currently available with method = analytic.')
  if(method == 'analytic' && type == 'quantiles'){
    pp.p.X = pp.p.S = matrix(nrow=length(grid), ncol = length(sigma.X))
    for(i in seq_along(grid)){
      p.X = sapply( L.X, function(x) pdist(q = grid[i], par = x, dist = dist.X ) )
      p.S = sapply( L.S, function(x) pdist(q = grid[i], par = x, dist = dist.S ) )
      pp.p.X[i,] = apply(p.X , 2, mean)
      pp.p.S[i,] = apply(p.S , 2, mean)
    }
    ret$p.X = pp.p.X
    ret$p.S = pp.p.S
    ret$p.Y = apply(Y.gen, 2, function(x) ecdf(x)(grid))
  }

  # Do simulation method
  if(method == 'simulation' && type == 'quantiles'){
    ret$p.X = apply(X.gen, 2, function(x) ecdf(x)(grid))
    ret$p.S = apply(S.gen, 2, function(x) ecdf(x)(grid))
    ret$p.Y = apply(Y.gen, 2, function(x) ecdf(x)(grid))
  }

  if(method == 'simulation' && type == 'percentiles'){
    ret$q.X = apply(X.gen, 2, function(x) quantile(x, probs = grid))
    ret$q.S = apply(S.gen, 2, function(x) quantile(x, probs = grid))
    ret$q.Y = apply(Y.gen, 2, function(x) quantile(x, probs = grid))
  }
  ret
}
