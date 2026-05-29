#' Compute information criteria for a fitted `bayestsm` model
#'
#' Computes DIC, WAIC1, and WAIC2 for a fitted `bayestsm` model from posterior
#' draws of the model parameters.
#'
#' @param mod A fitted model object returned by [bayestsm()].
#' @param samples Integer; number of posterior draws used to approximate the
#'   information criteria. Defaults to the number of post-burn-in draws in the
#'   first chain of `mod$par.X` after removing the first `warmup` iterations per chain.
#'   Has to be smaller than `nrow(mod$par.X[[1]])` minus `warmup`.
#' @param warmup The number of iterations discarded for warmup per chain.
#'   The default is `NULL` in which case, conservatively,
#'   the first half of posterior samples is discarded.
#' @param cores Optional integer; number of cores used for parallel evaluation
#'   of the observation-level likelihood contributions. Defaults to 2. To set to maximum available cores set to `NULL`.
#'
#' @return A numeric matrix with one row and three columns:
#' \describe{
#'   \item{`WAIC1`}{The first WAIC estimate based on the mean log pointwise
#'   posterior density and a bias correction using the average log likelihood.}
#'   \item{`WAIC2`}{The second WAIC estimate based on the mean log pointwise
#'   posterior density and a bias correction using the posterior variance of the
#'   log likelihood contributions.}
#'   \item{`DIC`}{The deviance information criterion.}
#' }
#'
#' @details
#' The function by removes the first `warmup` iteration from the MCMC
#' output stored in `mod$par.X` and `mod$par.S`, then samples `samples`
#' posterior draws from the remaining iterations. If `warmup` is not specified, the first half of draws is removed, which is a conservative choice.
#' For each sampled draw, the observation-level likelihood contributions are computed, and
#' these are combined into the reported information criteria.
#'
#' Lower values of DIC, WAIC1, and WAIC2 indicate better expected predictive fit,
#' but these criteria should primarily be compared between models fitted to the
#' same data.
#'
#' @seealso [bayestsm()]
#'
#' @examples
#' \dontrun{
#' # For a bayestsm model output saved in mod_slice (see ?bayestsm examples), run:
#' get_IC(mod_slice, warmup = 500)
#' }
#'
#' @export
get_IC = function(mod, samples = nrow(mod$par.X[[1]])-1, warmup = NULL, cores = 2){
  if (!inherits(mod, "bayestsm")) {
    stop("`mod` must be of class 'bayestsm'.", call. = FALSE)
  }
  if(is.null(cores)) cores = detectCores()

  # Parameter matrix pre-processing
  it = nrow(mod$par.X[[1]])
  if(is.null(warmup)) {
    warmup = round(0.5 * it)
    message('Warmup not specified. Discarding the first half of posterior samples per chain.')
  }
  m.X <- trim.mcmc(mod$par.X, burnin = warmup)
  m.S <- trim.mcmc(mod$par.S, burnin = warmup)
  m.X = as.matrix(m.X)
  m.S = as.matrix(m.S)
  if(samples >  nrow(m.X)) {stop ('more samples than mcmc draws selected')}

  ncol.X= ncol(m.X)
  ncol.S= ncol(m.S)
  m.X[,ncol.X] = log(m.X[,ncol.X])
  m.S[,ncol.S] = log(m.S[,ncol.S])
  m = cbind(m.X,m.S)
  m.s = m[sample(1:nrow(m), samples, replace=F),]

  cl    = makePSOCKcluster(cores)
  clusterSetRNGStream(cl)
  registerDoParallel(cl)
  s = round(seq (1, samples, length.out = cores+1))
  pst.mean = apply(m, 2, mean)

  run = foreach(chain_id = 1:cores, .packages = "BayesTSM") %dopar% {
    m.s.cores = m.s[s[chain_id]:(s[chain_id+1]-1),]
    out = apply(m.s.cores,1, function(x) obsLL(x, dat=mod$dat, Z.X = mod$Z.X, Z.S = mod$Z.S, dist.X=mod$dist.X,
                                               dist.S=mod$dist.S, log.scale = F, sum.ll=F) )
  }
  stopCluster(cl)

  run = do.call(cbind, run)

  lppd = sum (log(apply(run,1, mean)))
  q1   = sum (apply( log(run), 1, mean))
  q2   = sum( apply( log(run), 1, var) )
  q3   = mean( apply( log(run),2, sum) )
  q4   = sum(obsLL(pst.mean, dat=mod$dat, Z.X = mod$Z.X, Z.S = mod$Z.S, dist.X=mod$dist.X,
                   dist.S=mod$dist.S, log.scale = T, sum.ll=F))
  DIC  = -2* ( q4 - 2*(q4-q3) )
  WAIC1    = -2*(-lppd + 2*q1)
  WAIC2    = -2*(lppd - q2)
  mat=matrix(nrow=1, ncol=3)
  mat[1,]=c(WAIC1, WAIC2, DIC)
  colnames(mat) = c('WAIC1', 'WAIC2', 'DIC')
  mat
}

