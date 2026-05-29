#' Bind two \code{mcmc.list} objects by chain
#'
#' Helper function to row-bind corresponding chains from two MCMC lists.
#'
#' @param pr An \code{mcmc.list} object.
#' @param r An \code{mcmc.list} object with the same number of chains as \code{pr}.
#'
#' @return An \code{mcmc.list} object in which each chain is formed by
#'   row-binding the corresponding chains from \code{pr} and \code{r}.
#'
#' @keywords internal
#' @noRd
bind.mcmclists = function(pr, r){
  l = length(pr)
  b = list()
  b[[1]] = mcmc(rbind ( as.mcmc(pr[1]), as.mcmc(r[1]) ))
  if(l > 1){
  for( i in 2:l){
    b[[i]] = mcmc(rbind ( as.mcmc(pr[i]), as.mcmc(r[i]) ) )
  }}

  mcmc.list(b)
}
