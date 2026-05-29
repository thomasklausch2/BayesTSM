#' Trim and thin an `mcmc.list`
#'
#' Convenience function for trimming burn-in iterations and applying thinning to
#' an object of class `mcmc.list`.
#'
#' @param obj An object of class `mcmc.list`.
#' @param burnin Integer; first iteration to retain. Defaults to `1`.
#' @param end Integer; last iteration to retain. Defaults to the number of rows
#'   in the first chain.
#' @param thinning Integer; thinning interval. Defaults to `1`.
#'
#' @return An object of class `mcmc.list` containing the trimmed and thinned
#'   chains.
#'
#' @details
#' The function subsets each chain in `obj` to the iterations
#' `seq(burnin, end, thinning)` and reconstructs the result as an `mcmc.list`.
#'
#' Note that the argument name is `thinning` to match the existing function
#' definition.
#' @export
trim.mcmc = function(obj, burnin = 1, end = nrow (as.matrix(obj[1])), thinning = 1){
  mcmc.list(lapply(
      obj,
      function(x) mcmc(
                       x[seq(burnin,end,thinning),],
                       start = burnin,
                       end=end,
                       thin = thinning
                       )
    )
    )
}
