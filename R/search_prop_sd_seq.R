#' Heuristic search for a Metropolis-Hastings proposal standard deviation (sequential processing)
#'
#' Tunes the Metropolis-Hastings proposal standard deviation by iteratively
#' extending a [bayestsm] fitted model run and adjusting the proposal scale until the
#' observed acceptance rate falls within a target interval for a specified
#' number of successive iterations.
#'
#' @param m A fitted model object containing acceptance information and current
#'   proposal settings. This object is updated iteratively during the search.
#' @param mc Integer; number of additional MCMC iterations used at each tuning
#'   step. Defaults to `3000`.
#' @param succ.min Integer; number of successive successful tuning iterations
#'   required before the search stops. Defaults to `3`.
#' @param acc_bounds Numeric vector of length 2 giving the lower and upper
#'   target bounds for the acceptance rate of the `X` parameter update. Defaults
#'   to `c(0.2, 0.3)`.
#' @param silent When output should be produced (FALSE) or not (TRUE).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`prop.sd.X`}{The tuned proposal standard deviation for the `X`
#'   update.}
#'   \item{`ac.X`}{The final acceptance rate used for the stopping decision.}
#' }
#'
#' @details
#' Like [search_prop_sd()] but instead of parallel processing through [bayestsm] sequential processing with a `for` loop through [bayestsm_seq] is done. This can be helpful in Monte Carlo simulation studies where parallelization within `bayestsm` is not desirable.
#' @export
search_prop_sd_seq = function(m, mc = 3000, succ.min = 3,
                              acc_bounds =c(0.21,0.24), silent = F){
  if (!inherits(m, "bayestsm")) {
    stop("`object` must be of class 'bayestsm'.", call. = FALSE)
  }
  found.X = F
  it = 1; succ = 0;
  while(succ!=succ.min){
    if(!silent) message(sprintf('Iteration %d',it))
    if(it == 1) { ac = mean(m$ac)
    prop.sd = m$prop.sd}
    if(it >1) {
      m = bayestsm_seq( prev.run = m, mc = mc, prop.sd = prop.sd, silent = T)
      ac = mean(m$ac)
    }
    acc.bounds.mean = (acc_bounds[2]-acc_bounds[1])/2 + acc_bounds[1]
    if(!silent) message(sprintf('Averaged acceptance rate: %g', round(ac,3)))

    if((ac > acc_bounds[1] & ac < acc_bounds[2]) ){
      found.X = T
      if(!silent) message(sprintf('Success with sd=%g',round(prop.sd,3)))
    } else{
      if(ac < acc_bounds[1]) {
        dif =  1-(acc.bounds.mean - ac)/acc.bounds.mean
        prop.sd=prop.sd * dif
      }
      if(ac > acc_bounds[2]) {
        dif =  1+(ac - acc.bounds.mean)/acc.bounds.mean
        prop.sd=prop.sd * dif
      }
      found.X = F
      if(!silent) message(sprintf('sd set to %g', round(prop.sd, 3)))
    }
    it=it+1
    if(found.X ){ mc = mc*2; succ = succ +1; found.X =  F
    if(succ!=succ.min & !silent) message(sprintf('Repeating for %d posterior draws.', mc))}
  }
  message(sprintf('Search completed.'))
  ret= list()
  ret$prop.sd = prop.sd
  ret$ac      = ac
  ret
}
