#' Heuristic search for a Metropolis-Hastings proposal standard deviation
#'
#' Tunes the Metropolis-Hastings proposal standard deviation by iteratively
#' extending a bayestsm fitted model run and adjusting the proposal scale until the
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
#'   \item{`prop.sd`}{The tuned proposal standard deviation for the `X`
#'   update.}
#'   \item{`ac`}{The final acceptance rate used for the stopping decision.}
#' }
#'
#' @details
#' The function implements a simple stepwise heuristic for proposal tuning. It
#' starts from the current proposal standard deviation stored in `m$prop.sd`
#' and compares the observed Metropolis-Hastings acceptance rate with the target
#' interval `acc_bounds`.
#'
#' If the acceptance rate is below the lower target bound, the proposal standard
#' deviation is decreased. If the acceptance rate is above the upper target
#' bound, the proposal standard deviation is increased. Once the acceptance rate
#' falls within the target interval, the current tuning step is counted as a
#' success.
#'
#' After each successful step, the number of MCMC iterations used for the next
#' update is doubled. The search stops after `succ.min` successive successful
#' steps. This can help stabilize the tuning decision by checking the acceptance
#' behaviour over progressively longer runs.
#'
#' The function updates the model by calling `bayestsm()` with `prev.run = m`.
#' As currently implemented, only the proposal standard deviation for the `X`
#' update is tuned.
#'
#' This is a heuristic tuning tool intended to provide a reasonable proposal
#' scale before longer production runs. It does not guarantee optimal mixing or
#' efficiency.
#' @export
search_prop_sd = function(m, mc = 3000, succ.min = 3, acc_bounds =c(0.21,0.24), silent = F){
  if (!inherits(m, "bayestsm")) {
    stop("`m` must be of class 'bayestsm'.", call. = FALSE)
  }
  found.X = F
  it = 1; succ = 0;
  while(succ!=succ.min){
    if(!silent) message(sprintf('Iteration %d',it))
    if(it == 1) { ac = mean(m$ac)
                  prop.sd = m$prop.sd}
    if(it >1) {
      m = bayestsm( prev.run = m, mc = mc, prop.sd = prop.sd, silent = T)
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

