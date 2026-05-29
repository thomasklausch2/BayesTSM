#' Plot method for bayestsm objects
#'
#' @param x An object of class `"bayestsm"`.
#' @param warmup Number of initial MCMC iterations to discard.
#' @param ... Additional arguments passed to plotting functions.
#' @param plot.X Logical; whether trace plots for the parameters of the x-model should be made.
#' @param plot.S Logical; whether trace plots for the parameters of the s-model should be made.
#'
#' @method plot bayestsm
#' @export
plot.bayestsm <- function(x, warmup = 0, plot.X = TRUE, plot.S = TRUE, ...) {
  mod <- x
  if (!inherits(x, "bayestsm")) {
    stop("`x` must be of class 'bayestsm'.", call. = FALSE)
  }
  if(warmup>0){
  par.X <- trim.mcmc(mod$par.X, burnin = warmup)
  par.S <- trim.mcmc(mod$par.S, burnin = warmup)
  } else{
    par.X <- mod$par.X
    par.S <- mod$par.S
  }
  n_iter <- mod$convergence$n_iter
  if(n_iter > 2e4) {
    thinning = round(n_iter/ 2e4)
    par.X <- trim.mcmc(par.X, thinning = thinning)
    par.S <- trim.mcmc(par.S, thinning = thinning)
  }
  par.X <- mcmc.list( lapply(par.X, function(x) {
    colnames(x) <- paste( colnames( x ),  'of X')
    mcmc(x)
  }
  )
  )
  par.S <- mcmc.list( lapply(par.S, function(x) {
    colnames(x) <- paste( colnames( x ),  'of S')
    mcmc(x)
  }
  )
  )

  if(plot.X) graphics::plot(par.X)
  if(plot.S) graphics::plot(par.S)
}
