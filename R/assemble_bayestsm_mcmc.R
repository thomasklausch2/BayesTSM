#' Assemble Gibbs sampler output into mcmc objects
#'
#' Converts per-chain Gibbs sampler output into `coda::mcmc.list` objects for
#' the `X` and `S` parameter draws, concatenates latent draws and acceptance
#' indicators across chains, optionally appends a previous run, and constructs
#' post-burn-in thinned output.
#'
#' @param run List of Gibbs sampler outputs, one element per chain.
#' @param Z.X Covariate matrix for the `X` transition model.
#' @param Z.S Covariate matrix for the `S` transition model.
#' @param burnin Integer burn-in used for the post-burn-in output.
#' @param thinning Integer thinning interval used for the post-burn-in output.
#' @param mc Integer number of Gibbs iterations in the current run.
#' @param prev.run Optional previous fitted `bayestsm` object.
#' @param update.burnin Logical; whether to reset `burnin` to half of the total
#'   combined iterations when `prev.run` is supplied.
#'
#' @return A list containing `par.X.all`, `par.S.all`,
#'   `X.draw`, `S.draw`, `ac`, `ac.cur`, and possibly updated `burnin`.
#'
#' @noRd
assemble_bayestsm_mcmc <- function(run,
                                   Z.X,
                                   Z.S,
                                   burnin = 1,
                                   thinning,
                                   mc,
                                   prev.run = NULL,
                                   update.burnin = TRUE,
                                   rescale_intercept) {
  n.chains <- length(run)

  mcmc.par.X <- vector("list", n.chains)
  mcmc.par.S <- vector("list", n.chains)
  ac.list    <- vector("list", n.chains)
  X.list     <- vector("list", n.chains)
  S.list     <- vector("list", n.chains)

  x.names <- c("Intercept", colnames(Z.X), "sigma")
  s.names <- c("Intercept", colnames(Z.S), "sigma")

  for (i in seq_len(n.chains)) {
    par.X <- run[[i]]$par.X
    par.S <- run[[i]]$par.S

    par.X[, ncol(par.X)] <- exp(par.X[, ncol(par.X)])
    par.S[, ncol(par.S)] <- exp(par.S[, ncol(par.S)])
    par.X[,1] <- par.X[,1] + log(rescale_intercept)
    par.S[,1] <- par.S[,1] + log(rescale_intercept)

    colnames(par.X) <- x.names
    colnames(par.S) <- s.names

    mcmc.par.X[[i]] <- mcmc(par.X)
    mcmc.par.S[[i]] <- mcmc(par.S)

    ac.list[[i]] <- run[[i]]$ac
    X.list[[i]]  <- run[[i]]$X
    S.list[[i]]  <- run[[i]]$S
  }

  par.X.all <- mcmc.list(mcmc.par.X)
  par.S.all <- mcmc.list(mcmc.par.S)

  ac <- ac.list[[1]]
  if (n.chains > 1) {
    for (i in 2:n.chains) {
      ac <- cbind(ac, ac.list[[i]])
    }
  }

  X.draw <- unlist(X.list, use.names = FALSE)
  S.draw <- unlist(S.list, use.names = FALSE)

  par.X.all <- trim.mcmc(par.X.all, burnin = 1, thinning = thinning)
  par.S.all <- trim.mcmc(par.S.all, burnin = 1, thinning = thinning)

  ac.cur <- NULL
  if (!is.null(prev.run)) {
    par.X.all <- bind.mcmclists(prev.run$par.X, par.X.all)
    par.S.all <- bind.mcmclists(prev.run$par.S, par.S.all)

    ac.cur <- ac
    ac <- rbind(prev.run$ac, ac)

    mc.prev <- nrow(as.matrix(prev.run$par.X[1]))
    if (update.burnin) {
      burnin <- round((mc.prev + mc) / 2)
    }
  }

  return(list(
    par.X.all = par.X.all,
    par.S.all = par.S.all,
    # par.X.bi = par.X.bi,
    # par.S.bi = par.S.bi,
    X.draw = X.draw,
    S.draw = S.draw,
    ac = ac,
    ac.cur = ac.cur,
    burnin = burnin
  ))
}
