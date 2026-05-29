#' Run one Gibbs sampler chain for the progressive three-state model
#'
#' @noRd
full_conditional_Gibbs <- function(mc, d, d23, L, R, Z1.X, Z1.S,
                                p.X, p.S, p1.X, p1.S,
                                x.ini, s.ini,
                                beta.x.ini, sigma.x.ini,
                                beta.s.ini, sigma.s.ini,
                                prop.sd,
                                dist.X, dist.S,
                                beta.prior.X, beta.prior.S,
                                sig.prior.X, sig.prior.S,
                                fix.sigma.X, fix.sigma.S,
                                beta.prior,
                                MH = TRUE, slicesampler_stepsize = 1,
                                log_prior_fun) {
  if (length(prop.sd) == 1L) {
    prop.var.mat <- diag(prop.sd^2, p1.X + p1.S + 2)
  } else {
    prop.var.mat <- prop.sd
  }

  cur.par.Xreg <- matrix(NA_real_, ncol = p.X + 2, nrow = mc + 1)
  cur.par.Xreg[1, 1:(p.X + 1)] <- beta.x.ini
  cur.par.Xreg[1, p.X + 2] <- sigma.x.ini

  cur.par.Sreg <- matrix(NA_real_, ncol = p.S + 2, nrow = mc + 1)
  cur.par.Sreg[1, 1:(p.S + 1)] <- beta.s.ini
  cur.par.Sreg[1, p.S + 2] <- sigma.s.ini

  X <- x.ini
  S <- s.ini
  ac <- rep(NA, mc + 1)
  ac[1] <- 1

  for (i in seq_len(mc)) {
    par <- t(as.matrix(c(cur.par.Xreg[i, , drop = FALSE],
                         cur.par.Sreg[i, , drop = FALSE])))

    if(MH){
      mh <- mhstep_aft_XS(
        x = par,
        X = X,
        S = S,
        d23 = d23,
        Z.X = Z1.X,
        Z.S = Z1.S,
        tau.X = beta.prior.X,
        sig.prior.X = sig.prior.X,
        tau.S = beta.prior.S,
        sig.prior.S = sig.prior.S,
        dist.X = dist.X,
        dist.S = dist.S,
        prop.var = prop.var.mat,
        fix.sigma.S = fix.sigma.S,
        fix.sigma.X = fix.sigma.X,
        beta.prior = beta.prior,
        log_prior_fun = log_prior_fun
      )

      cur.par.Xreg[i + 1, ] <- mh$s[1:(p1.X + 1)]
      cur.par.Sreg[i + 1, ] <- mh$s[(p1.X + 2):(p1.X + p1.S + 2)]
      ac[i + 1] <- mh$ac
    } else{

      cur.par.Xreg[i + 1,] <- slice_passthrough(eta  = cur.par.Xreg[i, , drop = FALSE],
                                                w    = slicesampler_stepsize,
                                                y    = X,
                                                Z    = Z1.X,
                                                tau  = beta.prior.X,
                                                sig.prior = sig.prior.X,
                                                dist =  dist.X,
                                                beta.prior = beta.prior,
                                                log_prior_fun = log_prior_fun,
                                                fix.sigma = fix.sigma.X
                                                )

      cur.par.Sreg[i + 1,] <- slice_passthrough(eta  = cur.par.Sreg[i, , drop = FALSE],
                                                w    = slicesampler_stepsize,
                                                y    = S[d23],
                                                Z    = Z1.S[d23, , drop = FALSE],
                                                tau  = beta.prior.S,
                                                sig.prior = sig.prior.S,
                                                dist =  dist.S,
                                                beta.prior = beta.prior,
                                                log_prior_fun = log_prior_fun,
                                                fix.sigma = fix.sigma.S)

    }

    if (dist.X != "lognormal") {
      cur.par.X <- trans.par(Z1 = Z1.X, par = cur.par.Xreg[i + 1, ])
    } else {
      cur.par.X <- trans.par.ind.norm(
        Z1 = Z1.X,
        p = cur.par.Xreg[i + 1, 1:p1.X],
        v = cur.par.Xreg[i + 1, p1.X + 1]
      )
    }

    if (dist.S != "lognormal") {
      cur.par.S <- trans.par(Z1 = Z1.S, par = cur.par.Sreg[i + 1, ])
    } else {
      cur.par.S <- trans.par.ind.norm(
        Z1 = Z1.S,
        p = cur.par.Sreg[i + 1, 1:p1.S],
        v = cur.par.Sreg[i + 1, p1.S + 1]
      )
    }

    X <- pst.X(par = cur.par.X, S = S, d = d, L = L, R = R, dist = dist.X)
    S[d23] <- pst.S(
      par = cur.par.S[d23, , drop = FALSE],
      X = X[d23],
      d = d[d23],
      L = L[d23],
      R = R[d23],
      dist = dist.S
    )

    X[X == 0] <- 1e-300
    S[S == 0] <- 1e-300
  }

  cur.par.Xreg <- cur.par.Xreg[2:(mc + 1), , drop = FALSE]
  cur.par.Sreg <- cur.par.Sreg[2:(mc + 1), , drop = FALSE]

  list(
    cur.par.Xreg = cur.par.Xreg,
    cur.par.Sreg = cur.par.Sreg,
    X = X,
    S = S,
    ac = ac
  )
}
