#' Joint posterior density for `X` and `S` accelerated failure time models
#'
#' @keywords internal
#' @noRd
log_aft_compdatpst_XS <- function(
    eta, X, S,
    d23, Z.X, Z.S,
    tau.X, sig.prior.X,
    tau.S, sig.prior.S,
    dist.X, dist.S,
    beta.prior = "t",
    log_prior_fun
) {
  stopifnot(is.logical(d23), length(d23) == length(S))

  p1.X <- ncol(Z.X)
  p1.S <- ncol(Z.S)

  par.X <- eta[1:(p1.X + 1)]
  par.S <- eta[(p1.X + 2):(p1.X + p1.S + 2)]

  pst.X <- log_aft_compdatpst(
    eta = par.X,
    y = X,
    Z = Z.X,
    tau = tau.X,
    sig.prior = sig.prior.X,
    dist = dist.X,
    beta.prior = beta.prior,
    log_prior_fun = log_prior_fun
  )

  pst.S <- log_aft_compdatpst(
    eta = par.S,
    y = S[d23],
    Z = Z.S[d23, , drop = FALSE],
    tau = tau.S,
    sig.prior = sig.prior.S,
    dist = dist.S,
    beta.prior = beta.prior,
    log_prior_fun = log_prior_fun
  )

  pst.X + pst.S
}
