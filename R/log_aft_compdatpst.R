#' Log posterior of an AFT model
#' @keywords internal
#' @noRd
log_aft_compdatpst = function(eta, y, Z, tau = 4,
                              sig.prior = 1, dist = 'weibull', beta.prior = 't',
                              log_prior_fun){

  LL     = LL_aft(eta, y, Z, dist)
  log_pi = log_prior_fun (
    eta = eta,
    tau = tau,
    sig.prior = sig.prior,
    beta.prior = beta.prior
  )

  LL + log_pi
}
