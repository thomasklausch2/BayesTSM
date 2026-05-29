#' Return lop complete data posterior with one component updated, for slice sampler
#'
#' @noRd
log_aft_compdatpst_component <- function(v, j, eta, y, Z,
                                         tau = 4, sig.prior = 1, dist = 'weibull', beta.prior = 't',
                                         log_prior_fun){
  if(j > length(eta)) stop('j greater than length of eta')
  eta[j] = v
  log_aft_compdatpst(eta, y, Z, tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
}
