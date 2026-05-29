#' Sample one new parameter using unidimensional slice sampling
#'
#' @noRd
sample_slice_component = function(j, eta, w, y, Z,
                                  tau = 4, sig.prior = 1, dist = 'weibull', beta.prior = 't',
                                  log_prior_fun){
  v0 = eta[j]
  logy = log_aft_compdatpst_component(v = v0, j = j,
                                      eta = eta, y, Z,
                                      tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun) +
    log( runif(1) )

  l = v0 - w*runif(1)
  r = l+w

  logl = log_aft_compdatpst_component(v = l, j = j,
                                      eta = eta, y, Z,
                                      tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  while(logl > logy){
    l = l - w
    logl = log_aft_compdatpst_component(v = l, j = j,
                                        eta = eta, y, Z,
                                        tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  }

  logr = log_aft_compdatpst_component(v = r, j = j,
                                      eta = eta, y, Z,
                                      tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  while(logr > logy){
    r = r + w
    logr = log_aft_compdatpst_component(v = r, j = j,
                                        eta = eta, y, Z,
                                        tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  }

  v1 <- runif(1, l, r)
  logv1 = log_aft_compdatpst_component(v = v1, j = j,
                                       eta = eta, y, Z,
                                       tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  while(logv1 < logy){
    if(v1 < v0) l <- v1 else  r <- v1
    v1 <- runif(1, l, r)
    logv1 = log_aft_compdatpst_component(v = v1, j = j,
                                         eta = eta, y, Z,
                                         tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  }
  return(v1)
}
