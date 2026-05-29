#' Ons slice sampler passthrough through all parameters, constituting one draw of all parameters using slice sampling
#'
#' @noRd
slice_passthrough <- function(eta, w = 1, y, Z,
                              tau = 4, sig.prior = 1, dist = 'weibull', beta.prior = 't',
                              log_prior_fun, fix.sigma){

  sigma_pos = length(eta)
  for(i in seq_along(eta[1:(sigma_pos-1)])) eta[i] <-
      sample_slice_component(i, eta, w, y, Z,
                             tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  if( !fix.sigma ){
    eta[sigma_pos] <-
      sample_slice_component(sigma_pos, eta, w, y, Z,
                             tau, sig.prior, dist, beta.prior, log_prior_fun = log_prior_fun)
  }
  eta
}
