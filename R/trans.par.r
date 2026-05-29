#' Internal parameter transformations for latent time models
#'
#' Helper functions that transform regression-scale parameters to
#' distribution-specific parameters used internally by the latent time models.
#'
#' @param p1 Numeric vector of regression coefficients for the first latent time
#'   model.
#' @param p2 Numeric vector of regression coefficients for the second latent time
#'   model.
#' @param v1 Numeric scalar giving the log-scale parameter for the first latent
#'   time model.
#' @param v2 Numeric scalar giving the log-scale parameter for the second latent
#'   time model.
#' @param cr Numeric scalar; correlation parameter.
#' @param Z1 Numeric design matrix for the first latent time model.
#' @param Z2 Numeric design matrix for the second latent time model.
#' @param a Numeric vector of conditioning values on the log scale.
#' @param p Numeric vector of regression coefficients.
#' @param v Numeric scalar giving the log-scale parameter.
#' @param par Numeric vector containing regression coefficients followed by a
#'   log-scale parameter.
#'
#' @return
#' `trans.par.norm()` returns a two-column numeric matrix with conditional means
#' and conditional standard deviations under the normal-error model.
#' `trans.par.ind.norm()` returns a two-column numeric matrix with means and
#' standard deviations under the independent normal-error model.
#' `trans.par()` returns a two-column numeric matrix containing transformed
#' distribution parameters for the extreme-value or logistic-error model.
#'
#' @details
#' These functions are used internally to convert linear predictor and log-scale
#' parameterizations to the parameterizations required by the corresponding
#' latent time distributions.
#'
#' `trans.par.norm()` computes conditional normal parameters for one latent time
#' given the other under correlated normal errors.
#'
#' `trans.par.ind.norm()` computes mean and standard deviation parameters for an
#' independent normal-error model.
#'
#' `trans.par()` transforms regression coefficients and a log-scale parameter to
#' positive distribution parameters for the extreme-value and logistic-error
#' models.
#'
#' @name trans_par
#' @keywords internal
#' @noRd
trans.par.norm = function(p1, p2, v1, v2, cr, Z1, Z2, a){
  mu1 = Z1 %*% as.matrix(p1)
  mu2 = Z2 %*% as.matrix(p2)
  sd1 = exp(v1)
  sd2 = exp(v2)
  a  = logrob(a, tol=1e-300)
  m  = mu1 + sd1/sd2*cr*(a-mu2)
  v  = sqrt( (1-cr^2)*sd1^2 )
  cbind(m,v)
}

#' @rdname trans_par
#' @keywords internal
#' @noRd
trans.par.ind.norm = function(p, v, Z1){
  mu = Z1 %*% as.matrix(p)
  sd = exp(v)
  cbind(mu,sd)
}

#' @rdname trans_par
#' @keywords internal
#' @noRd
trans.par = function(Z1, par){
  p = length(par)
  p1 = exp(Z1 %*% as.matrix(par[1:(p-1)]))
  p2 = 1/exp(par[(p)])
  cbind(p1,p2)
}
