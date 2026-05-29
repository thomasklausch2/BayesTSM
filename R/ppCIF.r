#' Posterior Predictive CIFs for \eqn{X}, \eqn{S}, and \eqn{Y=X+S}
#'
#' Computes posterior predictive cumulative incidence functions (CIFs) for the
#' latent event times \eqn{X}, \eqn{S}, and \eqn{Y=X+S}. In the present
#' progressive three-state setting, these CIFs are identical to the posterior
#' predictive cumulative distribution functions (CDFs) of the corresponding
#' latent times.
#'
#' @param mod A fitted \code{bayestsm} model object.
#' @param fix_Z.X Optional numeric vector used to fix covariates in the
#'   \eqn{X}-submodel. Its length must equal \code{ncol(mod$Z.X)}. Entries correspond to the columns of \code{mod$Z.X}.
#'   Values set to \code{NA} mean that the corresponding column (variable) in \code{Z.X} is not fixed and, thus, empirically
#'   marginalized out. A valued entry specified by the user means that the CIF \code{X} is estimated conditional on that value.
#' @param fix_Z.S Optional numeric vector used to fix covariates in the
#'   \eqn{S}-submodel. See \code{fix_Z.X}.
#' @param pst.samples Posterior samples used to compute posterior predictive CIFs.
#' @param warmup The number of iterations discarded for warmup per chain. The default is `NULL` in which case, conservatively, the first half of posterior samples is discarded.
#' @param perc Numeric vector of probabilities in \eqn{[0,1)} used when
#'   \code{type = "percentiles"}. For each posterior draw, posterior predictive
#'   quantiles are evaluated at these probabilities.
#' @param quant Optional numeric vector of time points used when
#'   \code{type = "quantiles"}. At these values, the posterior predictive CDF is
#'   evaluated. If \code{NULL}, a grid is constructed automatically from
#'   \code{0} to the largest finite follow-up time in the data.
#' @param method Character string specifying the computational method.
#'   Currently supported values are \code{"analytic"} and
#'   \code{"simulation"}. See details.
#' @param type Character string specifying whether posterior predictive
#'   inference is returned as CDF values on a time grid
#'   (\code{type = "quantiles"}) or as quantiles on a probability grid
#'   (\code{type = "percentiles"}). The latter is available only with
#'   \code{method = "simulation"}.
#' @param cred_int The vector of the lower and upper bound of the credible (posterior)
#' intervals returned. Default (2.5%, 97.5%): c(0.025, 0.975).
#'
#' @details
#' Let \eqn{\theta^{(m)}} denote posterior draw \eqn{m=1,\dots,M} and let
#' \eqn{z_i} be the covariate vector for individual \eqn{i=1,\dots,n}. For a
#' given latent variable \eqn{T \in \{X,S\}}, the posterior predictive CDF at
#' covariate pattern \eqn{z_i} is
#' \deqn{
#' F_T(t \mid z_i, \theta^{(m)}) = P(T \le t \mid z_i, \theta^{(m)}).
#' }
#' If covariates are not fully fixed, the function returns a covariate-marginal
#' posterior predictive CDF obtained by averaging over the empirical
#' distribution of the covariates in the fitted sample:
#' \deqn{
#' \bar F_T(t \mid \theta^{(m)}) = \frac{1}{n}\sum_{i=1}^n F_T(t \mid z_i, \theta^{(m)}).
#' }
#' This yields a posterior sample of predictive CIFs, from which pointwise
#' posterior medians and credible intervals are reported.
#'
#' Two computational methods (\code{method}) are available which in practice give very similar results, provided \code{pst.samples} is chosen large.
#'
#' With \code{method = "analytic"}, closed-form evaluation of the CDF is used
#' whenever available. For \eqn{T \in \{X,S\}} and a grid of time points
#' \eqn{q_1,\dots,q_K}, this computes
#' \deqn{
#' \bar F_T(q_k \mid \theta^{(m)}) = \frac{1}{n}\sum_{i=1}^n F_T(q_k \mid z_i, \theta^{(m)}).
#' }
#' For \eqn{Y=X+S}, no closed-form expression is available in general, so its
#' posterior predictive CDF is still approximated by Monte Carlo simulation
#' under \code{method = "analytic"}.
#'
#' With \code{method = "simulation"}, posterior predictive samples are generated
#' for all requested quantities. For each posterior draw \eqn{\theta^{(m)}},
#' draws
#' \deqn{
#' X_i^{(m)} \sim p(X \mid z_i, \theta^{(m)}), \qquad
#' S_i^{(m)} \sim p(S \mid z_i, \theta^{(m)}), \qquad
#' Y_i^{(m)} = X_i^{(m)} + S_i^{(m)}
#' }
#' are generated and empirical CDFs or empirical quantiles are computed from
#' these simulated samples.
#'
#' The argument \code{type} determines whether inference is returned on a grid
#' of time points or on a grid of probabilities.
#'
#' If \code{type = "quantiles"}, the function evaluates the posterior
#' predictive CDF on the supplied grid \code{quant = (q_1,\dots,q_K)} and
#' returns
#' \deqn{
#' \bar F_T(q_k \mid \theta^{(m)}), \qquad k=1,\dots,K.
#' }
#' Thus, despite the name, \code{type = "quantiles"} means that the grid
#' consists of quantiles or time points at which the CDF is evaluated.
#'
#' If \code{type = "percentiles"}, the function evaluates the posterior
#' predictive quantile function on the supplied probability grid
#' \code{perc = (p_1,\dots,p_K)} and returns
#' \deqn{
#' Q_T(p_k \mid \theta^{(m)}) = \inf\{t : F_T(t \mid \theta^{(m)}) \ge p_k\},
#' \qquad k=1,\dots,K.
#' }
#' This option is available only with \code{method = "simulation"}.
#'
#' The arguments \code{fix_Z.X} and \code{fix_Z.S} allow selective fixing of
#' covariates when computing posterior predictive CIFs. These vectors must have
#' the same lengths as the numbers of columns in \code{mod$Z.X} and
#' \code{mod$Z.S}, respectively. Entries equal to \code{NA} indicate that the
#' corresponding covariate should remain unchanged and therefore still be
#' averaged over empirically. Non-missing entries replace the corresponding
#' covariate values for all individuals before computing the predictive CIF.
#'
#' This makes it possible to obtain subgroup-specific posterior predictive CIFs.
#' For example, a specific age value can be fixed while all remaining covariates
#' are left at \code{NA}, yielding a posterior predictive CIF for that age group
#' while marginalizing over the empirical distribution of the other covariates.
#' Formally, if the covariate vector is partitioned into fixed components
#' \eqn{z^{\mathrm{fix}}} and unfixed components \eqn{z_i^{\mathrm{free}}},
#' the function computes
#' \deqn{
#' \bar F_T(t \mid z^{\mathrm{fix}}, \theta^{(m)})
#' = \frac{1}{n}\sum_{i=1}^n
#' F_T\left(t \mid z^{\mathrm{fix}}, z_i^{\mathrm{free}}, \theta^{(m)}\right).
#' }
#' If all covariates are fixed, no marginalization over covariates remains and a
#' fully conditional posterior predictive CIF is obtained:
#' \deqn{
#' F_T(t \mid z^{\ast}, \theta^{(m)}).
#' }
#'
#' @return A list with components depending on \code{type}.
#'
#' If \code{type = "quantiles"}, the returned list contains:
#' \describe{
#'   \item{\code{med.p.x}, \code{med.p.s}, \code{med.p.y}}{Pointwise posterior
#'   medians of the predictive CDFs of \eqn{X}, \eqn{S}, and \eqn{Y=X+S} on the
#'   grid \code{grid}.}
#'   \item{\code{p.x.ci}, \code{p.s.ci}, \code{p.y.ci}}{Pointwise 95\%
#'   posterior credible intervals for the predictive CDFs.}
#'   \item{\code{grid}}{The grid of time points at which the predictive CDFs
#'   were evaluated.}
#'   \item{\code{type}}{The supplied \code{type} argument.}
#' }
#'
#' If \code{type = "percentiles"}, the returned list contains:
#' \describe{
#'   \item{\code{med.q.x}, \code{med.q.s}, \code{med.q.y}}{Pointwise posterior
#'   medians of the predictive quantile functions of \eqn{X}, \eqn{S}, and
#'   \eqn{Y=X+S} on the probability grid \code{grid}.}
#'   \item{\code{q.x.ci}, \code{q.s.ci}, \code{q.y.ci}}{Pointwise 95\%
#'   posterior credible intervals for the predictive quantile functions.}
#'   \item{\code{grid}}{The grid of probabilities at which the predictive
#'   quantile functions were evaluated.}
#'   \item{\code{type}}{The supplied \code{type} argument.}
#' }
#'
#' @seealso [bayestsm()] [plot.ppCIF()]
#'
#' @examples
#' \dontrun{
#' # Obtain quantiles for provided percentiles
#' # For a bayestsm model output saved in mod_slice (see ?bayestsm examples), run:
#' postCDF_perc = ppCIF(mod_slice, warmup = 500, perc = seq(0, 0.99, 0.01),
#'                      method = 'simulation', type = 'percentiles')
#'
#' # straight forward plot over provided percentiles - quantiles combintation
#' plot(postCDF_perc, xlim = c(0,40))
#'
#' # Change the credible interval returned to (1%, 99%)
#' postCDF_perc2 = ppCIF(mod_slice, warmup = 500, perc = seq(0, 0.99, 0.01),
#'                       method = 'simulation', type = 'percentiles',
#'                       cred_int = c(0.01, 0.99))
#' plot(postCDF_perc2, xlim = c(0,40))
#'
#' # Obtain percentiles for provided quantiles (default: grid between 0 and max. follow up)
#' postCDF_quant = ppCIF(mod_slice, warmup = 500,
#'                       method = 'simulation', type = 'quantiles')
#'
#' # plot similar to plot of postCDF_perc
#' plot(postCDF_quant, xlim = c(0,40))
#'
#' # Obtain percentiles for specific quantiles (e.g. probability for transition by 5 time units)
#' postCDF_quant2 = ppCIF(mod_slice, warmup = 500, quant = c(5, 10),
#'                        method = 'simulation', type = 'quantiles')
#'
#' postCDF_quant2
#'
#'
#' Alternatively, method = analytic can be used for quantiles with very similar results
#' postCDF_quant3 = ppCIF(mod_slice, warmup = 500,
#'                        method = 'analytic', type = 'quantiles')
#'
#' plot(postCDF_quant3, xlim = c(0,40))
#'
#' }
#'
#' @export
ppCIF = function(mod, fix_Z.X = NULL, fix_Z.S = NULL, pst.samples = 1e3, warmup = NULL,
                 perc = seq(0, 0.99, 0.01), quant = NULL, method = c('simulation','analytic'),
                 type = c('percentiles','quantiles'), cred_int = c(0.025, 0.975)) {

  # Messages and checks
  if (!inherits(mod, "bayestsm")) {
    stop("`mod` must be of class 'bayestsm'.", call. = FALSE)
  }
  method <- match.arg(method)
  type <- match.arg(type)
  if(pst.samples == 1) stop('Posterior samples have to be > 1, and ideally large.')
  if(method == "analytic")   {
    if(type == 'percentiles') stop('Inference via type = percentiles is not currently available with method = analytic.')
    message('Please note that with method = analytic the posterior predictive CDF of X+S=Y is still obtained by Monte Carlo simulation, because it is not in general anlytically tractable.')
  }
  if(method == "simulation") message('Obtaining the posterior predictive CDFs of X, S, and X+S=Y by Monte Carlo simulation.')

  # Covariate matrix preprocessing
  if( !is.null(mod$Z.X) ){
    Z.X = as.matrix(mod$Z.X)
    if( !is.null(fix_Z.X) ){
      if(length(fix_Z.X) != ncol(Z.X)) stop('Length of fix_Z.X has to be equal to ncol(Z.X)')
      Z.X[,!is.na(fix_Z.X)] = fix_Z.X[!is.na(fix_Z.X)]
    }
    Z1.X = cbind(1, Z.X)
  } else{
    if(!is.null(fix_Z.X)) warning('fix_Z.X is not null, while an intercept only model was passed. Ignoring fix_Z.X.')
    Z.X = NULL
    Z1.X = as.matrix(rep(1, nrow(mod$dat)))
  }

  if( !is.null(mod$Z.S) ){
    Z.S = as.matrix(mod$Z.S)
    if( !is.null(fix_Z.S) ){
      if(length(fix_Z.S) != ncol(Z.S)) stop('Length of fix_Z.S has to be equal to ncol(Z.S)')
      Z.S[,!is.na(fix_Z.S)] = fix_Z.S[!is.na(fix_Z.S)]
    }
    Z1.S = cbind(1, Z.S)
  } else{
    if(!is.null(fix_Z.S)) warning('fix_Z.S is not null, while an intercept only model was passed. Ignoring fix_Z.S.')
    Z.S = NULL
    Z1.S = as.matrix(rep(1, nrow(mod$dat)))
  }

  # Parameter matrix pre-processing
  it = nrow(mod$par.X[[1]])
  if(is.null(warmup)) {
    warmup = round(0.5 * it)
    message('Warmup not specified. Discarding the first half of posterior samples per chain.')
  }
  par.list.X <- trim.mcmc(mod$par.X, burnin = warmup)
  par.list.S <- trim.mcmc(mod$par.S, burnin = warmup)
  par.X <- do.call(rbind, lapply(par.list.X, as.matrix))
  par.S <- do.call(rbind, lapply(par.list.S, as.matrix))

  p.X <- ncol(par.X)
  p.S <- ncol(par.S)

  beta.X <- par.X[, seq_len(p.X - 1), drop = FALSE]
  beta.S <- par.S[, seq_len(p.S - 1), drop = FALSE]
  sigma.X <- par.X[, p.X]
  sigma.S <- par.S[, p.S]
  if(pst.samples > nrow(beta.X)) stop('pst.samples is chosen greater than the total number of posterior samples available. Update the model first to obtain more posterior samples.')

  # Get model distributions
  dist.X = mod$dist.X
  dist.S = mod$dist.S

  # Define quant, if null and needed
  if (is.null(quant) && type == 'quantiles'){
    message('Quantiles were not provided and therefore chosen automatically.')
    max.fu = max(mod$dat$R[is.finite(mod$dat$R)])
    quant = seq(0, max.fu, max.fu/100)
  }

  if (type == 'percentiles' && (1 %in% perc)){
    message('perc contained value(s) numerically equal to one (removed).')
    perc = perc[- which(perc == 1)]
  }

  # Define grid
  if(type == "percentiles") grid = perc
  if(type == "quantiles") grid = quant

  # Draw posterior samples
  s <- sample.int(nrow(par.X), pst.samples, replace = FALSE)

  Q <- get_pp_stats(
    beta.X = beta.X[s, , drop = FALSE],
    beta.S = beta.S[s, , drop = FALSE],
    sigma.X = sigma.X[s],
    sigma.S = sigma.S[s],
    Z1.X = Z1.X,
    Z1.S = Z1.S,
    method = method,
    type = type,
    dist.X = dist.X,
    dist.S = dist.S,
    grid = grid
  )

  # Posterior CDF statistics
  ret = list()

  if(type == "percentiles"){
    if(length(perc)>1){
    ret$med.q.x = apply(Q$q.X, 1, median)
    ret$med.q.s = apply(Q$q.S, 1, median)
    ret$med.q.y = apply(Q$q.Y, 1, median)
    ret$q.x.ci  = apply(Q$q.X, 1, quantile, cred_int)
    ret$q.s.ci  = apply(Q$q.S, 1, quantile, cred_int)
    ret$q.y.ci  = apply(Q$q.Y, 1, quantile, cred_int)
    } else{
      ret$med.q.x = median(Q$q.X)
      ret$med.q.s = median(Q$q.S)
      ret$med.q.y = median(Q$q.Y)
      ret$q.x.ci  = quantile(Q$q.X, cred_int)
      ret$q.s.ci  = quantile(Q$q.S, cred_int)
      ret$q.y.ci  = quantile(Q$q.Y, cred_int)
    }
  }

  if(type == "quantiles"){
    if(length(quant)>1){
    ret$med.p.x = apply(Q$p.X, 1, median)
    ret$med.p.s = apply(Q$p.S, 1, median)
    ret$med.p.y = apply(Q$p.Y, 1, median)
    ret$p.x.ci  = apply(Q$p.X, 1, quantile, cred_int)
    ret$p.s.ci  = apply(Q$p.S, 1, quantile, cred_int)
    ret$p.y.ci  = apply(Q$p.Y, 1, quantile, cred_int)
    } else{
      ret$med.p.x = median(Q$p.X)
      ret$med.p.s = median(Q$p.S)
      ret$med.p.y = median(Q$p.Y)
      ret$p.x.ci  = quantile(Q$p.X, cred_int)
      ret$p.s.ci  = quantile(Q$p.S, cred_int)
      ret$p.y.ci  = quantile(Q$p.Y, cred_int)
    }
  }

  ret$grid = grid
  ret$type = type
  class(ret) = "ppCIF"
  ret
}
