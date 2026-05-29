#' Bayesian Accelerated Failure-Time Three-State Model with Censoring After Intervention
#'
#' Estimates the Bayesian progressive three-state model for interval-censored
#' three-state screening outcomes with censoring after intervention
#' (Klausch et al., 2023). The model is fitted using a Gibbs sampler with data
#' augmentation for the latent event times `X` and `S`, and either slice-sampling
#' or Metropolis-Hastings updates for the transition-model parameters.
#'
#' @param d Numeric vector indicating the observed outcome category for each
#'   individual. Allowed values are `1`, `2`, and `3`.
#' @param L Numeric vector of left interval bounds. Minimum allowed value is `0`.
#' @param R Numeric vector of right interval bounds. Use `Inf` for right-censored
#'   observations. If `R = Inf`, `d` must be `1`. If `R < Inf`, `d` must be `2`
#'   or `3`.
#' @param Z.X Optional design matrix of covariates for the first transition model
#'   for latent time `X`.
#' @param Z.S Optional design matrix of covariates for the second transition model
#'   for latent time `S`. Defaults to `Z.X`.
#'
#' @param dist.X Character string specifying the distribution for the `X`
#'   transition model. Allowed values are `"weibull"`, `"lognormal"`, and
#'   `"loglog"`.
#' @param dist.S Character string specifying the distribution for the `S`
#'   transition model. Allowed values are `"weibull"`, `"lognormal"`, and
#'   `"loglog"`.
#'
#' @param beta.prior Character string specifying the default prior family for
#'   regression coefficients when `log_prior_fun = log_aft_prior`. Allowed values
#'   are `"t"` for a Student-\eqn{t} prior and `"norm"` for a normal prior.
#'   Default is `"t"`.
#' @param beta.prior.X Prior hyperparameter for the regression coefficients in
#'   the `X` transition model. If `log_prior_fun = log_aft_prior` and
#'   `beta.prior = "t"`, this is the degrees of freedom of the Student-\eqn{t}
#'   prior. If `beta.prior = "norm"`, this is the standard deviation of the
#'   normal prior. For custom prior functions, this value is passed as `tau`.
#'   Default is `4`.
#' @param beta.prior.S Prior hyperparameter for the regression coefficients in
#'   the `S` transition model. Interpreted in the same way as `beta.prior.X`.
#'   For custom prior functions, this value is passed as `tau`. Defaults to
#'   `beta.prior.X`.
#' @param sig.prior.X Prior hyperparameter for the scale parameter `sigma_X` of
#'   the `X` transition model. When `fix.sigma.X = FALSE`, this value is passed
#'   as `sig.prior` to the prior function for the `X` model. When
#'   `fix.sigma.X = TRUE`, it is interpreted as the fixed value of `sigma_X`.
#'   Default is `1`.
#' @param sig.prior.S Prior hyperparameter for the scale parameter `sigma_S` of
#'   the `S` transition model. Interpreted in the same way as `sig.prior.X`.
#'   When `fix.sigma.S = TRUE`, this is the fixed value of `sigma_S`. Default is
#'   `1`.
#' @param fix.sigma.X Logical; whether to fix the scale parameter `sigma_X` at
#'   `sig.prior.X`.
#' @param fix.sigma.S Logical; whether to fix the scale parameter `sigma_S` at
#'   `sig.prior.S`.
#'
#' @param mc Integer; number of Gibbs iterations per chain.
#' @param chains Integer; number of MCMC chains.
#' @param thinning Integer; thinning interval applied when constructing the
#'   post-sampling output.
#' @param warmup Numeric; number of original Gibbs iterations to omit as warmup
#'   before assessing convergence. Internally this is converted to stored
#'   post-thinning iterations using `round(warmup / thinning)`. `warmup` should be significantly smaller than `mc`; default `mc * 0.05`.
#' @param prop.sd Proposal scale for the random-walk Metropolis-Hastings update.
#'   Used only when `MH = TRUE`. This may be a proposal standard deviation or
#'   proposal covariance specification used for the parameter updates. If `NULL`
#'   and `MH = TRUE`, the function runs a short preliminary sampler and searches
#'   for a useful proposal scale automatically.
#'
#' @param update_till_convergence Logical; whether to automatically extend the
#'   MCMC run until the convergence criteria based on `min_R` and `min_eff` are
#'   met, or until `maxit` is reached.
#' @param mc_update Integer; number of Gibbs iterations per chain to add when
#'   updating a previous run or when automatic convergence updating is enabled.
#'   Defaults to `mc`.
#' @param min_R Numeric; target upper threshold for convergence based on the
#'   potential scale reduction factor. Used when `update_till_convergence = TRUE`.
#'   Default is `1.1`.
#' @param min_eff Numeric; target minimum effective sample size. Used when
#'   `update_till_convergence = TRUE`. Default is `chains * 100`.
#' @param maxit Optional integer; maximum number of automatic update cycles when
#'   `update_till_convergence = TRUE`. If `NULL`, no explicit maximum number of
#'   update cycles is imposed by this argument.
#'
#' @param prev.run Optional fitted object returned by a previous call to
#'   `bayestsm()`. If supplied, the sampler continues from the last sampled
#'   parameter values and reuses the stored data, priors, proposal settings,
#'   distributions, fixed-scale settings, and other model options from the
#'   previous run.
#'
#' @param MH Logical; whether to use random-walk Metropolis-Hastings updates for
#'   the transition-model parameters. If `FALSE`, slice-sampling updates are
#'   used. Default is `FALSE`.
#' @param rescale_times Logical; whether to rescale `L` and `R` by the median of
#'   the finite values of `R` before sampling. Returned data are stored on the
#'   original time scale. Default is `TRUE`.
#' @param slicesampler_stepsize Numeric step-out size used by the slice sampler
#'   when `MH = FALSE`. Default is `1`.
#' @param silent Logical; whether to suppress progress messages.
#' @param log_prior_fun Function used to evaluate the log-prior density for the
#'   AFT transition-model parameters. The default is `log_aft_prior`. Custom
#'   prior functions should accept the arguments `eta`, `tau`, and `sig.prior`;
#'   see [log_aft_prior()] for details.
#' @param do_cpp Logical; whether to use the C++ implementation of the Gibbs
#'   sampler. If `FALSE`, the R implementation is used. The samplers are equivalent, but the C++ implementation is faster.
#' @param seed_chains Integer; A vector of length `length(chains)` with seeds to which the MCMC chains will be initialized using `set.seed`. If `NULL`, random seeds are generated.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`par.X`}{An `mcmc.list` containing the stored sampled draws of the
#'   `X` transition-model parameters for all chains. The final column is returned
#'   on the natural scale as `sigma`.}
#'   \item{`par.S`}{An `mcmc.list` containing the stored sampled draws of the
#'   `S` transition-model parameters for all chains. The final column is returned
#'   on the natural scale as `sigma`.}
#'   \item{`X`}{Numeric vector containing sampled latent `X` values, concatenated
#'   over chains.}
#'   \item{`S`}{Numeric vector containing sampled latent `S` values, concatenated
#'   over chains.}
#'   \item{`ac`}{Acceptance indicators for the Metropolis-Hastings parameter
#'   updates.}
#'   \item{`ac.cur`}{Acceptance indicators from the current extension run only,
#'   returned when `prev.run` is supplied.}
#'   \item{`dat`}{A `data.frame` containing the data used in fitting, with
#'   columns `d`, `L`, and `R` on the original time scale.}
#'   \item{`priors`}{A list containing the prior hyperparameters
#'   `beta.prior.X`, `beta.prior.S`, `sig.prior.X`, and `sig.prior.S`.}
#'   \item{`thinning`}{The thinning interval used to construct `par.X` and
#'   `par.S`. For long MCMC chains, larger `thinning` values reduce memory
#'   burden.}
#'   \item{`Z.S`}{The covariate matrix used for the `S` transition model.}
#'   \item{`Z.X`}{The covariate matrix used for the `X` transition model.}
#'   \item{`prop.sd`}{The proposal tuning specification used for the
#'   Metropolis-Hastings updates.}
#'   \item{`dist.X`}{The distribution used for the `X` transition model.}
#'   \item{`dist.S`}{The distribution used for the `S` transition model.}
#'   \item{`fix.sigma.X`}{Logical; whether `sigma_X` was fixed during fitting.}
#'   \item{`fix.sigma.S`}{Logical; whether `sigma_S` was fixed during fitting.}
#'   \item{`warmup`}{The number of original Gibbs iterations omitted as warmup
#'   before convergence assessment.}
#'   \item{`beta.prior`}{The prior family used by the default prior function,
#'   either `"t"` or `"norm"`.}
#'   \item{`seed_chains`}{Integer vector of random seeds used for the chains.}
#'   \item{`runtime`}{Total runtime accumulated over the initial run and any
#'   extension runs.}
#'   \item{`update_till_convergence`}{Logical; whether automatic convergence
#'   updating was requested.}
#'   \item{`min_R`}{The convergence threshold used for the potential scale
#'   reduction factor.}
#'   \item{`min_eff`}{The minimum effective sample size threshold used for
#'   convergence assessment.}
#'   \item{`silent`}{Logical; whether progress messages were suppressed.}
#'   \item{`do_cpp`}{Logical; whether the C++ sampler was used.}
#'   \item{`MH`}{Logical; whether Metropolis-Hastings updates were used.}
#'   \item{`slicesampler_stepsize`}{Step-out size used by the slice sampler.}
#'   \item{`log_prior_fun`}{Prior function used for the AFT transition-model
#'   parameters.}
#'   \item{`rescale_times`}{Logical; whether to rescale input times by the median of finite `R`.}
#'   \item{`n_update`}{Integer; counter of the number of times the model was updated.}
#'   \item{`seed_chains`}{Integer; seeds set at initialization of the Gibbs sampler.}
#' }
#'
#' @details
#' The function fits a Bayesian accelerated failure-time model for progressive
#' three-state screening outcomes with censoring after intervention. The time
#' from state 1 to state 2 is denoted `X`, and the time from state 2 to state 3
#' is denoted `S`. In each Gibbs iteration, the latent transition times are
#' updated by data augmentation conditional on the current parameter values, and
#' the transition-model parameters are updated using either slice-sampling or
#' random-walk Metropolis-Hastings steps.
#'
#' The observed data are supplied through the event indicator `d` and interval
#' bounds `L` and `R`, with optional covariate matrices `Z.X` and `Z.S`. The
#' distributions of the latent times `X` and `S` are specified by `dist.X` and
#' `dist.S`.
#'
#' By default, `L` and `R` are divided by the median of the finite right interval
#' bounds before sampling. The fitted object stores the observed data on the
#' original time scale. Klausch et al. (2023) did not rescale time. However, in
#' an AFT model the intercept is measured on the log-time scale, so a fixed
#' default prior on the intercept can have different shrinkage behavior under
#' arbitrary choices of the time unit. Rescaling by a typical observed time scale
#' makes the default prior less sensitive to the units in which time is supplied.
#' This behavior can be disabled with `rescale_times = FALSE`.
#'
#' The default prior is evaluated by [log_aft_prior()]. For each transition model,
#' the AFT parameter vector contains the intercept first, followed by the slope
#' coefficients, with `log(sigma)` in the final position. With
#' `log_prior_fun = log_aft_prior` and `beta.prior = "t"`, independent
#' Student-\eqn{t} priors are assigned to the regression coefficients. The
#' degrees of freedom are given by `beta.prior.X` for the `X` model and
#' `beta.prior.S` for the `S` model. With `beta.prior = "norm"`, independent
#' normal priors are assigned to the regression coefficients, and
#' `beta.prior.X` and `beta.prior.S` are interpreted as prior standard
#' deviations.
#'
#' The scale parameters `sigma_X` and `sigma_S` can either be estimated or fixed.
#' When `fix.sigma.X = FALSE` and `fix.sigma.S = FALSE`, the respective scale
#' parameters are estimated. Under the default prior function, normal prior
#' kernels restricted to positive values are used for `sigma_X` and `sigma_S`,
#' with prior standard deviations `sig.prior.X` and `sig.prior.S`, respectively.
#' Equivalently, these are half-normal priors up to an additive constant in the
#' log-prior. The Jacobian adjustment for the `log(sigma)` parameterization is
#' included by the prior function. When `fix.sigma.X = TRUE` or
#' `fix.sigma.S = TRUE`, the corresponding scale parameter is fixed throughout
#' sampling, with `sig.prior.X` and/or `sig.prior.S` giving the fixed value. For
#' example, `fix.sigma.S = TRUE` with `sig.prior.S = 1` fixes `sigma_S` at 1
#' during estimation. If, in addition, `dist.S = "weibull"` is chosen, the
#' latent `S` transition follows the exponential special case of the Weibull
#' model.
#'
#' Custom prior functions can be supplied through `log_prior_fun`. The same prior
#' function is currently used for the `X` and `S` transition models, but the
#' prior parameters may differ between the two models. Specifically,
#' `beta.prior.X` and `sig.prior.X` are passed as `tau` and `sig.prior` to the
#' prior function for the `X` model, while `beta.prior.S` and `sig.prior.S` are
#' passed as `tau` and `sig.prior` to the prior function for the `S` model.
#' Additional prior parameters, if needed, should currently be hardcoded inside
#' the custom prior function. See [log_aft_prior()] for details on the required
#' parameter vector and on how to modify the default prior.
#'
#' Klausch et al. (2023) used random-walk Metropolis-Hastings updates for the
#' transition-model parameters. This can be requested with `MH = TRUE`. The
#' current default is `MH = FALSE`, which uses a slice sampler instead. In
#' practice, the slice sampler has shown faster convergence. The default
#' step-out size `slicesampler_stepsize = 1` usually works well, but it can be
#' increased or decreased if the slice sampler needs further tuning.
#'
#' Passing a previous fitted object through `prev.run` continues an earlier MCMC
#' run rather than starting from scratch. In that case, the function initializes
#' from the last sampled parameter values of the previous run and reuses the
#' stored data, covariates, priors, proposal settings, distributions, fixed-scale
#' settings, prior function, sampler choice, and controller settings.
#'
#' If `update_till_convergence = TRUE`, the function checks convergence after
#' the initial run and repeatedly extends the sampler by `mc_update` iterations
#' per chain until the convergence criteria based on `min_R` and `min_eff` are
#' satisfied, or until `maxit` is reached. Progress is printed after each
#' extension run of `mc_update` iterations per chain.
#'
#' @seealso [log_aft_prior()] [ppCIF()] [search_prop_sd()] [trim.mcmc()]
#'
#' @references
#' Klausch, T., Akwiwu, E. M. U., van de Wiel, M. A., Coupé, V. M. H., and
#' Berkhof, J. (2023). A Bayesian accelerated failure time model for interval
#' censored three-state screening outcomes. \emph{The Annals of Applied
#' Statistics}, 17(2), 1285--1306. \doi{10.1214/22-AOAS1669}
#'
#' @examples
#' library(BayesTSM)
#'
#' # Generate data
#' dat = gendat(
#'               n = 1000,  # Sample size
#'               p = 2,     # Number of normal distributed covariates
#'               r = 0,     # Correlation of the covariates
#'               sigma.X = 0.3,         # True scale parameter
#'               mu.X    = 2,           # True intercept parameter
#'               beta.X  = c(0.5,0.5),  # True slope parameters
#'               sigma.S = 0.5,         # True scale parameter
#'               mu.S    = 1,           # True intercept parameters
#'               beta.S  = c(0.5,0.5),  # True slope parameters
#'               dist.X  = 'weibull',   # Distribution of X
#'               dist.S  = 'weibull',   # Distribution of S
#'               v.min   = 1,           # Min time between screening moments
#'               v.max   = 5,           # Max time between screening moments
#'               Tmax    = 2e2,         # Max number of screening times
#'               mean.rc = 10           # Mean time to right censoring
#'             )
#'
#' # Run bayestsm Gibbs sampler with data augmentation and slice sampling of the parameters
#' mod_slice   = bayestsm(
#'              d              = dat$d,
#'              L              = dat$L,
#'              R              = dat$R,
#'              Z.X            = dat[,c('Z.1','Z.2')],
#'              Z.S            = dat[,c('Z.1','Z.2')],
#'              mc             = 1e4,
#'              warmup         = 5e2,              # discarded before assessing convergence
#'              thinning        = 10,               # The chain can be thinned to save memory
#'              chains         = 2,                # In practice use more than 2 chains
#'              update_till_convergence = FALSE,
#'              MH             = FALSE,            # set to TRUE for Metropolis sampling
#'              dist.X         = 'weibull',        # Correctly specified distributions
#'              dist.S         = 'weibull'
#'              )
#'
#' # Run bayestsm Gibbs sampler with data augmentation and Metropolis sampling of the parameters
#' mod_MH      = bayestsm(
#'              d              = dat$d,
#'              L              = dat$L,
#'              R              = dat$R,
#'              Z.X            = dat[,c('Z.1','Z.2')],
#'              Z.S            = dat[,c('Z.1','Z.2')],
#'              mc             = 1e4,
#'              warmup         = 5e2,              # discarded before assessing convergence
#'              thinning        = 10,               # The chain can be thinned to save memory
#'              chains         = 2,                # In practice use more than 2 chains
#'              update_till_convergence = FALSE,
#'              MH             = TRUE,             # set to TRUE for Metropolis sampling
#'              dist.X         = 'weibull',        # Correctly specified distributions
#'              dist.S         = 'weibull'
#'              )
#'
#' @export
bayestsm = function(
    d = NULL,
    L = NULL,
    R = NULL,
    Z.X = NULL,
    Z.S = Z.X,

    dist.X = "weibull",
    dist.S = "weibull",

    beta.prior = "t",
    beta.prior.X = 4,
    beta.prior.S = beta.prior.X,
    sig.prior.X = 1,
    sig.prior.S = 1,
    fix.sigma.X = FALSE,
    fix.sigma.S = FALSE,

    mc = 1e4,
    chains = 2,
    thinning = 1,
    warmup = mc * 0.05,
    prop.sd = NULL,

    update_till_convergence = FALSE,
    mc_update = mc,
    min_R = 1.01,
    min_eff = chains * 100,
    maxit = NULL,

    prev.run = NULL,

    MH     = FALSE,
    rescale_times = TRUE,
    slicesampler_stepsize = 1,
    silent = FALSE,
    log_prior_fun = log_aft_prior,
    do_cpp = TRUE,
    seed_chains = NULL
){

  total_runtime = 0
  # Load parameters from previous run
  if(!is.null(prev.run)){
    chains = length(prev.run$par.X)
    dims.X <- dim(as.matrix(prev.run$par.X[[1]]))
    dims.S <- dim(as.matrix(prev.run$par.S[[1]]))

    start.val.X <- matrix(ncol = dims.X[2], nrow = chains)
    start.val.S <- matrix(ncol = dims.S[2], nrow = chains)

    for(i in 1:chains){
      start.val.X[i,] = as.matrix(prev.run$par.X[i])[dims.X[1],]
      start.val.S[i,] = as.matrix(prev.run$par.S[i])[dims.S[1],]
    }
    d = prev.run$dat$d
    L = prev.run$dat$L
    R = prev.run$dat$R
    Z.X = prev.run$Z.X
    Z.S = prev.run$Z.S
    if(is.null(prop.sd)){ prop.sd =  prev.run$prop.sd }
    beta.prior.X = prev.run$priors$beta.prior.X
    beta.prior.S = prev.run$priors$beta.prior.S
    sig.prior.X = prev.run$priors$sig.prior.X
    sig.prior.S = prev.run$priors$sig.prior.S
    thinning = prev.run$thinning
    X.prev = prev.run$X
    S.prev = prev.run$S
    dist.X = prev.run$dist.X
    dist.S = prev.run$dist.S
    fix.sigma.X = prev.run$fix.sigma.X
    fix.sigma.S = prev.run$fix.sigma.S
    beta.prior = prev.run$beta.prior
    mc = mc_update
    update_till_convergence = prev.run$update_till_convergence
    min_R = prev.run$min_R
    min_eff = prev.run$min_eff
    total_runtime = prev.run$runtime
    silent = prev.run$silent
    do_cpp = prev.run$do_cpp
    n_update = prev.run$n_update
    seed_chains = prev.run$seed_chains

    MH = prev.run$MH
    rescale_times = prev.run$rescale_times
    slicesampler_stepsize = prev.run$slicesampler_stepsize
    log_prior_fun = prev.run$log_prior_fun

    if(!update_till_convergence & !silent) message('Updating previous MCMC run.') #Overriding data, chains, thinning, burn-in, and priors
  }
  timestamp1 <- Sys.time()
  # Make basic data checks and pre-process
  inp <- check_bayestsm_inputs(
    d = d,
    L = L,
    R = R,
    Z.X = Z.X,
    Z.S = Z.S,
    w = 1,
    prop.sd = prop.sd,
    prev.run = prev.run
  )

  d = inp$d
  L = inp$L
  R = inp$R
  Z.X = inp$Z.X
  Z.S = inp$Z.S
  Z1.X = inp$Z1.X
  Z1.S = inp$Z1.S
  p.X = inp$p.X
  p.S = inp$p.S
  p1.X = inp$p1.X
  p1.S = inp$p1.S
  prop.sd = inp$prop.sd

  # Definitions and scaling
  n   = length(d)
  d23 = d == 2 | d == 3
  n23 = sum(d23)
  if(rescale_times) medLR = median(c(L[is.infinite(R)], R[is.finite(R)])) else medLR = 1
  Lorig = L
  Rorig = R
  L   = L/medLR
  R   = R/medLR

  # Determine useful prod.sd if not provided
  if(is.null(prop.sd) && MH){
    if(!silent) message(sprintf('No proposal sd provided. Searching.'))
    mod_ini   = bayestsm( d              = d,
                          L              = L,
                          R              = R,
                          Z.X            = Z.X,
                          Z.S            = Z.S,
                          mc             = 3e3,
                          chains         = chains,
                          prop.sd        = 0.01,
                          MH             = TRUE,
                          silent = T,
                          seed_chains = seed_chains)

    s = search_prop_sd(mod_ini, mc = 5e3)
    prop.sd = s$prop.sd
  }
  if(!MH) prop.sd = 1

  # Get random seeds
  if(is.null(seed_chains)) seed_chains = sample(1:10^5,chains, replace=F)

  ## Start run
  if(is.null(prev.run) | (!is.null(prev.run) & !update_till_convergence)){
    if(!silent) message(sprintf("Starting Gibbs sampler with %d chains and %d iterations.", chains, mc))
  }
  cl    = makePSOCKcluster(chains)
  clusterSetRNGStream(cl)
  registerDoParallel(cl)

    run <- foreach::foreach(
    chain_id = seq_len(chains),
    .packages = c("BayesTSM", "MCMCpack", "mvtnorm")
    ) %dopar% {

    # Set seed of chain
      if (is.null(prev.run)) {
        set.seed(seed_chains[chain_id])
      } else {
        .Random.seed <- prev.run$rng_state[[chain_id]]
      }

    # Search starting values
    if(is.null(prev.run)){
      ini <- find_starting_values(
        L = L,
        R = R,
        d = d,
        Z.X = Z.X,
        Z.S = Z.S,
        Z1.X = Z1.X,
        Z1.S = Z1.S,
        p.X = p.X,
        p.S = p.S,
        dist.X = dist.X,
        dist.S = dist.S,
        beta.prior.X = beta.prior.X,
        beta.prior.S = beta.prior.S,
        sig.prior.X = sig.prior.X,
        sig.prior.S = sig.prior.S,
        fix.sigma.X = fix.sigma.X,
        fix.sigma.S = fix.sigma.S,
        beta.prior = beta.prior,
        log_prior_fun = log_prior_fun
      )
      x.ini <- ini$x.ini
      s.ini <- ini$s.ini
      beta.x.ini <- ini$beta.x.ini
      sigma.x.ini <- ini$sigma.x.ini
      beta.s.ini <- ini$beta.s.ini
      sigma.s.ini <- ini$sigma.s.ini
      par.x <- ini$par.x
      par.s <- ini$par.s
      log.pst.X <- ini$log.pst.X
      log.pst.S <- ini$log.pst.S
      counter <- ini$counter
      #list2env(ini, envir = environment())
    }

    if(!is.null(prev.run)){
      if (p.X > 0) beta.x.ini <- start.val.X[chain_id, 2:(p.X + 1)] else beta.x.ini <- NULL
      if (p.S > 0) beta.s.ini <- start.val.S[chain_id, 2:(p.S + 1)] else beta.s.ini <- NULL
      sigma.x.ini = log(start.val.X[chain_id,(p.X+2)])
      sigma.s.ini = log(start.val.S[chain_id,(p.S+2)])
      mu.x.ini = start.val.X[chain_id,1] - log(medLR) #rescaling
      mu.s.ini = start.val.S[chain_id,1] - log(medLR) #rescaling
      x.ini = X.prev[ ((n*(chain_id-1))+1): (n*chain_id) ]
      s.ini = S.prev[ ((n*(chain_id-1))+1): (n*chain_id) ]
      beta.x.ini = c(mu.x.ini, beta.x.ini)
      beta.s.ini = c(mu.s.ini, beta.s.ini)
    }

    # Run Gibbs chain conditional on starting values
    if(do_cpp){
    gibbs_out <- full_conditional_Gibbs_cpp(
      mc = mc,
      d = d,
      d23 = d23,
      L = L,
      R = R,
      Z1.X = Z1.X,
      Z1.S = Z1.S,
      p.X = p.X,
      p.S = p.S,
      p1.X = p1.X,
      p1.S = p1.S,
      x.ini = x.ini,
      s.ini = s.ini,
      beta.x.ini = beta.x.ini,
      sigma.x.ini = sigma.x.ini,
      beta.s.ini = beta.s.ini,
      sigma.s.ini = sigma.s.ini,
      prop.sd = prop.sd,
      dist.X = dist.X,
      dist.S = dist.S,
      beta.prior.X = beta.prior.X,
      beta.prior.S = beta.prior.S,
      sig.prior.X = sig.prior.X,
      sig.prior.S = sig.prior.S,
      fix.sigma.X = fix.sigma.X,
      fix.sigma.S = fix.sigma.S,
      beta.prior = beta.prior,
      MH = MH,
      slicesampler_stepsize = slicesampler_stepsize,
      log_prior_fun = log_prior_fun
    )}
    else{
    gibbs_out <- full_conditional_Gibbs(
      mc = mc,
      d = d,
      d23 = d23,
      L = L,
      R = R,
      Z1.X = Z1.X,
      Z1.S = Z1.S,
      p.X = p.X,
      p.S = p.S,
      p1.X = p1.X,
      p1.S = p1.S,
      x.ini = x.ini,
      s.ini = s.ini,
      beta.x.ini = beta.x.ini,
      sigma.x.ini = sigma.x.ini,
      beta.s.ini = beta.s.ini,
      sigma.s.ini = sigma.s.ini,
      prop.sd = prop.sd,
      dist.X = dist.X,
      dist.S = dist.S,
      beta.prior.X = beta.prior.X,
      beta.prior.S = beta.prior.S,
      sig.prior.X = sig.prior.X,
      sig.prior.S = sig.prior.S,
      fix.sigma.X = fix.sigma.X,
      fix.sigma.S = fix.sigma.S,
      beta.prior = beta.prior,
      MH = MH,
      slicesampler_stepsize = slicesampler_stepsize,
      log_prior_fun = log_prior_fun
    )}

    out=list()
    out$X     = gibbs_out$X
    out$S     = gibbs_out$S
    out$par.X = gibbs_out$cur.par.Xreg
    out$par.S = gibbs_out$cur.par.Sreg
    out$ac    = gibbs_out$ac
    out$rng_state = .Random.seed
    out
    }
    stopCluster(cl)

  # Householding
  rng_state <- lapply(run, `[[`, "rng_state")
  assembled <- assemble_bayestsm_mcmc(
    run = run,
    Z.X = Z.X,
    Z.S = Z.S,
    burnin = 1,
    thinning = thinning,
    mc = mc,
    prev.run = prev.run,
    update.burnin = F,
    rescale_intercept = medLR
  )

  par.X  <- assembled$par.X.all
  par.S  <- assembled$par.S.all
  X.draw    <- assembled$X.draw
  S.draw    <- assembled$S.draw
  ac        <- assembled$ac

  if (!is.null(prev.run)) {
    ac.cur <- assembled$ac.cur
  }

  # save prior info
  dat = data.frame(d=d, L=Lorig, R=Rorig)
  priors = list()
  priors$beta.prior.X = beta.prior.X
  priors$beta.prior.S = beta.prior.S
  priors$sig.prior.X = sig.prior.X
  priors$sig.prior.S = sig.prior.S

  # Exporting
  timestamp2 <- Sys.time()
  ret = list()
  ret$par.X = par.X
  ret$par.S = par.S
  ret$X = X.draw
  ret$S = S.draw
  ret$ac = ac
  if(!is.null(prev.run)){
    ret$ac.cur = ac.cur}
  ret$dat = dat
  ret$priors = priors
  ret$thinning = thinning
  ret$Z.S = Z.S
  ret$Z.X = Z.X
  ret$prop.sd = prop.sd
  ret$dist.X = dist.X
  ret$dist.S = dist.S
  ret$fix.sigma.X = fix.sigma.X
  ret$fix.sigma.S = fix.sigma.S
  ret$warmup = warmup
  ret$beta.prior = beta.prior
  ret$seed_chains = seed_chains
  ret$runtime = total_runtime + (timestamp2 - timestamp1)
  ret$update_till_convergence = update_till_convergence
  ret$min_R = min_R
  ret$min_eff = min_eff
  ret$silent = silent
  ret$do_cpp = do_cpp
  ret$MH = MH
  ret$slicesampler_stepsize = slicesampler_stepsize
  ret$log_prior_fun = log_prior_fun
  ret$rescale_times = rescale_times
  if(is.null(prev.run)) n_update = 0
  ret$n_update      = n_update + 1
  ret$seed_chains   = seed_chains
  ret$rng_state     = rng_state

  # Convergence checks and model updating
  conv <- handle_bayestsm_convergence(
    ret = ret,
    update_till_convergence = update_till_convergence,
    chains = chains,
    min_R = min_R,
    min_eff = min_eff,
    maxit = if(!is.null(maxit)) round(maxit/thinning) else NULL,
    burnin = round(warmup / thinning),
    mc_update = mc_update,
    silent = silent,
    fix.sigma.X = fix.sigma.X,
    fix.sigma.S = fix.sigma.S
  )

  ret <- conv$ret

  if (!update_till_convergence) {
    conv <- handle_bayestsm_convergence(
      ret = ret,
      update_till_convergence = T,
      chains = chains,
      min_R = min_R,
      min_eff = min_eff,
      maxit = if(!is.null(maxit)) round(maxit/thinning) else NULL,
      burnin = round(warmup / thinning),
      mc_update = mc_update,
      silent = silent,
      fix.sigma.X = fix.sigma.X,
      fix.sigma.S = fix.sigma.S
    )

    ret <- conv$ret
    ret$call <- match.call()
    class(ret) <- c("bayestsm", class(ret))
    return(ret)
  }

  while (isTRUE(conv$needs_update)) {
    n_iter_before <- nrow(ret$par.X[[1L]])

    # Important:
    # This object is only used for a single extension run.
    # The automatic convergence controller must not recurse.
    # It is also made silent so that only the controller prints progress messages.
    prev_for_update <- ret
    prev_for_update$update_till_convergence <- FALSE
    prev_for_update$silent <- TRUE

    ret <- bayestsm(
      prev.run = prev_for_update,
      warmup = prev_for_update$warmup,
      mc_update = conv$mc_update_next,
      min_R = min_R,
      min_eff = min_eff,
      maxit = if(!is.null(maxit)) round(maxit/thinning) else NULL,
      update_till_convergence = FALSE,
      do_cpp = do_cpp,
      silent = TRUE,
    )

    # Restore external-facing settings.
    ret$update_till_convergence <- TRUE
    ret$min_R <- min_R
    ret$min_eff <- min_eff
    ret$silent <- silent
    ret$warmup <- warmup

    rm(prev_for_update)
    invisible(gc(verbose = FALSE))

    n_iter_after <- nrow(ret$par.X[[1L]])

    if (n_iter_after <= n_iter_before) {
      stop(
        paste0(
          "The convergence update did not add any stored MCMC iterations. ",
          "Increase mc_update or decrease thinning."
        ),
        call. = FALSE
      )
    }

    conv <- handle_bayestsm_convergence(
      ret = ret,
      update_till_convergence = TRUE,
      chains = chains,
      min_R = min_R,
      min_eff = min_eff,
      maxit = if(!is.null(maxit)) round(maxit/thinning) else NULL,
      burnin = round(warmup / thinning),
      mc_update = mc_update,
      silent = silent,
      fix.sigma.X = fix.sigma.X,
      fix.sigma.S = fix.sigma.S
    )

    ret <- conv$ret
  }

  ret$call <- match.call()
  class(ret) <- c("bayestsm", class(ret))
  return(ret)
}
