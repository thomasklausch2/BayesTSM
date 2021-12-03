### Demonstration of bts_survreg - Bayesian 3 state survival regression
# This vignette uses simulated data to demonstrate use of the bts_survreg function
library(doParallel)
library(MCMCpack)
library(coda)
library(mvtnorm)
library(MHadaptive)
library(mcmcr)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(actuar)

# We source startup.r which loads all required functions
source('controllers/startup.r')

# Function gendat generates 3-state interval censored time to event (screening) data
# The parameters of the interval censoring mechanism are described in the main paper, section 5, simulation and in the sipplemental material
set.seed(022021)
dat = gendat(n = 2000,  # Sample size
             p = 2,     # Number of normal distributed covariates
             r = 0,     # Correlation of the covariates
             sigma.X = 0.3,           # True scale parameter
             mu.X    = 2,           # True intercept parameter
             beta.X  = c(0.5,0.5),  # True slope parameters of all covariates of X
             sigma.S = 0.2,         # True scale parameter
             mu.S    = 1,         # True intercept parameters
             beta.S  = c(0.5,0.5),  # True slope parameter of all covariates of S
             dist.X  = 'weibull',   # Distribution of X, alternatives are loglog and lognormal
             dist.S  = 'weibull',   # Distribution of S, alternatives are loglog and lognormal
             v.min   = 1,           # Minimum time between screening moments 
             v.max   = 5,           # Maximum time between screening moments
             Tmax    = 2e2,         # Maximum number of screening times (this is set to high value; but too high values increase computation time)
             mean.rc = 10           # Mean time to right censoring (parameter of exponential distribution)
             )


# screen the simulated data
head(dat) # Contains also true (latent) X and S (not observed in reality)
#summary(dat)
prop.table(table(dat$d)) # Distribution of censoring events

# Function bts_survreg is used to estimate the model (parallel processing requires at least 3 free CPUs)
mod         = bts_survreg(d              = dat$d, # censoring indicator, assumes values 1, 2, 3
                          L              = dat$L, # Time of left censoring; assumes time of last visit of d=1
                          R              = dat$R, # Time of right censoring; assumes inf if d=1
                          Z.X            = dat[,c('Z.1','Z.2')], # Covariates of X
                          Z.S            = dat[,c('Z.1','Z.2')], # Covariates of S
                          mc             = 2e4,       # MCMC draws (can be updated later, see below), half will be dropped for burn-in (other values can also be specified using brunin argument)
                          chains         = 3,         # Number of parallel MCMC chains
                          thin           = 5,         # The function returns thinned and unthinned chains; there the thinning interval may be given
                          do.seperate.MH = F,         # Whether the metropolis step should be done jointly (F) or seperately for parameters of X and S
                          prop.sd.X      = 0.006,      # The proposal standard deviation of the normal distribution used for the metropolis step (if do.seperate.MH = T also prop.sd.S should be tuned)
                          beta.prior.X   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (X)
                          beta.prior.S   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (S)
                          sig.prior.X    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (X)
                          sig.prior.S    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (S)
                          dist.X         = 'weibull', # Distribution of X
                          dist.S         = 'weibull', # Distribution of S
                          fix.sigma.X    = F,         # Should sigma of X be fixed at its prior value?
                          fix.sigma.S    = F,         # Should sigma of S be fixed at its prior value? (use e.g. to specify exponential distribution)
                          parallel       = T,         # TRUE if chains should be run on seperate CPUs
                          beta.prior     = 't')       # Alternatively a normal prior for beta can be specified with 'norm'

# Object returns the following values in a list
names(mod)
head(mod$par.X.all[1]) # mcmc.list for all draws of parameters of X
head(mod$par.S.all[1]) # mcmc.list for all draws of parameters of S
head(mod$par.X.bi[1])  # mcmc.list for all draws of parameters of X after burn-in and thinning
head(mod$par.S.bi[1])  # mcmc.list for all draws of parameters of S after burn-in and thinning
head(mod$X)            # Last sample of X
head(mod$S)            # Last sample of S
mod$ac.X               # Matrix with ncol = chains and nrow = mc giving the metropolis acceptance (1 if accepeted)
mod$ac.S               # Matrix with ncol = chains and nrow = mc giving the metropolis acceptance for parameters of S (only if do.seperate.MH==T)
head(mod$dat)          # Copy of data
mod$priors             # Named list with hyperprior values
mod$thining            # The thinning interval chosen
head(mod$Z.X)          # Copy of the covariates of X
head(mod$Z.S)          # Copy of the covariates of S
mod$prop.sd.X          # Copy of prop.sd.X
mod$prop.sd.S          # Copy of prop.sd.S
mod$dist.X             # Copy of dist.X
mod$dist.S             # Copy of dist.S
mod$fix.sigma.X        # Copy of fix.sigma.X
mod$fix.sigma.S        # Copy of fix.sigma.S
mod$burnin             # Copy of burning
mod$beta.prior         # Copy of beta.prior

# Some basic posterior statistics
apply(mod$ac.X,2,mean) # Metropolis acceptance rate by chain
summary(mod$par.X.bi)  # Standard MCMC summaries of mcmc.list using MCMCpack
summary(mod$par.S.bi)  # Standard MCMC summaries of mcmc.list using MCMCpack
plot(mod$par.X.bi, smooth=F)  # Standard MCMC chain plots of mcmc.list using MCMCpack
plot(mod$par.S.bi, smooth=F)  # Standard MCMC chain plots of mcmc.list using MCMCpack
gelman.diag(mod$par.X.bi)     # Gelman convergence diagnostics
gelman.diag(mod$par.S.bi)     # Gelman convergence diagnostics

# Updating MCMC chain (this can take a bit longer)
mod.updated = bts_survreg(prev.run = mod, mc = 1e5) # An update for another 10^5 mc draws
gelman.diag(mod.updated$par.X.bi)     # Gelman convergence diagnostics
gelman.diag(mod.updated$par.S.bi)     # Gelman convergence diagnostics
plot(mod.updated$par.X.bi, smooth=F)  # Standard MCMC chain plots of mcmc.list using MCMCpack
plot(mod.updated$par.S.bi, smooth=F)  # Standard MCMC chain plots of mcmc.list using MCMCpack

# The proposal SD of the Metropolis samples has to be tuned manually
# Alternatively we use a heuristic search in the main paper, described in detail in the supplemental material
# First, a short run of the main model is required with a default starting rate (here .1 which is too high)
mod.ini     = bts_survreg(d              = dat$d, # censoring indicator, assumes values 1, 2, 3
                          L              = dat$L, # Time of left censoring; assumes time of last visit of d=1
                          R              = dat$R, # Time of right censoring; assumes inf if d=1
                          Z.X            = dat[,c('Z.1','Z.2')], # Covariates of X
                          Z.S            = dat[,c('Z.1','Z.2')], # Covariates of S
                          mc             = 5e3,       # MCMC draws (can be updated later, see below), half will be dropped for burn-in (other values can also be specified using brunin argument)
                          chains         = 3,         # Number of parallel MCMC chains
                          thin           = 5,         # The function returns thinned and unthinned chains; there the thinning interval may be given
                          do.seperate.MH = F,         # Whether the metropolis step should be done jointly (F) or seperately for parameters of X and S
                          prop.sd.X      = .1,        # The proposal standard deviation of the normal distribution used for the metropolis step (if do.seperate.MH = T also prop.sd.S should be tuned)
                          beta.prior.X   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (X)
                          beta.prior.S   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (S)
                          sig.prior.X    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (X)
                          sig.prior.S    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (S)
                          dist.X         = 'weibull', # Distribution of X
                          dist.S         = 'weibull', # Distribution of S
                          fix.sigma.X    = F,         # Should sigma of X be fixed at its prior value?
                          fix.sigma.S    = F,         # Should sigma of S be fixed at its prior value? (use e.g. to specify exponential distribution)
                          parallel       = T,         # TRUE if chains should be run on seperate CPUs
                          beta.prior     = 't')       # Alternatively a normal prior for beta can be specified with 'norm'

# Second, we run search.prop.sd implementing the heuristic search
s   = search.prop.sd(mod.ini, acc.bounds.X =c(0.21,0.25)) # acc.bounds.X gives acceptable acceptance rate bounds as described in the supplemental material
s$prop.sd.X # The tuned acceptance probability


# Once the model is fit we wish to look at the predictive cumulative incidence functions
# Marginal predictive CIF can be obtained via the quantile function for given percentiles perc
cif.q = get.pCIF.q(mod = mod,           # The fitted model of bts_survreg
                pst.samples = 10^3,     # The number of posterior draws used in calculating the predictive CIF (computation time and accuracy scale with this number)
                perc = seq(0, 1, 0.01)  # Vector of percentiles / probabilities at which we want to evaluate the inverse of the CDF (length scales computation time)
              )
# Or for given quantiles
cif.p = get.pCIF.p(mod = mod,           # The fitted model of bts_survreg
                pst.samples = 10^3,     # The number of posterior draws used in calculating the predictive CIF (computation time and accuracy scale with this number)
                q    = seq(0, 20, 0.01) # Vector of quantiles ('x-values') at which we want to evaluate the CDF (length scales computation time)
)

# Plots of obtained quantiles against specified percentiles
par(mfrow=c(1,3))
plot(cif.q$med.cdf.x, cif.q$perc, ty = 'l', main = 'Time X', xlim=c(0,20))
lines(ecdf(dat$X),col=2) # True empirical CDF (data)
lines( cif.q$med.cdf.x.ci[1,], cif.q$perc  ,col=3, lty=2)
lines( cif.q$med.cdf.x.ci[2,], cif.q$perc  ,col=3, lty=2)

plot(cif.q$med.cdf.s, cif.q$perc, ty = 'l', main = 'Time S', xlim=c(0,20))
lines(ecdf(dat$S),col=2) # True empirical CDF (data)
lines( cif.q$med.cdf.s.ci[1,], cif.q$perc  ,col=3, lty=2)
lines( cif.q$med.cdf.s.ci[2,], cif.q$perc  ,col=3, lty=2)

plot(cif.q$med.cdf.xs, cif.q$perc, ty = 'l', main = 'Time X+S', xlim=c(0,20))
lines(ecdf(dat$X+dat$S),col=2) # True empirical CDF (data)
lines( cif.q$med.cdf.xs.ci[1,], cif.q$perc  ,col=3, lty=2)
lines( cif.q$med.cdf.xs.ci[2,], cif.q$perc  ,col=3, lty=2)

# Plots of obtained percentiles against specified quantiles (approx. equal to previous plot)
par(mfrow=c(1,3))
plot(cif.p$q, cif.p$med.cdf.x, ty = 'l', main = 'Conditional CIF of X', ylim=c(0,1))
lines(ecdf(dat$X),col=2) # True empirical CDF (data)
lines(cif.p$q, cif.p$med.cdf.x.ci[1,], col=3, lty=2)
lines(cif.p$q, cif.p$med.cdf.x.ci[2,], col=3, lty=2)

plot(cif.p$q, cif.p$med.cdf.s, ty = 'l', main = 'Conditional CIF of S', ylim=c(0,1))
lines(ecdf(dat$S),col=2) # True empirical CDF (data)
lines(cif.p$q, cif.p$med.cdf.s.ci[1,], col=3, lty=2)
lines(cif.p$q, cif.p$med.cdf.s.ci[2,], col=3, lty=2)

plot(cif.p$q, cif.p$med.cdf.xs, ty = 'l', main = 'Conditional CIF of X+S', ylim=c(0,1))
lines(ecdf(dat$X+dat$S),col=2) # True empirical CDF (data)
lines(cif.p$q, cif.p$med.cdf.xs.ci[1,], col=3, lty=2)
lines(cif.p$q, cif.p$med.cdf.xs.ci[2,], col=3, lty=2)

# Predictive CIF for specified quantile (e.g. 10)
get.pCIF.p(mod = mod, pst.samples = 10^3, q = c(0,10))

# Predictive (inverse) CIF for specified percentile (e.g. median)
get.pCIF.q(mod = mod, pst.samples = 10^3, p = c(0,0.5))


# Conditional predictive CIF 

ppd1 = get.ppd.cond( mod, 
                     Z = c(0.5,0.5),         # The values for a new z to condition on 
                     pst.samples = 1000,     # The number of posterior draws used in calculating the ppd (computation time and accuracy scale with this number)
                     grid = seq(0,100,0.5)   # The quantiles grid at which to evaluate the CDF; note that currently, contrary to get.ppd, get.ppd.cond evaluates the CDF at pre-specified quantiles
                    )

par(mfrow=c(1,3))
plot(ppd1$grid, ppd1$med.cdf.x, ty = 'l', main = 'Conditional CIF of X', ylim=c(0,1))
lines(ppd1$grid, ppd1$med.cdf.x.ci[1,], col=3, lty=2)
lines(ppd1$grid, ppd1$med.cdf.x.ci[2,], col=3, lty=2)

plot(ppd1$grid, ppd1$med.cdf.s, ty = 'l', main = 'Conditional CIF of S', ylim=c(0,1))
lines(ppd1$grid, ppd1$med.cdf.s.ci[1,], col=3, lty=2)
lines(ppd1$grid, ppd1$med.cdf.s.ci[2,], col=3, lty=2)

plot(ppd1$grid, ppd1$med.cdf.xs, ty = 'l', main = 'Conditional CIF of X+S', ylim=c(0,1))
lines(ppd1$grid, ppd1$med.cdf.xs.ci[1,], col=3, lty=2)
lines(ppd1$grid, ppd1$med.cdf.xs.ci[2,], col=3, lty=2)


# Also we can get the information criteria
ic.fit = get.IC(mod, 
                samples = 1e3, # The number of posterior draws used in calculating the ppd (computation time and accuracy scale with this number)
                cores = NULL  # Number of cores used in parallel computation, if NULL all cores are used (currently only parallel computation is possible)
                )
ic.fit


# Now a sensitivity analysis for a latent confounder is conducted (note you may need to choose a wider grid in practice)
source('source/sensAnalysis.step2.r')
b.UX.grid = c(-1,1)
b.US.grid = c(1,1)
m.sens = list()
for(i in 1:length(b.UX.grid)){
  m.sens[[i]] = sensAnalysis.step2(mod, b.UX = b.UX.grid[i], b.US = b.US.grid[i], mc  = 1e4, chains = 3, parallel = T)
}




