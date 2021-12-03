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

# We source startup.r which loads all required functions for the Bayesian model
source('controllers/startup.r')
# In addition we manually load the EM add-on
source('source/EM/px.r')
source('source/EM/ps.r')
source('source/EM/truncbounds.s.r')
source('source/EM/truncbounds.x.r')
source('source/EM/Ex.int.r')
source('source/EM/Es.int.r')
source('source/EM/EX.truncnorm.r')
source('source/EM/map.r')
source('source/EM/em.lognorm.r')
source('source/EM/fastlognormalconv.r')
#


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
             dist.X  = 'lognormal',   # Distribution of X, alternatives are loglog and lognormal
             dist.S  = 'lognormal',   # Distribution of S, alternatives are loglog and lognormal
             v.min   = 1,           # Minimum time between screening moments 
             v.max   = 5,           # Maximum time between screening moments
             Tmax    = 2e2,         # Maximum number of screening times (this is set to high value; but too high values increase computation time)
             mean.rc = 10           # Mean time to right censoring (parameter of exponential distribution)
)

# screen the simulated data
head(dat) # Contains also true (latent) X and S (not observed in reality)
#summary(dat)
prop.table(table(dat$d)) # Distribution of censoring events

# Running the EM algorithm
fit <- em.lognorm(d = dat$d, # Screening event
                  L = dat$L, # Left bound
                  R = dat$R, # Right bound
                  Z.X = dat[,c('Z.1','Z.2')], #Covariates X
                  Z.S = dat[,c('Z.1','Z.2')], #Covariates S
                  sig.x = 1, #Starting value sigma of X
                  sig.s = 1, #Starting value sigma of S
                  beta.x = c(0,0,0), #Starting value betas of X (including the intercept)
                  beta.s = c(0,0,0), #Starting value betas of S (including the intercept)
                  n.rej = 1e2,     #Number of importance samples used for the numerical integration in the E-step (may need to be set higher)
                  tol = 1e-4,      #Convergence limit (may need to be set lower)
                  max.it = 1e3,    # Maximum number of iterations
                  silent = F       # Should output be shown while running?
)
  
fit$beta.x #final estimates beta X
fit$beta.s #final estimates beta S
fit$sig.x  #final estimates sigma X
fit$sig.s  #final estimates sigma S
fit$beta.x.seq # whole sequence of beta of X
fit$beta.s.seq # whole sequence of beta of S
fit$sig.x.seq # whole sequence of sigma of X
fit$sig.s.seq # whole sequence of sigma of S
fit$ll        # likelihood of last run
fit$ll.seq    # likelihood sequence of all runs
fit$convergence # Was convergence reaches before max.it?
fit$er          # Difference between ll of last run and second to last run
fit$tol         # Convergence limit

