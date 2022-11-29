
rm(list=ls()) # clear workspace

# Load packages -----------------------------------------------------------

source('fx.R')  # for simulation and DIC
library(mnormt) # for simulations
library(brms)   # for regressions

# Simulate true model -----------------------------------------------------

set.seed(2022)

# Set constants
N = J = c(10,50,100) # cluster and obs per cluster sizes to try
TAU0  = c(0.05,1,5)  # degrees of clustering (holding sigma constant)
S     = 100          # number of iterations per sample per model
nchain = 1           # number of chains per MCMC sample per model
niter = 3000         # number of MCMC samples per model
burn  = 1000         # number of burn-in MCMC samples per model
thin  = 1            # amount of thinning per MCMC sample per model

# specify experimental conditions
conds = expand.grid(N=N, J=J,TAU0=TAU0)

# Set population parameters for two predictor model (same in each sample)
GAMMA = c('(Intercept)' = 1.1, 'X1'=1.66, 'X2'=1.2) # fixed effects
TAU   = c('id.X1' = 0.75, 'id.X1.(Intercept)'=0)    # random variance (excluding TAU0, which changes each sample)
SIGMA = 3 # residual variance

# Fit initial model to dummy data, then just call update below to avoid recompiling time.
# Fit five models:
# A. Excludes random slopes for X    (random restricted)
# B. Excludes X2 as a predictor      (fixed restricted)
# C. Estimate random effects for X2  (overfit random)
# D. Estimate fixed effects for X3   (overfit fixed)
# E. True model

# DM    = expand.grid(1:10,1:10) # design matrix
# DM$X1 = rnorm(nrow(DM))        # standard normal predictor X1
# DM$X2 = rnorm(nrow(DM),2,1)    # predictor X2 mu = 2, s =1
# DM$X3 = rnorm(nrow(DM),0,2)    # non-predictor mu = 0, s =2
# DM$y  = rnorm(nrow(DM))        # standard normal random outcome
# colnames(DM) = c('obs','id','X1','X2','X3', 'y')
# 
# modA = brm(y ~ X1 + X2 +(1|id),       data=DM, chains = 1, iter=niter, warmup=burn, thin=thin, silent = 2)
# modB = brm(y ~ X1 +(X1||id),          data=DM, chains = 1, iter=niter, warmup=burn, thin=thin, silent = 2)
# modC = brm(y ~ X1 + X2 +(X1+X2||id),  data=DM, chains = 1, iter=niter, warmup=burn, thin=thin, silent = 2)
# modD = brm(y ~ X1 + X2 + X3 + (1|id), data=DM, chains = 1, iter=niter, warmup=burn, thin=thin, silent = 2)
# modE = brm(y ~ X1 + X2 +(X1||id),     data=DM, chains = 1, iter=niter, warmup=burn, thin=thin, silent = 2)
# 
# # save models for future use
# saveRDS(modA, file='compiled_mods/modA')
# saveRDS(modB, file='compiled_mods/modB')
# saveRDS(modC, file='compiled_mods/modC')
# saveRDS(modD, file='compiled_mods/modD')
# saveRDS(modE, file='compiled_mods/modE')

modA = readRDS('compiled_mods/modA')
modB = readRDS('compiled_mods/modB')
modC = readRDS('compiled_mods/modC')
modD = readRDS('compiled_mods/modD')
modE = readRDS('compiled_mods/modE')

# Create output matrix
columns = c('n','j','t0','s', 
            'dicA','waicA','looA',
            'dicB','waicB','looB',
            'dicC','waicC','looC',
            'dicD','waicD','looD',
            'dicE','waicE','looE'
)

out = matrix(NA, nrow(conds)*S, length(columns), dimnames = list(NULL, columns))
idx = 1        # row counter
save_every = 1 # save output every N rows

# Go! ---------------------------------------------------------------------

for(i in seq(nrow(conds))) {
  
  cat('\n\n**************** Simulating condition', i, '/',nrow(conds),'****************\n')
  
  # Step 1. simulate S samples of data
  n  = conds$N[i]
  j  = conds$J[i]
  t0 = conds$TAU0[i]
  
  DM    = expand.grid(1:j,1:n) # design matrix
  DM$X1 = rnorm(nrow(DM))      # standard normal predictor X1
  DM$X2 = rnorm(nrow(DM),2,1)  # predictor X2 mu = 2, s =1
  DM$X3 = rnorm(nrow(DM),0,2)  # non-predictor mu = 0, s =2
  
  colnames(DM) = c('obs','id','X1','X2','X3')
  
  # specify parameters
  p = list(GAMMA=GAMMA, 
           TAU=matrix(c(t0,TAU[[2]],TAU[[2]],TAU[[1]]),2), 
           SIGMA=SIGMA)

  formula = y ~ X1 + X2 + (X1|id)
  y = simulate_lmer(formula, DM, p, nsim=S)
  
  # at high N,J, check to see pars are recovered
  # DM$y = y[,1]
  # summary(lmer(y ~ X1 + X2 + (X1|id), data=DM)) 
  
  pb  = txtProgressBar(1, S) # keep track of progress
  for(s in 1:S) {
    setTxtProgressBar(pb,s) # update progress
    
    DM$y = y[,s]
    
    # Step 2. Fit (update) models
    suppressMessages(
      {
        modA = update(modA, newdata=DM, recompile = F, refresh=0)
        modB = update(modB, newdata=DM, recompile = F, refresh=0)
        modC = update(modC, newdata=DM, recompile = F, refresh=0)
        modD = update(modD, newdata=DM, recompile = F, refresh=0)
        modE = update(modE, newdata=DM, recompile = F, refresh=0)
      }
    )
    
    # Step 3. Compute fit indices for each model
    suppressWarnings(
      {
        dicA  = DIC_brms(modA)['DIC']
        waicA = WAIC(modA)$estimates['waic','Estimate']
        looA  = LOO(modA)$estimates['looic','Estimate']
        
        dicB  = DIC_brms(modB)['DIC']
        waicB = WAIC(modB)$estimates['waic','Estimate']
        looB  = LOO(modB)$estimates['looic','Estimate']
        
        dicC  = DIC_brms(modC)['DIC']
        waicC = WAIC(modC)$estimates['waic','Estimate']
        looC  = LOO(modC)$estimates['looic','Estimate']
        
        dicD  = DIC_brms(modD)['DIC']
        waicD = WAIC(modD)$estimates['waic','Estimate']
        looD  = LOO(modD)$estimates['looic','Estimate']
        
        dicE  = DIC_brms(modE)['DIC']
        waicE = WAIC(modE)$estimates['waic','Estimate']
        looE  = LOO(modE)$estimates['looic','Estimate']
      }
    )
    
    
    # Step 4. Save output for this iteration
    out[idx, ] = c(n, j, t0, s,
                   dicA, waicA, looA,
                   dicB, waicB, looB,
                   dicC, waicC, looC,
                   dicD, waicD, looD,
                   dicE, waicE, looE)
    
    idx = idx + 1 # update counter
  }
  # save temporarily, just in case something crashes
  if(i%%save_every==0) write.csv(data.frame(out), paste0('out/.sim_out_tmp_i=',i,'.csv'))
}

out = as.data.frame(out) # for convenience
head(out)
write.csv(out,file = 'out/sim_out.csv')



