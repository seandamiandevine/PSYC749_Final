

sim1 = function(formula, DM, pars, seed=2022) {
  #' Simulate one lmer response from a formula and design matrix
  #' @param formula: R formula object, in lme4 style
  #' @param DM: Design matrix, with all relevant variables present
  #' @param seed: Random seed.
  #'  
  #' @returns vector with simulated outcome
  
  
  set.seed(seed)  
  
  # 1. Extract pars
  gam = pars$GAMMA
  v   = pars$TAU
  sig = pars$SIGMA
  
  # 2. Extract vars from formula
  preds = attributes(terms(formula))$term.labels
  f     = preds[-length(preds)]
  r     = gsub(' ', '', strsplit(preds[length(preds)],'\\|')[[1]])
  g     = r[length(r)]
  r     = r[-length(r)]
  
  dmat = cbind(DM[,g], 1, DM[,f])
  n    = length(unique(dmat[,1]))
  
  # 3. Sample random effects (level-2)
  Uj = rmnorm(n, mean = rep(0, length(r)+1), varcov = v) # intercept, slopes
  
  # 4. Sample residuals (level-1)
  Rij = rnorm(nrow(dmat), 0, sqrt(sig))
  
  # 5. Compute betas as GAMMA + Uj
  # First, add zero column for fixed effects only
  gam_m = matrix(gam, n,length(gam), byrow = TRUE)
  if(ncol(Uj)!=ncol(gam_m)) {
    diff = ncol(gam_m) - ncol(Uj)
    if(diff < 0) stop('VarCov mispecified')
    for(i in 1:diff) Uj = cbind(Uj, 0)
  }
  b = gam_m + Uj
  
  # 6. Compute response
  y = rowSums(dmat[,-1]*b[dmat[,1],]) + Rij
  
  y
  
}

simulate_lmer = function(formula, DM, pars, nsim, seed=2022) {
  #' Wrapper for sim1(), which reproduces the process `nsim` times
  
  set.seed(seed)
  sapply(1:nsim, function(i) sim1(formula, DM, pars, runif(1,0,1e6)))
  
}

DIC_brms = function(mod, version=1) {
  
  #' Convenience function to compute the DIC for a brms object
  #' Computed at the mean posterior parameter estimates, as implemented in Gelman et al., 2014, Stats. & Comp.
  #' @param mod: a fitted brms object
  #' @param version: how to compute pDIC (see Gelman et al., 2014, Stats. & Comp.)
  #' @return [vector] expected log predictive density, effective number of parameters, dic value
  
  # p(y | MAP)
  y    = mod$data[,1]
  yhat = predict(mod)[,'Estimate']
  sig  = summary(mod)$spec_pars[,'Estimate']
  pMAP = sum(dnorm(y, yhat, sig, log=T))
  
  # pDIC
  ll = log_lik(mod) # point-wise likelihood
  
  if(version==1) {
    mll   = mean(rowSums(ll))    # mean deviance
    pDIC  = 2*(pMAP-mll)
  } else {
    pDIC = 2*var(rowSums(ll))
  }
  
  # elpd
  elpd = pMAP - pDIC
  
  # compute dic
  return(c(elpd = elpd, pDIC=pDIC, DIC=-2*pMAP + 2*pDIC))
  
}
