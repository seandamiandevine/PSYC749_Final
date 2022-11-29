
rm(list=ls())
library(lme4)
library(mnormt)
source('fx.R')

# Set constants -----------------------------------------------------------

set.seed(2022) 

# Simulation structure
# N,J should be large to properly recovery parameters
N = 100  # number of clusters
J = 250 # obs/cluster
S = 100 # number of samples


# Design matrix
DM    = expand.grid(obs=1:J,id=1:N) # design matrix
DM$X1 = rnorm(nrow(DM))             # standard normal predictor X1
DM$X2 = rnorm(nrow(DM),2,1)         # predictor X2 mu = 2, s =1
DM$X3 = rnorm(nrow(DM),0,2)         # non-predictor mu = 0, s =2

formula = y ~ X1 + X2 + (X1|id)
recov   = matrix(NA, S, 17, 
                 dimnames = list(NULL, c('iter',
                                         'G00','G10','G20','TAU0','TAU1','SIGMA',
                                         'g00','g10','g20','tau0','tau1','sigma',
                                         'dgam','dtau','dsig','dtot')))

pb = txtProgressBar(1, S)
for(i in 1:S) {
  
  setTxtProgressBar(pb, i)
  
  # 1. Simulate from population parameters
  GAMMA = rnorm(3, 0, 5)           # fixed effects (g00, g10, g20)
  TAU   = matrix(c(rgamma(1,2,2),
                   0,0,
                   rgamma(1,2,2)),2)  # VarCor Matrix
  SIGMA = rgamma(1,3,3)               # residual variance
  p = list(GAMMA=GAMMA, TAU=TAU, SIGMA=SIGMA)
  
  # 1. Fit true model with lme4
  DM$y = sim1(formula, DM, p, seed = round(runif(1,0,1e5)))
  mod  = lmer(y ~ X1 + X2 + (X1|id), data=DM)
  
  # 2. Extract values of interest
  gam = unname(fixef(mod))
  tau = unname(attr(VarCorr(mod)$id,'stddev')^2)
  sig = unname(sigma(mod)^2)
  
  # 3. Compute difference from true value
  dgam = sum(abs(gam - GAMMA))
  dtau = sum(abs(tau - diag(TAU)))
  dsig = abs(sig - SIGMA)
  dtot = dgam + dtau + dsig
  
  recov[i,] = c(i, GAMMA, diag(TAU), SIGMA, gam, tau, sig, dgam, dtau, dsig, dtot)
  
}

write.csv(as.data.frame(recov), 'supp/par_recov.csv')

# Visualize ---------------------------------------------------------------

# recov = read.csv('supp/par_recov.csv')

pars = c('G00','G10','G20','TAU0','TAU1','SIGMA')

pdf('supp/recov.pdf', 10,10)
layout(matrix(1:6,3, 2,byrow=T))

for(par in pars) {
  
  plot(recov[,par], recov[,tolower(par)], ylab='Estimated', xlab='True', main=par,
       cex=1.5, cex.main=1.5, cex.lab=1.5, cex.axis=1.5, pch=3, col=scales::alpha('red',.5))
  
  r = cor.test(recov[,par], recov[,tolower(par)])
  legend('topleft',bty='n',legend=paste0('r = ',round(r$estimate,4), '\np = ',round(r$p.value,4)), 
         cex = 1.5)
  
}

dev.off()

