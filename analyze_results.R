
rm(list=ls())

# Load simulated data -----------------------------------------------------

dat = read.csv('out/sim_out.csv')

# double check all simulations are saved
table(dat$n, dat$j, dat$t0) # should all be 100


# Visualize accuracy ---------------------------------------------------------------

cols = c('darkblue','darkred','darkgreen')
save_ind = F # save individual plots (T/F)

if(!save_ind) {
  pdf('figs/acc_matplot_full.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

## DIC -----

dic = dat[,c('n','j','t0',paste0('dic',LETTERS[1:5]))]
dic$chosen = apply(dic[,4:ncol(dic)], 1, which.min)
dic$acc = as.numeric(dic$chosen==5)

pcor = tapply(dic$acc, list(dic$t0, dic$j, dic$n), mean)

for(i in 1:3) {
  if(save_ind) pdf(paste0('figs/dic_acc_matplot_n=',dimnames(pcor)[[3]][i],'.pdf'), 6, 6)
  matplot(pcor[,,i], pch=18, lty=1, type='b', ylim=c(0,1), lwd=3, 
          col=cols, cex=1.5,
          xlab = expression(tau[0]^2), xaxt='n', 
          ylab = 'P(Correct Model Selected)', 
          main = paste0('L2 Sample Size = ', dimnames(pcor)[[3]][i]),
          cex.lab=1.4, cex.axis=1.5, cex.main=1.5)
  axis(1, at = 1:3, labels=rownames(pcor[,,i]), cex.axis=1.5)
  legend('bottomright', bty='n', lty=1, pch=18, col=cols, title='L1 Sample Size', 
         legend=colnames(pcor[,,i]), lwd=2,cex=1.5 )
  # abline(h=1/5, lty=2)
  if(save_ind) dev.off()
}

dic$n_c = dic$n - median(dic$n)
dic$j_c = dic$j - median(dic$j)
dic$t0_c = dic$t0 - median(dic$t0)

mdic = glm(acc ~ n_c*j_c*t0_c, data=dic, family='binomial')
capture.output(summary(mdic), file = 'out/acc_dic_lm.txt')

## WAIC -----

waic = dat[,c('n','j','t0',paste0('waic',LETTERS[1:5]))]
waic$chosen = apply(waic[,4:ncol(waic)], 1, which.min)
waic$acc = as.numeric(waic$chosen==5)

pcor = tapply(waic$acc, list(waic$t0, waic$j, waic$n), mean)

for(i in 1:3) {
  if(save_ind) pdf(paste0('figs/waic_acc_matplot_n=',dimnames(pcor)[[3]][i],'.pdf'), 6, 6)
  matplot(pcor[,,i], pch=18, lty=1, type='b', ylim=c(0,1), lwd=3, 
          col=cols, cex=1.5,
          xlab = expression(tau[0]^2), xaxt='n', 
          ylab = 'P(Correct Model Selected)', 
          main = paste0('L2 Sample Size = ', dimnames(pcor)[[3]][i]),
          cex.lab=1.4, cex.axis=1.5, cex.main=1.5)
  axis(1, at = 1:3, labels=rownames(pcor[,,i]), cex.axis=1.5)
  legend('bottomright', bty='n', lty=1, pch=18, col=cols, title='L1 Sample Size', 
         legend=colnames(pcor[,,i]), lwd=2,cex=1.5 )
  # abline(h=1/5, lty=2)
  if(save_ind) dev.off()
}

waic$n_c = waic$n - median(waic$n)
waic$j_c = waic$j - median(waic$j)
waic$t0_c = waic$t0 - median(waic$t0)

mwaic = glm(acc ~ n_c*j_c*t0_c, data=waic, family='binomial')
capture.output(summary(mwaic), file = 'out/acc_waic_lm.txt')

## LOO-IC -----

loo = dat[,c('n','j','t0',paste0('loo',LETTERS[1:5]))]
loo$chosen = apply(loo[,4:ncol(loo)], 1, which.min)
loo$acc = as.numeric(loo$chosen==5)

pcor = tapply(loo$acc, list(loo$t0, loo$j, loo$n), mean)

for(i in 1:3) {
  if(save_ind) pdf(paste0('figs/loo_acc_matplot_n=',dimnames(pcor)[[3]][i],'.pdf'), 6, 6)
  matplot(pcor[,,i], pch=18, lty=1, type='b', ylim=c(0,1), lwd=3, 
          col=cols, cex=1.5,
          xlab = expression(tau[0]^2), xaxt='n', 
          ylab = 'P(Correct Model Selected)', 
          main = paste0('L2 Sample Size = ', dimnames(pcor)[[3]][i]),
          cex.lab=1.4, cex.axis=1.5, cex.main=1.5)
  axis(1, at = 1:3, labels=rownames(pcor[,,i]), cex.axis=1.5)
  legend('bottomright', bty='n', lty=1, pch=18, col=cols, title='L1 Sample Size', 
         legend=colnames(pcor[,,i]), lwd=2,cex=1.5 )
  # abline(h=1/5, lty=2)
  if(save_ind) dev.off()
}

if(!save_ind) dev.off()

loo$n_c = loo$n - median(loo$n)
loo$j_c = loo$j - median(loo$j)
loo$t0_c = loo$t0 - median(loo$t0)

mloo = glm(acc ~ n_c*j_c*t0_c, data=loo, family='binomial')
summary(mloo)
capture.output(summary(mloo), file = 'out/acc_loo_lm.txt')
  
# Visualize Chosen Models -------------------------------------------------

## DIC ----

if(!save_ind) {
  pdf('figs/choose_bar_dic.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(dic[,paste0('dic',LETTERS[i])], list(waic$t0, waic$j, waic$n), mean)
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/dic_choose_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(0,1), xaxt='n',
                xlab = 'Model', ylab='P(Model Selected)',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()
  
## WAIC ----

if(!save_ind) {
  pdf('figs/choose_bar_waic.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(waic$chosen, list(waic$t0, waic$j, waic$n), function(j) mean(j==i))
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/waic_choose_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(0,1), xaxt='n',
                xlab = 'Model', ylab='P(Model Selected)',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()

## LOO-IC ----

if(!save_ind) {
  pdf('figs/choose_bar_loo.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(loo$chosen, list(loo$t0, loo$j, loo$n), function(j) mean(j==i))
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/loo_choose_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(0,1), xaxt='n',
                xlab = 'Model', ylab='P(Model Selected)',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()


# Visualize Raw Scores ----------------------------------------------------


## DIC ----

if(!save_ind) {
  pdf('figs/val_bar_dic.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(dic[,paste0('dic',LETTERS[i])], list(waic$t0, waic$j, waic$n), sum)
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/dic_val_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(min(pchoose[,i,j,]-1000), max(pchoose[,i,j,]+1000)),xpd=F, xaxt='n',
                xlab = 'Model', ylab='Cum. DIC',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()

## WAIC ----
if(!save_ind) {
  pdf('figs/val_bar_waic.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(waic[,paste0('waic',LETTERS[i])], list(waic$t0, waic$j, waic$n), sum)
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/waic_val_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(min(pchoose[,i,j,]-1000), max(pchoose[,i,j,]+1000)),xpd=F, xaxt='n',
                xlab = 'Model', ylab='Cum. WAIC',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()

## LOO-IC ----

if(!save_ind) {
  pdf('figs/val_bar_loo.pdf', 12, 12)
  layout(matrix(1:9, 3, 3, byrow = T))
}

pchoose = array(NA, c(3,3,3,5), dimnames = list(c(.05,1,5), c(10,50,100), c(10,50,100), LETTERS[1:5]))
for(i in 1:5) {
  pchoose[,,,i] = tapply(loo[,paste0('loo',LETTERS[i])], list(waic$t0, waic$j, waic$n), sum)
}

for(i in 1:3) {
  for(j in 1:3) {
    if(save_ind) pdf(paste0('figs/loo_val_bar_n=',dimnames(pcor)[[2]][i],'_j=',dimnames(pcor)[[3]][j],'.pdf'), 6, 6)
    
    b = barplot(pchoose[,i,j,], beside=T, ylim=c(min(pchoose[,i,j,]-1000), max(pchoose[,i,j,]+1000)),xpd=F, xaxt='n',
                xlab = 'Model', ylab='Cum. LOO-IC',
                main = paste0('L1 Sample Size = ', dimnames(pchoose)[[2]][i], '\nL2 Sample Size = ',dimnames(pchoose)[[3]][j]), 
                cex.lab=1.4, cex.axis=1.5, cex.main=1.5, 
                legend.text = T, 
                args.legend = list(x='topleft', bty='n', title=expression(tau[0]^2),cex=1.5)
    )
    
    axis(1, at=colMeans(b), labels=LETTERS[1:5], cex.axis=1.5)
    
    if(save_ind) dev.off()
  }
}

if(!save_ind) dev.off()
