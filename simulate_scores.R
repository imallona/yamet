#!/usr/bin/env R
##
## Cometh agreement scoring simulations
##
## Izaskun Mallona
## 12th nov 2018

library('latex2exp')

WD <- file.path('/home/imallona', 'cg_shadows', 'data')
fn <- file.path(WD, 'sorjuela.CG.2.tsv')

d <- read.table(fn, header = TRUE, nrows = 500000)
summary(d)
d$coverage <- d$MM + d$UM + d$MU + d$UU
d$beta <- (d$MM + 0.5*d$UM + 0.5*d$MU) / d$coverage

d$score1 <- (d$MM + d$UU) / d$coverage
d$score2 <- (d$MM) / d$coverage
d$score3 <- (2*d$MM + d$UU) / d$coverage

targets <- c('coverage', 'beta', 'score1', 'score2', 'score3')


png(file.path(WD, 'exploratory_1.png'), width = 700, height = 700)

plot(d[,targets])
dev.off()

pdf(file.path(WD, 'exploratory_2.pdf'), width = 14)


## scatterplotMatrix(d[,targets], main="", diagonal = 'density', plot.points = FALSE)

NPLOT = 150


## score3


## plot(y = jitter(d$score3[1:NPLOT]),
##      x = jitter(d$beta[1:NPLOT]),
##      main = TeX('$score3= \\frac{2*MM + UU}{cov} vs \\beta$'),
##      xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
##      ylab = 'score3 with jitter') 


plot(y = d$score3[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score3= \\frac{2*MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score3[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score3 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main =  TeX('$score3= \\frac{2*MM + UU}{cov}$'))




## on score1
plot(y = jitter(d$score1[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     main = TeX('$score1= \\frac{MM + UU}{cov} vs \\beta$'),
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     ylab = 'score1 with jitter') 


plot(y = d$score1[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score1',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score1= \\frac{MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score1[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score1 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main = TeX('$score1= \\frac{MM + UU}{cov}$'))





## on score2
plot(y = jitter(d$score2[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     main = TeX('$score2= \\frac{MM}{cov} vs \\beta$'),
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     ylab = 'score2 with jitter') 


plot(y = d$score2[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score2',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score2= \\frac{MM}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score2[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score2 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main = TeX('$score2= \\frac{MM}{cov}$'))


dev.off()






## lineplots depicting score 1, 3 and beta



plot(y = d$score1[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score1',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT),
     ylim = c(0,1))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = d$score2[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 1, col = 'darkgreen')

legend("topright",
       legend = c(TeX('$score1= \\frac{MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  '$score2 = \\frac{MM}{cov}$'),
       col = c('black', 'blue', 'darkgreen'),
       lwd = c(2,2,2),
       lty = c(1,1,1),
       ## xjust = 1, yjust = 1,
       title = "Scores")

     
## plot(y = jitter(d$score1[1:NPLOT] - d$beta[1:NPLOT]), x = jitter(d$beta[1:NPLOT]))


## plot(y = d$score1[1:NPLOT], x = 1:NPLOT, type = 'l')
## lines(y = d$beta[1:NPLOT], x = 1:NPLOT, type = 'l',  col = 'blue')
    



## set.seed(1)
## n = 500
## x = round(rbeta(n, 0.5, 0.2), 2)

## hist(x, prob = TRUE, main = 'beta')
## rug(x)
## curve(dbeta(x, .5, .4), lwd=2, col="blue", add=TRUE)

## ## cov <- floor(abs(rnorm(n,20,10)))
    
## ## fd <- data.frame(beta = x,
## ##                  cov = cov)
## ## fd$cov_meth <- floor(fd$beta * fd$cov)

## ## ## a third of the cases are MM tuples
## ## fd$mm <- fd$um <- fd$mu <- floor(fd$cov_meth/3)
## ## fd$uu <- fd$cov - (fd$mm + fd$um + fd$mu)

## ## summary(fd)

## ## hist(fd$cov_meth)
## ## plot(fd)



    
## ## fd <- data.frame(beta = x,
## ##                  mm = floor(abs(rnorm(n,5,2))),
## ##                  um = floor(abs(rnorm(n,5,2))),
## ##                  mu =  floor(abs(rnorm(n,5,2))))

## ## fd$meth <- (2*fd$mm) + fd$um + fd$mu
## ## fd$cov <- fd$meth / fd$beta
## ## fd$uu <- fd$cov - (fd$mm + fd$mu + fd$mm)

## ## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

## ## plot(fd)

## ## ## the beta vs coverage is shitty


    
## fd <- data.frame(beta = x,
##                  cov = floor(abs(rnorm(n,25,3))),
##                  mm = floor(abs(rnorm(n,5,2))),
##                  uu =  floor(abs(rnorm(n,5,2))))

                 

## fd$missing_meth <- fd$cov* fd$beta
## fd$mu <- round(fd$missing_meth/2)
## fd$um <- round(fd$missing_meth/2)
    
## ## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

## plot(fd)

## ## the beta vs coverage is shitty

## ## what if sampling them?
