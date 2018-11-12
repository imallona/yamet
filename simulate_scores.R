#!/usr/bin/env R
##
## Cometh agreement scoring simulations
##
## Izaskun Mallona
## 12th nov 2018

set.seed(1)
n = 500
x = round(rbeta(n, 0.5, 0.2), 2)

hist(x, prob = TRUE, main = 'beta')
rug(x)
curve(dbeta(x, .5, .4), lwd=2, col="blue", add=TRUE)

## cov <- floor(abs(rnorm(n,20,10)))
    
## fd <- data.frame(beta = x,
##                  cov = cov)
## fd$cov_meth <- floor(fd$beta * fd$cov)

## ## a third of the cases are MM tuples
## fd$mm <- fd$um <- fd$mu <- floor(fd$cov_meth/3)
## fd$uu <- fd$cov - (fd$mm + fd$um + fd$mu)

## summary(fd)

## hist(fd$cov_meth)
## plot(fd)



    
## fd <- data.frame(beta = x,
##                  mm = floor(abs(rnorm(n,5,2))),
##                  um = floor(abs(rnorm(n,5,2))),
##                  mu =  floor(abs(rnorm(n,5,2))))

## fd$meth <- (2*fd$mm) + fd$um + fd$mu
## fd$cov <- fd$meth / fd$beta
## fd$uu <- fd$cov - (fd$mm + fd$mu + fd$mm)

## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

## plot(fd)

## ## the beta vs coverage is shitty


    
fd <- data.frame(beta = x,
                 cov = floor(abs(rnorm(n,25,3))),
                 mm = floor(abs(rnorm(n,5,2))),
                 uu =  floor(abs(rnorm(n,5,2))))

                 

fd$missing_meth <- fd$cov* fd$beta
fd$mu <- round(fd$missing_meth/2)
fd$um <- round(fd$missing_meth/2)
    
## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

plot(fd)

## the beta vs coverage is shitty

## what if sampling them?
