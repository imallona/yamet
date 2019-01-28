#!/bin/bash

library(entropy)


load('/home/imallona/Desktop/test.RData')
d <- foo


## new approach


d$pUU <- (d$cov - (d$cov * d$beta))/d$cov
d$pMM <- 1 - d$pUU

head(d)

## d$max_entropy <- apply(data.frame(pUU = 2 * d$pUU * d$coverage,
##                                   pUM = (d$pUU * d$coverage) + (d$pMM * d$coverage),
##                                   pMU = (d$pUU * d$coverage) + (d$pMM * d$coverage),
##                                   pMM = 2 * d$pMM * d$coverage),
##                        1, entropy)

d$max_entropy <- apply(data.frame(pUU = d$pUU^2,
                                  pUM = d$pUU * d$pMM,
                                  pMU = d$pUU * d$pMM,
                                  pMM = d$pMM^2),
                       1, entropy)


head(d)


## plot(d$beta, d$entropy)

## plot(y = d$entropy, x= d$max_entropy)


## d$adj_entropy <- (d$entropy - d$max_entropy)
targets <- c('coverage',
             'distance',
             'beta',
             'entropy',
             'max_entropy')

plot(d[,targets])


## now min entropy (to standardize with min/max values)


d$min_entropy_halfs <- NULL

d$min_entropy_down <- apply(data.frame(pUU = d$pUU,
                                  pUM = 0,
                                  pMU = 0,
                                  pMM = d$pMM),
                       1, entropy)

head(d[d$beta == 0.6,])
summary(d$pUU)
summary(d$pMM)

d$min_entropy_halfs <- apply(data.frame(pUU = d$pUU^2,                                        
                                        pUM = d$pMM^2,
                                        pMU = 0,
                                        pMM = 0),
                             1, entropy)
d$min_entropy_halfs[is.na(d$min_entropy_halfs)] <- 0

summary(d$min_entropy_down)
summary(d$min_entropy_halfs)

targets <- c('beta',
             'entropy',
             'max_entropy',
             'min_entropy_down',
             'min_entropy_halfs')


plot(d[,targets], xlim = c(0,1.5), ylim = c(0, 1.5))
