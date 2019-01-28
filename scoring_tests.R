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

## d$min_entropy_halfs <- apply(data.frame(pUU = 0,                                        
##                                         pUM = abs(d$pUU - 0.5),
##                                         pMU = abs(d$pMM - 0.5),
##                                         pMM = 0),
##                              1, entropy)
## d$min_entropy_halfs[is.na(d$min_entropy_halfs)] <- 0


## d$min_entropy_halfs <- apply(data.frame(pUU = 0,                                        
##                                         pUM = abs(d$pUU - 0.5)*d$pUU * abs(d$pMM - 0.5)*d$pMM ,
##                                         pMU = abs(d$pUU - 0.5)*d$pUU * abs(d$pMM - 0.5)*d$pMM ,
##                                         pMM = 0),
##                              1, entropy)
## d$min_entropy_halfs[is.na(d$min_entropy_halfs)] <- 0


## d$min_entropy_halfs <- apply(data.frame(pUU = 0,                                        
##                                         pUM = abs(d$pUU - 0.5)*d$pUU * abs(d$pMM - 0.5)*d$pMM ,
##                                         pMU = abs(d$pUU - 0.5)*d$pUU * abs(d$pMM - 0.5)*d$pMM ,
##                                         pMM = 0),
##                              1, entropy)
## d$min_entropy_halfs[is.na(d$min_entropy_halfs)] <- 0


## summary(d$beta)
## summary(d$min_entropy_down)
## summary(d$min_entropy_halfs)

## targets <- c('beta',
##              'entropy',
##              'max_entropy',
##              'min_entropy_down',
##              'min_entropy_halfs')

max(d$min_entropy_down)
head(d[d$entropy < d$min_entropy_down,])


plot(d[,targets], xlim = c(0,1.5), ylim = c(0, 1.5))



d$test<- apply(data.frame(pUU = 0,
                          pUM = 0,
                          pMU = 2* (d$pUU * d$pMM),
                          pMM = d$pMM^2),
               1, entropy)

d$test[is.na(d$test)] <- 0

## d$test2 <- apply(data.frame(pUU = d$pUU,
##                           pUM = 0,
##                           pMU = d$pUU * d$pMM,
##                           pMM = 0),
##                1, entropy)

d$test2 <- apply(data.frame(pUU = d$pUU^2,
                          pUM = 0,
                          pMU = 2 * (d$pUU * d$pMM),
                          pMM = 0),
               1, entropy)


d$test2[is.na(d$test2)] <- 0

d$test3 <- apply(data.frame(pUU = d$pUU,
                          pUM = abs(d$pUU-d$pMM),
                          pMU = 0,
                          pMM = 0),
               1, entropy)


d$test3[is.na(d$test3)] <- 0
summary(d$test3)

d$test4 <- apply(data.frame(pUU = max(d$pUU, d$pMM),
                          pUM = 0,
                          pMU = abs(d$pMM - d$pUU)*4*abs(d$pMM - d$pUU),
                          pMM = 0),
               1, entropy)


## d$test4 <- apply(data.frame(pUU = d$pUU,
##                           pUM = 0,
##                           pMU = (abs(d$pMM - d$pUU)*2)^2,
##                           pMM = 0),
##                  1, entropy)


d$test4[is.na(d$test4)] <- 0
summary(d$test4)



d$test5 <- apply(data.frame(pUU = d$pUU,
                          pUM = abs(d$pMM - d$pUU)*2,
                          pMU = abs(d$pMM - d$pUU)*2,
                          pMM = d$pMM),
                 1, entropy)


d$test5[is.na(d$test5)] <- 0
summary(d$test5)






targets <- c('beta',
             'entropy',
             'max_entropy',
             'min_entropy_down',
             'test',
             'test2',
             'test3')

## summary(d$entropy)

## plot(d[,targets], xlim = c(0,1.5), ylim = c(0, 1.5))

## plot(d$beta, d$entropy)
## lines(d$beta, d$max_entropy, col = 'brown', type = 'p', pch = 19, cex = 0.5)
## lines(d$beta, d$min_entropy_down, col = 'blue', type = 'p', pch = 19, cex = 0.5)
## lines(d$beta, d$test, col = 'darkred', type = 'p', pch = 19, cex = 0.5)
## lines(d$beta, d$test2, col = 'darkgreen', type = 'p', pch = 19, cex = 0.5)
## lines(d$beta, d$test3, col = 'gold', type = 'p', pch = 19, cex = 0.5)
## lines(d$beta, d$test4, col = 'orange', type = 'p', pch = 19, cex = 0.5)
## plot(d[,targets], xlim = c(0,1.5), ylim = c(0, 1.5))


plot(d$beta, d$entropy)
lines(d$beta, d$max_entropy, col = 'brown', type = 'p', pch = 19, cex = 0.5)
lines(d$beta, d$min_entropy_down, col = 'blue', type = 'p', pch = 19, cex = 0.5)
lines(d$beta, d$test4, col = 'orange', type = 'p', pch = 19, cex = 0.5)


## so let's standardize by the minimum

d$lowerbound <- apply(d[,c('min_entropy_down', 'test4')], 1, min)
summary(d$lowerbound)


plot(d$beta, d$entropy)
lines(d$beta, d$max_entropy, col = 'brown', type = 'p', pch = 19, cex = 0.5)
lines(d$beta, d$lowerbound, col = 'blue', type = 'p', pch = 19, cex = 0.5)



d$standardized_entropy <- (d$entropy - d$lowerbound) / (d$max_entropy - d$lowerbound)

d$standardized_entropy[is.na(d$standardized_entropy)] <- 0
summary(d$standardized_entropy)

targets <- c('beta',
             'entropy',
             'max_entropy',
             'lowerbound',
             'standardized_entropy')
summary(d[,targets])

plot(d[,targets])

cor(d$beta, d$standardized_entropy)
cor(d$beta, d$entropy)
