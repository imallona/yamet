#!/usr/bin/env R
##
## not snakemaked
##
## requires intersect_to_correlate.sh to be run first

library(ggplot2)

fns <- list.files('.', pattern = 'intersect_*')

## fns <- grep('*VFH*VFH', '*VFH*RVP*', '*VFH*ELB*, '' 

d <- list()

for (fn in fns) {
    d[[fn]] <- read.table(fn, header = TRUE)
}

## par(mfrow = c(3,3))
## for (item in d) {
##     plot(item[,1], item[,3])
## }

## for (hela in c('RVP', 'VFH')) {
##     for (gm in c('ELB', 'RGH')) {
##         for (imr in 'TXF') {

##         }
##     }
## }

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

hela <-  c('RVP', 'VFH')
gm <- c('ELB', 'RGH')
imr <- 'TXF'


## entropies

par(mfrow = c(5,5), mai = c(0.01, 0.01, 0.01, 0.01), pty = 's')
for (x in  c(imr, gm, hela)) {
    for (y in c(imr, gm, hela)) {
        curr <- d[[grep(sprintf('.*%s.*%s.*', x, y), fns)]]
        plot(curr[,1], curr[,3], axes=FALSE, frame.plot=TRUE, xlab = '', ylab = '',
             pch = 21,
             col = add.alpha('black', 0.2))
        
    }
}

## methyaltions


par(mfrow = c(5,5), mai = c(0.01, 0.01, 0.01, 0.01), pty = 's')
for (x in  c(imr, gm, hela)) {
    for (y in c(imr, gm, hela)) {
        curr <- d[[grep(sprintf('.*%s.*%s.*', x, y), fns)]]
        plot(curr[,2], curr[,4], axes=FALSE, frame.plot=TRUE, xlab = '', ylab = '',
             pch = 21,
             col = add.alpha('black', 0.2))
        
    }
}


## foo <- do.call(rbind.data.frame, d)
