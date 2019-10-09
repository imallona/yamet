#!/usr/bin/env R
##
## not snakemaked
##
## requires intersect_to_correlate.sh to be run first

## library(ggplot2)
library(Cairo)


add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

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
hela <-  c('RVP', 'VFH')
gm <- c('ELB', 'RGH')
imr <- 'TXF'


samples <- c(imr, hela, gm)
annotation <- c( 'IMR90', 'HeLa-S3', 'HeLa-S3', 'GM23248', 'GM23248')

 ## bottom, left, top and right

CairoPNG('correlation_encode_entropies.png', width= 600, height = 600)

par(mfrow = c(5,5),
    mai = c(0.01, 0.01, 0.01, 0.01),
    pty = 's',
    mar=c(1, 1, 1, 1))

i <- 1

for (x in 1:length(samples)) {
    for (y in 1:length(samples)) {
        idx = i %% length(samples)

        curr <- d[[grep(sprintf('.*%s.*%s.*', samples[x], samples[y]), fns)]]
       
        print(sprintf('i = %s, x= %s, y = %s, idx = %s', i, x, y, idx))
        if (x == y) {
            ## plot(1)
            ## plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            ## text(x = 1, y = 1, sprintf('%s\n%s', annotation[x],samples[x]), cex = 1.4)

            plot(curr[,1], curr[,3], type='p', xlab = '', ylab = '', col ='white')
            
            text(x = 0.7, y = 0.7, sprintf('%s\n%s', annotation[x],samples[x]), cex = 1.4)
            
        } else if (y > x) {
           
            
            plot(curr[,1], curr[,3], axes=FALSE, frame.plot=TRUE, xlab = '', ylab = '',
                 pch = 19,
                 col = add.alpha('black', 0.05),
                 cex = 0.5)            
        } else {
            plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            text(x = 1, y = 1,
                 sprintf('rho = %s\nn =  %s',
                         round(cor(curr[,1], curr[,3], method = 'spearman'), 3),
                         nrow(curr)),
                 cex = 1.4)
        }
        i = i +1
    }
}

## mtext("standardized entropy (stdH)", side = 1, line = -2, adj = 3.5, cex = 1.2)
mtext("Shannon's entropy (H)", side = 1, line = -2, outer = TRUE,  cex = 1.2)
## mtext("standardized entropy (stdH)", side = 2, line = 40, adj =40)
mtext("Shannon's entropy (H)", side = 2, line = -2, outer = TRUE, cex = 1.2)

dev.off()



CairoPNG('correlation_encode_betas.png', width= 600, height = 600)

par(mfrow = c(5,5),
    mai = c(0.01, 0.01, 0.01, 0.01),
    pty = 's',
    mar=c(1, 1, 1, 1))

i <- 1

for (x in 1:length(samples)) {
    for (y in 1:length(samples)) {
        idx = i %% length(samples)

        curr <- d[[grep(sprintf('.*%s.*%s.*', samples[x], samples[y]), fns)]]
       
        print(sprintf('i = %s, x= %s, y = %s, idx = %s', i, x, y, idx))
        if (x == y) {
            ## plot(1)
            ## plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            plot(curr[,2], curr[,4], type='p', xlab = '', ylab = '', col ='white')
            
            text(x = 0.5, y = 0.5, sprintf('%s\n%s', annotation[x],samples[x]), cex = 1.4)
            
        } else if (y > x) {
           
            
            plot(curr[,2], curr[,4], axes=FALSE, frame.plot=TRUE, xlab = '', ylab = '',
                 pch = 19,
                 col = add.alpha('black', 0.05),
                 cex = 0.5)            
        } else {
            plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            text(x = 1, y = 1,
                 sprintf('rho = %s\nn =  %s',
                         round(cor(curr[,2], curr[,4], method = 'spearman'), 3),
                         nrow(curr)),
                 cex = 1.4)
        }
        i = i +1
    }
}

## mtext("standardized entropy (stdH)", side = 1, line = -2, adj = 3.5, cex = 1.2)
mtext("DNA methylation (beta)", side = 1, line = -2, outer = TRUE,  cex = 1.2)
## mtext("standardized entropy (stdH)", side = 2, line = 40, adj =40)
mtext("DNA methylation (beta)", side = 2, line = -2, outer = TRUE, cex = 1.2)

dev.off()



## cross beta vs entropy across replicates


CairoPNG('correlation_encode_betas_vs_entropies.png', width= 600, height = 600)

par(mfrow = c(5,5),
    mai = c(0.01, 0.01, 0.01, 0.01),
    pty = 's',
    mar=c(1, 1, 1, 1))

i <- 1

for (x in 1:length(samples)) {
    for (y in 1:length(samples)) {
        idx = i %% length(samples)

        curr <- d[[grep(sprintf('.*%s.*%s.*', samples[x], samples[y]), fns)]]
       
        print(sprintf('i = %s, x= %s, y = %s, idx = %s', i, x, y, idx))
        if (x == y) {
            ## plot(1)
            ## plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            plot(x = curr[,1], y = curr[,2], type='p', xlab = '', ylab = '', col ='white')
            
            text(x = 0.7, y = 0.5, sprintf('%s\n%s', annotation[x],samples[x]), cex = 1.4)
            
        } else if (y > x) {
           
            
            plot(curr[,1], curr[,2], axes=FALSE, frame.plot=TRUE, xlab = '', ylab = '',
                 pch = 19,
                 col = add.alpha('black', 0.05),
                 cex = 0.5)            
        } else {
            plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0,2), ylim = c(0,2))
            text(x = 1, y = 1,
                 sprintf('rho = %s\nn =  %s',
                         round(cor(curr[,1], curr[,2], method = 'spearman'), 3),
                         nrow(curr)),
                 cex = 1.4)
        }
        i = i +1
    }
}

## mtext("standardized entropy (stdH)", side = 1, line = -2, adj = 3.5, cex = 1.2)
mtext("Shannon's Entropy (H)", side = 1, line = -2, outer = TRUE,  cex = 1.2)
## mtext("standardized entropy (stdH)", side = 2, line = 40, adj =40)
mtext("DNA methylation (beta)", side = 2, line = -2, outer = TRUE, cex = 1.2)

dev.off()



