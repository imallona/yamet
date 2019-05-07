#!/bin/env R
##
## sample-based meth report
## Izaskun Mallona
##
## Thu May  2 2019

suppressPackageStartupMessages({
    library(data.table)
    library('latex2exp')
    library(entropy)
    ## library(cramer)
    ## library(Cairo)
    ## library(knitr)
    ## library(ggplot2)
    ## library(RColorBrewer)
    ## library(gridExtra)
    library(dplyr)
    ## library(pheatmap)
    ## library(ggExtra)
    ## library(tidyverse)
    ## library(viridis)
})



args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(args)
print(output)


ac  <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    m
}

########


## fn <-  file.path(snakemake@wildcards$sample,
##                              sprintf('%s_cov_%s_hmm_%s.test',
##                                      snakemake@wildcards$sample,
##                                      snakemake@wildcards$cov,
##                                      snakemake@wildcards$hmm))

fn <- colored_entropy
getwd()
print(fn)

system(sprintf('file %s', fn))

d <- data.table::fread(file.path(getwd(), fn), sep = '\t')


head(d)

hmm_colors <- c('13_ReprPC', '1_TssA', '4_Tx', '14_ReprPCWk',
                '2_TssAFlnk',
                '15_Quies',
                '5_TxWk',
                '7_Enh',
                '8_ZNF/Rpts',
                '3_TxFlnk',
                '6_EnhG',
                '9_Het',
                '11_BivFlnk',
                '10_TssBiv',
                '12_EnhBiv')

colnames(d) <- c('chr', 'start', 'end', 'name', 'entropy', 'strand', 'beta', 'hmm')

## d$beta <- d$beta/1000 # not needed anylonger, now scores are unmodified

d$coverage <- sapply(strsplit(as.character(gsub('[^0-9;]', '', d$name)), ';'),
       function(x) sum(as.numeric(x)))

## d$m <- beta2m(d$beta)
## d$distance <- d$end - d$start
## d$log10dist <- log10(d$distance)


## reports

wd <- file.path(dirname(colored_entropy), 'reports')
dir.create(wd, showWarnings = FALSE)


## report 1 max entropies   ##############################################################

if (nrow(d) > 10000) {
    set.seed(4)
    
    idx <- sample(1:nrow(d), 10000)

    sampled <- d[idx,]
} else
    sampled <- d


sampled$pU <- (sampled$cov - (sampled$cov * sampled$beta))/sampled$cov
sampled$pM <- 1 - sampled$pU

sampled$max_entropy <- apply(data.frame(pUU = sampled$pU^2,
                                        pUM = sampled$pU * sampled$pM,
                                        pMU = sampled$pU * sampled$pM,
                                        pMM = sampled$pM^2),
                             1, entropy)


sampled$min_entropy_full <- apply(data.frame(pUU = sampled$pU,
                                             pUM = 0,
                                             pMU = 0,
                                             pMM = sampled$pM),
                                  1, entropy)


sampled$min_entropy_beta <- apply(data.frame(pUU = 0,
                                             pUM = abs(sampled$pM - sampled$pU)/2 ,
                                             pMU =  0,
                                             pMM = sampled$pM),
                                  1, entropy)

sampled$min_entropy_beta[sampled$beta >= 0.5] <- apply(data.frame(
    pU = 0,
    pUM = abs(d[sampled$beta >= 0.5,]$pM - d[sampled$beta >= 0.5,]$pU)/2 ,
    pMU =  0,
    pMM = d[sampled$beta >= 0.5,]$pU),
    1,
    entropy)


sampled$lower_bound <- apply(d[,c('min_entropy_full', 'min_entropy_beta')], 1, min)

png(file.path(wd, sprintf('std_bounds.png')))
par(cex.axis = 1.4,
    cex.lab = 1.4,
    cex.main = 1.4,
    cex.sub = 1.4,
    pty = "s",
    mar=c(5.1,4.1,4.1,2.1),
    oma = c(4, 4, 1, 1))

plot(sampled$beta, sampled$entropy, pch = 19, col = ac('black', 0.5),
     xlab = sprintf('%s %s methylation (beta value)', samples_dict[item, 'description'],
                    item),
     ylab = "Shannon's entropy (H)" )
lines(sampled$beta, sampled$max_entropy, col = 'darkred', type = 'p', pch = 4, cex = 1)
lines(sampled$beta, sampled$lower_bound, col = 'darkblue', type = 'p', pch = 4, cex = 1)

legend('topright',
       col = c('black', 'darkred', 'darkblue'),
       pch = c(19, 4, 4),
       legend = c('observation', TeX('$H_{max}$'), TeX('$H_{min}$')))

dev.off()

rm(sampled)

## report 1 max entropies end  ##############################################################

## report 2 start

## report 2 end


###

ls()
str(snakemake)
str(Snakemake)
## save.image('test.RData')

getwd()

## snakemake@source()



date()
sessionInfo()
.libPaths()

## devtools::session_info()


## snakemake success run
write.table(x = rnorm(1),
            file = output)


