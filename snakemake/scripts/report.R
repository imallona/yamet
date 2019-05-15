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


ac  <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    m
}

##' Gets the upper and lower entropy bounds, as used for normalization
##'
##' 
##' @title 
##' @param d and annotated 2-tuples dataframe
##' @param sample a boolean, whether to sample to 10k rows or not (TRUE by default)
##' @return 
##' @author Izaskun Mallona
get_entropy_bounds <- function(d, sample = TRUE) {
    if (nrow(d) > 10000 & sample) {
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
        pUM = abs(sampled[sampled$beta >= 0.5,]$pM - sampled[sampled$beta >= 0.5,]$pU)/2 ,
        pMU =  0,
        pMM = sampled[sampled$beta >= 0.5,]$pU),
        1,
        entropy)

    sampled$lower_bound <- apply(sampled[,c('min_entropy_full', 'min_entropy_beta')], 1, min)
    
    return(sampled)
}

##' Plots the theoretical upper and lower entropy bounds (the ones used for normalization)
##'
##' .. content for \details{} ..
##' @title 
##' @param sampled an input dataframe
##' @param wd working directory
##' @param fn output fn (a PNG)
##' @param sample_name the sample tag to annotate the file
##' @return 
##' @author Izaskun Mallona
plot_entropy_bounds <- function(sampled, wd, fn = 'std_bounds.png', sample_name = '') {
    png(file.path(wd, fn))
    
    par(cex.axis = 1.4,
        cex.lab = 1.4,
        cex.main = 1.4,
        cex.sub = 1.4,
        pty = "s",
        mar=c(5.1,4.1,4.1,2.1),
        oma = c(4, 4, 1, 1))

    plot(sampled$beta, sampled$entropy, pch = 19, col = ac('black', 0.5),
         xlab = sprintf('%s methylation (beta value)', sample_name),
         ylab = "Shannon's entropy (H)" )
    lines(sampled$beta, sampled$max_entropy, col = 'darkred', type = 'p', pch = 4, cex = 1)
    lines(sampled$beta, sampled$lower_bound, col = 'darkblue', type = 'p', pch = 4, cex = 1)

    legend('topright',
           col = c('black', 'darkred', 'darkblue'),
           pch = c(19, 4, 4),
           legend = c('observation', TeX('$H_{max}$'), TeX('$H_{min}$')))

    dev.off()
}


##' Standardize Shannon entropies by the low and upper bounds according to the meth level
##' @title Standardize entropies
##' @param d the methtuple data
##' @return a dataframe with standardized entropy
##' @author Izaskun Mallona
standardize_entropy <- function(d) {
    d <- get_entropy_bounds(d, sample = FALSE)
    
    d$standardized_entropy <- (d$entropy - d$lower_bound) / (d$max_entropy - d$lower_bound)
    d$standardized_entropy[is.na(d$standardized_entropy)] <- 0

    return(d)
}

plot_standardized_entropy <- function(d, wd, fn = 'std_entropies.png', sample_name = "") {
    png(file.path(wd, fn))
    par(cex.axis = 1.4,
        cex.lab = 1.4,
        cex.main = 1.4,
        cex.sub = 1.4,
        pty = "s",
        mar=c(5.1,4.1,4.1,2.1),
        oma = c(4, 4, 1, 1))

    plot(d$beta, d$standardized_entropy, pch = 19, col = ac('black', 0.5),
         xlab = sprintf('%s methylation (beta value)', sample_name),
         ylab = "standardized entropy (stdH)" )

    dev.off()
}


plot_marginals_standardized <- function(d, wd, fn = 'std_vs_entropy.png', sample_name = '')  {

    p=ggplot(d[[item]], aes(x=beta, y=standardized_entropy, color=entropy)) +
        geom_point() +
        theme(legend.position ='top') + xlab(sprintf('%s methylation (beta)', item)) +
        ylab('standardized entropy (stdH)')
    
    ## marginal density
    p2 <- ggMarginal(p, type="histogram")

    p <- ggplot(d[[item]], aes(x=entropy, y=standardized_entropy, color=beta)) +
        geom_point() +
        theme(legend.position ='top')  + xlab('entropy (H)') +
        ylab('standardized entropy (stdH)')

    p3 <- ggMarginal(p, type="histogram")

    p4 <-grid.arrange(p2, p3, ncol=2)
    p4
    ggsave(sprintf(file.path(wd, fn)), plot = p4, width = 10, height=5)
}

########


print('Data load')
fn <- colored_entropy

## d <- data.table::fread(fn) # crashes
d <- read.table(fn)

########

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


print('report 1 reports/std_bounds.png')     ###############################

plot_entropy_bounds(sampled = get_entropy_bounds(d), wd = wd)

print('report 2 reports/std_entropies.png')  ###############################

d <- standardize_entropy(d = d)
plot_standardized_entropy(d = d, wd = wd)

print('report 3 reports/std_vs_entropy.png') ###############################

plot_marginals_standardized(d, wd = wd, fn = 'std_vs_entropy.png', sample_name = '') 

date()
sessionInfo()
.libPaths()

## devtools::session_info()


## snakemake success run
write.table(x = rnorm(1),
            file = output)


