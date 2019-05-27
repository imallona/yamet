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
    library(ggplot2)
    library(RColorBrewer)
    library(gridExtra)
    library(dplyr)
    ## library(pheatmap)
    library(ggExtra)
    ## library(tidyverse)
    library(viridis)
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

    p=ggplot(d, aes(x=beta, y=standardized_entropy, color=entropy)) +
        geom_point() +
        theme(legend.position ='top') + xlab(sprintf('%s methylation (beta)', sample_name)) +
        ylab('standardized entropy (stdH)')
    
    ## marginal density
    p2 <- ggMarginal(p, type="histogram")

    p <- ggplot(d, aes(x=entropy, y=standardized_entropy, color=beta)) +
        geom_point() +
        theme(legend.position ='top')  + xlab('entropy (H)') +
        ylab('standardized entropy (stdH)')

    p3 <- ggMarginal(p, type="histogram")

    p4 <-grid.arrange(p2, p3, ncol=2)
    p4
    ggsave(sprintf(file.path(wd, fn)), plot = p4, width = 10, height=5)
}

## bycolor


add_in_boundary_boolean <- function(d, colors) {

    for (color in colors) {
        d[,color] <- NA
        d[,color] <- grepl(color, as.character(d$hmm))
    }

    d$in_boundary <- apply(d[,as.character(colors)], 1, function(x) sum(x) > 1)
    
    return(d)
}

bycolor1 <- function(d, wd, colors, sample_id ='') {

    entropy_max <- max(d$entropy)
    
    targets <-  c('coverage',
                  'log10dist',
                  'beta',
                  'entropy')

    for (color in colors) {

        curr <- d[d[,color],]
        png(file.path(wd, make.names(sprintf('exploratory_1_%s.png', color))),
            width = 700,
            height = 700)

        plot(curr[,targets],
             pch = 20,
             col = ifelse(curr[,"in_boundary"], 'red', 'black'),
             main = color)
        
        dev.off()

        png(file.path(wd, make.names(sprintf('exploratory_2_%s.png', color))),
            width = 400,
            height = 400)
        
        plot(x = curr$beta,
             y = curr$entropy,
             xlim = c(0,1),
             ylim = c(0, entropy_max),
             pch = 20,
             col = ac(ifelse(curr[,"in_boundary"], 'red', 'black'), 0.5),
             main = color)
        
        dev.off()

        p <- ggplot(curr[,c('beta', 'entropy')], aes(beta, entropy))

        ## h3 <- p + stat_bin2d(bins=25) + scale_fill_gradientn(colours=r) + xlim(-0.1, 1.1) +
        ##     ylim(-0.1, entropy_max + 0.1) + ggtitle(color)
        
        h3 <- p + stat_bin2d(bins=25) + scale_fill_gradientn(colours=terrain.colors(20)) + xlim(-0.1, 1.1) +
            ylim(-0.1, entropy_max + 0.1) + ggtitle(color) + theme(text = element_text(size=20)) +
            xlab('DNA methylation (beta)') + ylab("Shannon's entropy (H)")
                
        ggsave(h3, file = file.path(wd, make.names(sprintf('exploratory_3_%s.png', color))))

    }
}

bycolor2 <- function(d, wd, colors, sample_id ='') {
    
    d$unmethylated <- d$beta < 0.2
    d$methylated <- d$beta >= 0.8
    d$zero_entropy <- d == 0

    non_boundary <- list(median = tapply(d[!d$in_boundary, 'entropy'],
                                         as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                         function(x) quantile(x, probs = 0.5)),
                         lower = tapply(d[!d$in_boundary, 'entropy'],
                                        as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                        function(x) quantile(x, probs = 0.25)),
                         upper = tapply(d[!d$in_boundary, 'entropy'],
                                        as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                        function(x) quantile(x, probs = 0.755)))

    ## sort by median, plot
    sorted <- names(sort(non_boundary$median))

    png(file.path(wd, 'quantiles_entropy_no_boundaries.png'), height = 480, width = 550)

    par(cex.axis = 1.8,
        cex.lab = 1.8,
        cex.main = 1.8,
        cex.sub = 1.8,
        pty = "s",
        mar=c(5.1,4.1,4.1,2.1),
        oma = c(1.2, 0, 0, 0),
        xpd = TRUE)

    plot(non_boundary$lower[sorted], pch = 20,
         ylim = c(0, max(unlist(non_boundary))),
         ylab = "Shannon's entropy (H)" ,
         xaxt = "n",
         xlab = "",
         main = '')

    points(non_boundary$median[sorted], pch = 20, col = 'red')
    points(non_boundary$upper[sorted], pch = 20)

    segments(y0 = non_boundary$lower[sorted],
             y1 = non_boundary$upper[sorted],
             x0 = 1:length(non_boundary$lower),
             x1 = 1:length(non_boundary$lower))

    legend('topright', col = c('black', 'red', 'black'), pch = 20, c('3Q', 'median', '1Q'),
           inset=c(-0.25,0))

    axis(1, at = 1:length(sorted), labels = sorted, las = 2)

    dev.off()

    png(file.path(wd, 'quantiles_methylation_no_boundaries.png'))

    par(cex.axis = 1.8,
        cex.lab = 1.8,
        cex.main = 1.8,
        cex.sub = 1.8,
        pty = "s",
        mar=c(5.1,4.1,4.1,2.1),
        oma = c(1.2, 0, 0, 0))

    plot(non_boundary$lower[sorted], pch = 20,
         ylim = c(0, max(unlist(non_boundary))),
         ylab = 'methylation (beta)',
         xaxt = "n",
         xlab = "",
         main = '')

    points(non_boundary$median[sorted], pch = 20, col = 'red')
    points(non_boundary$upper[sorted], pch = 20)
    segments(y0 = non_boundary$lower[sorted],
             y1 = non_boundary$upper[sorted],
             x0 = 1:length(non_boundary$lower),
             x1 = 1:length(non_boundary$lower))

    axis(1, at = 1:length(sorted), labels = sorted, las = 2)

    dev.off()

}

bycolor3 <- function(d, wd, colors, sample_id ='') {    
    curr <- d[!d$in_boundary,]
    curr$hmm <- as.factor(as.character(curr$hmm))
    ## categorical methylation status
    curr$meth_cat <- 'mid'
    curr$meth_cat[curr$unmethylated] <- 'low'
    curr$meth_cat[curr$methylated] <- 'high'
    curr$meth_cat <- factor(curr$meth_cat, levels = c('low', 'mid', 'high'))

    h <- ggplot(aes(y = entropy, x = hmm, fill = meth_cat), data = curr) +
        geom_boxplot(outlier.alpha = 0.5, color = 'darkblue') +
        scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
        ylab("Shannon's entropy (H)") +
        xlab('chromatin state (chromHMM)') +
        theme_bw() + 
        theme(text = element_text(size = 15),
              axis.text.x = element_text(angle = 90, hjust = 1))


    ggsave(h, file = file.path(wd, 'entropy_no_boundary_boxplot_stratified_by_meth.png'),
           height = 4, width = 8)


}

bycolor4 <- function(d, wd, colors, sample_id ='') {
    curr <- d[!d$in_boundary,]
    curr$hmm <- as.factor(as.character(curr$hmm))

    curr$meth_cat <- 'mid'
    curr$meth_cat[curr$unmethylated] <- 'low'
    curr$meth_cat[curr$methylated] <- 'high'
    curr$meth_cat <- factor(curr$meth_cat, levels = c('low', 'mid', 'high'))

    
    h1 <- ggplot(curr, aes(factor(hmm),fill = meth_cat)) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +    
    geom_bar(stat="count", width = 0.6, position = position_dodge(width=0.7), color = 'black') + 
    scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
    theme_bw() +
    xlab('chromatin state (chromHMM)') +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotation_logticks(sides = 'l') 


    ggsave(h1, file = file.path(wd, 'counts_no_boundary_boxplot_stratified_by_meth.png'),
           height = 4, width = 8)
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

colnames(d)
table(d$hmm)

d$log10dist <- log10(d$end - d$start)



d <- add_in_boundary_boolean(d = d, colors = hmm_colors)
str(d)
table( ac(ifelse(d[,"in_boundary"], 'red', 'black'), 0.5))

print('report 4 untested')

## bycolor1(d = d ,
##          wd = wd, colors = hmm_colors, sample_id ='')

print('report 5 untested')

bycolor2(d = d,
         wd = wd, colors = hmm_colors, sample_id ='')

print('report 6 untested')

bycolor3(d = d,
         wd = wd, colors = hmm_colors, sample_id ='')


print('report 7 untested')

bycolor4(d = d,
         wd = wd, colors = hmm_colors, sample_id ='')




date()
sessionInfo()
.libPaths()

## devtools::session_info()


## snakemake success run
write.table(x = rnorm(1),
            file = output)


