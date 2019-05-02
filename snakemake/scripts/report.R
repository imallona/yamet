#!/bin/env R
##
## Izaskun Mallona
##
## Thu May  2 2019

suppressPackageStartupMessages({
    library(data.table)
    library('latex2exp')
    library(entropy)
    library(cramer)
    library(Cairo)
    library(knitr)
    library(ggplot2)
    library(RColorBrewer)
    library(gridExtra)
    library(data.table)
    library(dplyr)
    library(pheatmap)
    library(ggExtra)
    library(tidyverse)
    library(viridis)
})


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

d <- list()

d$entropy <- fread(file.path(snakemake@widlcards$sample,
                             sprintf('%s_cov_%s_entropy_hmm_%s.bed.gz',
                                     snakemake@widlcards$sample,
                                     snakemake@wildcards$cov,
                                     snakemake@wildcards$hmm)))


d$meth <- fread(file.path(snakemake@widlcards$sample,
                          sprintf('%s_cov_%s_methylation_hmm_%s.bed.gz',
                                  snakemake@widlcards$sample,
                                  snakemake@wildcards$cov,
                                  snakemake@wildcards$hmm)))


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

d$merged <- data.frame(pos1 = d$entropy$V2,
                           pos2 =  d$entropy$V3,
                           name = d$entropy$V4,
                           entropy = d$entropy$V5,
                           strand = d$entropy$V6,
                           hmm = d$entropy$V7,
                           beta = d$meth$V5)

d$merged$beta <- d$merged$beta/1000

d$merged$coverage <- sapply(strsplit(as.character(gsub('[^0-9;]', '', d$merged$name)), ';'),
       function(x) sum(as.numeric(x)))

d$merged$m <- beta2m(d$merged$beta)
d$merged$distance <- d$merged$pos2 - d$merged$pos1
d$merged$log10dist <- log10(d$merged$distance)



###

ls()
str(snakemake)
str(Snakemake)
## save.image('test.RData')

getwd()

## snakemake@source()

fn <-  file.path(snakemake@wildcards$sample,
                             sprintf('%s_cov_%s_hmm_%s.test',
                                     snakemake@wildcards$sample,
                                     snakemake@wildcards$cov,
                                     snakemake@wildcards$hmm))
print(fn)

write.table(x = rnorm(1),
            file = fn)

date()
sessionInfo()
.libPaths()

## devtools::session_info()

