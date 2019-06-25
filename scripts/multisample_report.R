#!/bin/env R
##
## sample-based meth report
## Izaskun Mallona
##
## Thu May  29 2019

suppressPackageStartupMessages({
    library(data.table)
    ## library('latex2exp')
    ## library(entropy)
    ## library(cramer)
    ## library(Cairo)
    ## library(knitr)
    library(ggplot2)
    library(RColorBrewer)
    ## library(gridExtra)
    ## library(dplyr)
    ## library(pheatmap)
    ## library(ggExtra)
    ## library(tidyverse)
    ## library(viridis)
    library(PMCMRplus) ## posthoc KW
})

    
print(args)

for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

## the argument should be the commaseparated glob of tsvs from run_sample_report
annotated <- strsplit(annotated, ',')


## manual processing start


fns <- list.files('.', pattern = "*integrate.tsv.gz")
samples <- sapply(strsplit(basename(fns), '_'), function(x) return(x[1]))
hmms <-  sapply(strsplit(basename(fns), '_'), function(x) return(x[8]))

d <- list()
for (item in uniqeue(samples)) d[[item]] <- list()
for (i in 1:length(fns)) {
    d[[samples[i]]][[hmms[i]]] <- read.table(fns[i], header = FALSE)
    d[[samples[i]]][[hmms[i]]]$hmm <- hmms[i]
     d[[samples[i]]][[hmms[i]]]$sample <- samples[i]
    
}


## manual processing end

###


date()
sessionInfo()
.libPaths()

## devtools::session_info()


## snakemake success run
write.table(x = rnorm(1),
            file = output)


