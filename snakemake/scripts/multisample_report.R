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
    library(PMCMRplus) ## posthoc KW
})

    
print(args)

for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

## the argument should be the commaseparated glob of tsvs from run_sample_report
annotated <- strsplit(annotated, ',')

###


date()
sessionInfo()
.libPaths()

## devtools::session_info()


## snakemake success run
write.table(x = rnorm(1),
            file = output)


