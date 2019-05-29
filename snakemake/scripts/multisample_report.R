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

    
## print(args)

## for (i in seq_len(length(args))) {
##     eval(parse(text = args[[i]]))
## }

## ## the argument should be the commaseparated glob of tsvs from run_sample_report
## annotated <- strsplit(annotated, ',')


## manual processing start


fns <- list.files('.', pattern = "*integrate.tsv.gz", recursive = TRUE)
samples <- sapply(strsplit(basename(fns), '_'), function(x) return(x[1]))
hmms <-  sapply(strsplit(basename(fns), '_'), function(x) return(x[8]))

d <- list()
for (item in uniqeue(samples)) d[[item]] <- list()
for (i in 1:length(fns)) {
    d[[samples[i]]][[hmms[i]]] <- read.table(fns[i], header = FALSE)
    d[[samples[i]]][[hmms[i]]]$segmentation <- hmms[i]
    d[[samples[i]]][[hmms[i]]]$sample <- samples[i]
    
}


for (item in names(d)) {
    d[[item]] <- do.call(rbind.data.frame, d[[item]])
}
d <- do.call(rbind.data.frame, d)

dim(d)

set.seed(-2)
sampling <- sample(1:nrow(d), 100000)

sampled <- d[sampling,]


print('parse this from the yaml!')
## whether the HMM and the tissue sort of `match` (same origin)
sampled$consistent <- FALSE
sampled$consistent[sampled$sample == 'ENCFF112TXF' & sampled$segmentation == 'E017'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF957OIM' & sampled$segmentation == 'E118'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF572KNK' & sampled$segmentation == 'E118'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF193RVP' & sampled$segmentation == 'E117'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF845VFH' & sampled$segmentation == 'E117'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF079RGH' & sampled$segmentation == 'E126'] <- TRUE
sampled$consistent[sampled$sample == 'ENCFF119ELB' & sampled$segmentation == 'E126'] <- TRUE

sampled$description <- NA
sampled$description[sampled$sample == 'ENCFF112TXF'] <- 'IMR90'
sampled$description[sampled$sample == 'ENCFF957OIM'] <- 'HepG2'
sampled$description[sampled$sample == 'ENCFF572KNK'] <- 'HepG2'
sampled$description[sampled$sample == 'ENCFF193RVP'] <- 'HeLa-S3'
sampled$description[sampled$sample == 'ENCFF845VFH' ] <-'HeLa-S3'
sampled$description[sampled$sample == 'ENCFF079RGH' ] <- 'GM23248'
sampled$description[sampled$sample == 'ENCFF119ELB' ] <- 'GM23248'

colnames(sampled)[1:5] <- c('entropy', 'beta', 'standardized_entropy', 'in_boundary', 'hmm')

for (item in c('entropy', 'beta', 'standardized_entropy')) {
    sampled[,item] <- as.numeric(as.character(sampled[,item]))
}

for (item in c('segmentation', 'hmm', 'in_boundary')) {
    sampled[,item] <- as.factor(sampled[,item])
}

sampled$in_boundary <- as.logical(sampled$in_boundary)

summary(sampled)

kt <- kruskal.test(standardized_entropy~hmm, data = sampled[!sampled$in_boundary,])

kt

dunn <- kwAllPairsDunnTest(standardized_entropy~hmm, data = sampled[!sampled$in_boundary,])
plot(dunn)

dunn_consistent <- kwAllPairsDunnTest(standardized_entropy~hmm, data = sampled[!sampled$in_boundary & sampled$consistent,])


dunn_unconsistent <- kwAllPairsDunnTest(standardized_entropy~hmm, data = sampled[!sampled$in_boundary & !sampled$consistent,])

par(mfrow = c(2,2))
plot(dunn)
plot(dunn_consistent)
plot(dunn_unconsistent)

dunn_consistent$statistic - dunn_unconsistent$statistic
    
plot(density(na.omit(as.numeric(dunn_consistent$statistic - dunn_unconsistent$statistic))))


p1 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=hmm, y=standardized_entropy, fill=description)) + 
    geom_boxplot()

ggsave('p1.png', plot = p1)
       
p2 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=hmm, y=standardized_entropy, fill=description)) + 
    geom_boxplot()


ggsave('p2.png', plot = p2)

p3 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=hmm, y=standardized_entropy, fill=description)) +
  geom_boxplot()+
  facet_wrap(.~sample, ncol =3)+
  labs(x="sample")+
    theme(axis.text.x=element_text(angle = 90))

ggsave('p3.png', plot = p3)

p4 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=hmm, y=standardized_entropy, fill=description)) +
    geom_boxplot()+
    facet_wrap(~ sample*consistent, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90)) 

ggsave('p4.png', plot = p4)

p5 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=sample, y=standardized_entropy, fill=consistent)) +
    geom_boxplot()+
    facet_wrap(~ hmm, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90))

ggsave('p5.png', plot = p5)


p6 <- ggplot(data = sampled[!sampled$in_boundary,],
       aes(x=description, y=standardized_entropy, fill=consistent)) +
    geom_boxplot()+
    facet_wrap(~ hmm, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90)) 

ggsave('p6.png', plot = p6)

## do consistent elements show higher variability?
                                                  
## manual processing end

###


## date()
## sessionInfo()
## .libPaths()

## ## devtools::session_info()


## ## snakemake success run
## write.table(x = rnorm(1),
##             file = output)


