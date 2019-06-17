#!/bin/env R
##
## sample-based meth report
## Izaskun Mallona
##
## Thu May  29 2019

stop('get those genes that are intermediate methylation only, avoid extremes!')

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
    library(nnet) # multinomial
    library(iCOBRA)
    library(caret)
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
for (item in unique(samples)) d[[item]] <- list()
for (i in 1:length(fns)) {
    d[[samples[i]]][[hmms[i]]] <- read.table(fns[i], header = TRUE)
    d[[samples[i]]][[hmms[i]]]$segmentation <- hmms[i]
    d[[samples[i]]][[hmms[i]]]$sample <- samples[i]
    
}


for (item in names(d)) {
    d[[item]] <- do.call(rbind.data.frame, d[[item]])
}
d <- do.call(rbind.data.frame, d)

## colnames(d)[1:5] <- c('entropy', 'beta', 'standardized_entropy', 'in_boundary', 'hmm')


dim(d)

set.seed(-2)
sampling <- sample(1:nrow(d[!d$in_boundary,]), 100000)

sampled <- d[!d$in_boundary,][sampling,]

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
    sampled[,item] <- as.factor(as.character(sampled[,item]))
}

## sampled$in_boundary <- as.logical(sampled$in_boundary)

summary(sampled)

kt <- kruskal.test(standardized_entropy~hmm, data = sampled)

kt

dunn <- kwAllPairsDunnTest(standardized_entropy~hmm, data = sampled)
## png(file.path(wd, 'dunn_std_entropy.pdf'))
## plot(dunn)
## dev.off()

dunn_consistent <- kwAllPairsDunnTest(standardized_entropy~hmm,
                                      data = sampled[sampled$consistent,])

dunn_unconsistent <- kwAllPairsDunnTest(standardized_entropy~hmm,
                                        data = sampled[ !sampled$consistent,])


png(file.path(wd, 'dunn_std_entropy.png'), width = 600, height = 600)
par(mfrow = c(2,2))
plot(dunn, main = 'overall', las = 2)
plot(dunn_consistent, main = 'match', las = 2)
plot(dunn_unconsistent, main = 'discordant', las = 2)
dev.off()


dunn_consistent$statistic - dunn_unconsistent$statistic
    
plot(density(na.omit(as.numeric(dunn_consistent$statistic - dunn_unconsistent$statistic))))

p1 <- ggplot(data = sampled,
       aes(x=hmm, y=standardized_entropy, fill=description)) + 
    geom_boxplot()

ggsave('p1.png', plot = p1)
       
p2 <- ggplot(data = sampled,
       aes(x=hmm, y=standardized_entropy, fill=description)) + 
    geom_boxplot()


ggsave('p2.png', plot = p2)

p3 <- ggplot(data = sampled,
       aes(x=hmm, y=standardized_entropy, fill=description)) +
  geom_boxplot()+
  facet_wrap(.~sample, ncol =3)+
  labs(x="sample")+
    theme(axis.text.x=element_text(angle = 90))

ggsave('p3.png', plot = p3)

p4 <- ggplot(data = sampled,
       aes(x=hmm, y=standardized_entropy, fill=description)) +
    geom_boxplot()+
    facet_wrap(~ sample*consistent, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90)) 

ggsave('p4.png', plot = p4)

p5 <- ggplot(data = sampled,
       aes(x=sample, y=standardized_entropy, fill=consistent)) +
    geom_boxplot()+
    facet_wrap(~ hmm, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90))

ggsave('p5.png', plot = p5)


p6 <- ggplot(data = sampled,
       aes(x=description, y=standardized_entropy, fill=consistent)) +
    geom_boxplot()+
    facet_wrap(~ hmm, ncol = 3) +
    theme(axis.text.x = element_text(angle = 90)) 

ggsave('p6.png', plot = p6)


## which does best, std entropy, entropy, beta, or mixes?

## train and test setup: 75-25, and the incoherent as extra test

train <- list()
test <- list()
tmp <- list()

tmp$coherent <- sampled[sampled$consistent,]
test$incoherent <- sampled[!sampled$consistent,]

for (item in names(tmp)) tmp[[item]]$hmm <- as.factor(as.character(tmp[[item]]$hmm))

smp_size <- floor(0.75 * nrow(tmp$coherent))

set.seed(555)
train_ind <- sample(seq_len(nrow(tmp$coherent)), size = smp_size)

train$coherent <- tmp$coherent[train_ind, ]
test$coherent <- tmp$coherent[-train_ind, ]

# model 

fit1c <- multinom(as.factor(hmm) ~ as.numeric(standardized_entropy),
                  data = train$coherent)

fit2c <- multinom(as.factor(hmm) ~ as.numeric(entropy),
                  data = train$coherent)

fit3c <- multinom(as.factor(hmm) ~ as.numeric(beta),
                  data = train$coherent)

fit4c <- multinom(as.factor(hmm) ~ as.numeric(beta) + as.numeric(entropy),
                  data = train$coherent)

fit5c <- multinom(as.factor(hmm) ~ as.numeric(beta) + as.numeric(standardized_entropy),
                  data = train$coherent)


rank(c(fit1c$AIC, fit2c$AIC, fit3c$AIC,fit4c$AIC))

pvals <- list()

for (item in c('fit1c', 'fit2c',  'fit3c', 'fit4c', 'fit5c')) {
    curr <- get(item)
    z <- summary(curr)$coefficients/summary(curr)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1)) * 2

    pvals[[item]] <-  (1 - pnorm(abs(z), 0, 1)) * 2
}

confmats <- list()
for (item in c('fit1c', 'fit2c',  'fit3c', 'fit4c', 'fit5c')) {
    
    confmats[[item]] <- caret::confusionMatrix(data = predict(get(item), test$coherent),
                                               reference = test$coherent$hmm)
}

sapply(confmats, function(x) print(x$overall))

## what if predicting wrong hmms?

confmats_control <- list()
for (item in c('fit1c', 'fit2c',  'fit3c', 'fit4c', 'fit5c')) {
    
    confmats_control[[item]] <- caret::confusionMatrix(data = predict(get(item), test$incoherent),
                                               reference = test$incoherent$hmm)
}

sapply(confmats_control, function(x) print(x$overall))


# Not a big difference... predictive power is low and beta works better!


train$coherent$is_enh <- grepl('Enh', train$coherent$hmm )
train$coherent$is_tss <- grepl('Tss', train$coherent$hmm )
train$coherent$is_biv <- grepl('Biv', train$coherent$hmm )
train$coherent$is_tx <- grepl('Tx', train$coherent$hmm )


test$coherent$is_enh <- grepl('Enh', test$coherent$hmm )
test$coherent$is_tss <- grepl('Tss', test$coherent$hmm )
test$coherent$is_biv <- grepl('Biv', test$coherent$hmm )
test$coherent$is_tx <- grepl('Tx', test$coherent$hmm )


test$incoherent$is_enh <- grepl('Enh', test$incoherent$hmm )
test$incoherent$is_tss <- grepl('Tss', test$incoherent$hmm )
test$incoherent$is_biv <- grepl('Biv', test$incoherent$hmm )
test$incoherent$is_tx <- grepl('Tx', test$incoherent$hmm )


## now get's better
for (item in c('is_enh', 'is_tss', 'is_biv', 'is_tx')){
    mnom <- multinom(as.factor(train$coherent[,item]) ~ as.numeric(standardized_entropy),
                     data = train$coherent)
    print(item)
    print(caret::confusionMatrix(data = predict(mnom, test$coherent),
                           reference = as.factor(test$coherent[,item]))$overall)
}

## but does not hold for the incoherent scenario
for (item in c('is_enh', 'is_tss', 'is_biv', 'is_tx')){
    mnom <- multinom(as.factor(train$coherent[,item]) ~ as.numeric(standardized_entropy),
                     data = train$coherent)
    print(item)
    print(caret::confusionMatrix(data = predict(mnom, test$incoherent),
                           reference = as.factor(test$incoherent[,item]))$overall)
}



## what if using dnameth? excellent as well, this might not make much sense
for (item in c('is_enh', 'is_tss', 'is_biv', 'is_tx')){
    mnom <- multinom(as.factor(train$coherent[,item]) ~ as.numeric(beta),
                     data = train$coherent)
    print(item)
    print(caret::confusionMatrix(data = predict(mnom, test$coherent),
                           reference = as.factor(test$coherent[,item]))$overall)
}



## that's weird Fri May 31 14:34:01 CEST 2019



## predicted_probabilities

## fitted(fit4c))
## fit4c_cobra <- COBRAData(score = fitted(fit4c), truth = sampled$hmm)



## library(mlogit)

## d1 <- mlogit.data(sampled[sampled$in_bounday,],
##                  shape = "wide", choice = "hmm", varying = c(1:3,8))

## m <- mlogit(depvar ~ ic + oc | 0, H)
## summary(m)

## f1 <- mlogit(hmm ~ standardized_entropy, data = sampled)



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


