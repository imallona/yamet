#!/usr/bin/env Rscript
##
## Tuple lengths all samples report
##
## 09 dec 2019
## Izaskun Mallona

suppressPackageStartupMessages({
    library("optparse")
    library('ggplot2')
    library(data.table)
    library(tidyr)
    library(dplyr)
    library(gridExtra)
    ## library(viridis)
    library(GGally)
    library(gtools)
})

option_list = list(
    make_option(c('-c', '--coverage'),
                type = 'numeric', default = 10,
                help = 'min coverage (from file name)'),
    make_option(c("-i", "--input_path"),
                type = "character", default = NULL, 
                help = "path to sample's directory containing inputs", metavar = "character"),
    make_option(c("-o", "--output_path"), type = "character",
                default = NULL, 
                help = "path to output directory", metavar = "character"),
    make_option(c('-t', '--nthreads'),
                type = 'numeric',
                default = 1,
                help = 'number of cores (for datatables)'))


## opt <- list(coverage = 20,
##             input_path = '.',
##             output_path = 'test',
##             nthreads = 10)

            
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat('Running with options:\n')
print(str(opt))

setDTthreads(threads = opt$nthreads)

dir.create(opt$output_path, showWarnings = FALSE)

h_files <- list.files(opt$input_path, sprintf('.*ggpairs_entropy_data.rds'), recursive = TRUE)

h <- list()
h <- lapply(h_files, readRDS)
names(h) <- h_files

for (i in 1:length(h)) {
    h[[i]]$sample <- gsub('_ggpairs_entropy_data.rds', '',  basename(names(h)[i]))
}


h <- gather(do.call(rbind.data.frame, h), tuple_length, entropy, CG_2:CG_9)

p1 <- ggplot(h, aes(x=tuple_length, y=entropy, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(plot = p1, filename = file.path(opt$out,
                                       sprintf('entropy_vs_tuplelength.png')),
       device = 'png')

rm(p1)

p2 <- ggplot(h, aes(x=tuple_length, y=entropy, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 
    

ggsave(plot = p2, filename = file.path(opt$out,
                                       sprintf('entropy_vs_tuplelength_log10.png')),
       device = 'png')

rm(p2)

p3 <- ggplot(h, aes(x=sample, y=entropy)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~tuple_length) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 


ggsave(plot = p3, filename = file.path(opt$out,
                                       sprintf('entropy_vs_tuplelength_sample_log10.png')),
       device = 'png')

rm(p3)

## adding DNA meth data

m_files <- list.files(opt$input_path, sprintf('.*ggpairs_methylation_data.rds'), recursive = TRUE)

m <- list()
m <- lapply(m_files, readRDS)
names(m) <- m_files

for (i in 1:length(m)) {
    m[[i]]$sample <- gsub('_ggpairs_methylation_data.rds', '',  basename(names(m)[i]))
}


m <- gather(do.call(rbind.data.frame, m), tuple_length, methylation, CG_2:CG_9)



p1 <- ggplot(m, aes(x=tuple_length, y=methylation, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(plot = p1, filename = file.path(opt$out,
                                       sprintf('methylation_vs_tuplelength.png')),
       device = 'png')

rm(p1)

p2 <- ggplot(m, aes(x=tuple_length, y=methylation, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 
    

ggsave(plot = p2, filename = file.path(opt$out,
                                       sprintf('methylation_vs_tuplelength_log10.png')),
       device = 'png')

rm(p2)

p3 <- ggplot(m, aes(x=sample, y=methylation)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~tuple_length) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 


ggsave(plot = p3, filename = file.path(opt$out,
                                       sprintf('methylation_vs_tuplelength_sample_log10.png')),
       device = 'png')

rm(p3)


## integrating meth and entropy

colnames(m)[4] <- 'value'
m$variable <- 'methylation'
colnames(h)[4] <- 'value'
h$variable <- 'entropy'
## m$quantile <- quantcut(m$value, q = 4)

m$category <- ifelse(is.na(m$value), NA,
              ifelse(m$value >= 0.75, 'x>=0.75',
              ifelse(m$value >= 0.5, '0.5<x<0.75',
              ifelse(m$value > 0.25, '0.25<x<0.5',
                     'x<0.25'))))

## for entropy
h$category <- ifelse(is.na(h$value), NA,
              ifelse(h$value >= 2.5, 'x>=2.5',
              ifelse(h$value >= 2, '2<x<2.5',
              ifelse(h$value > 0.25, '0.25<x<2',
                     'x<0.25'))))

## still, is the coupling of meth category and entropy category what matters! and sample
m$id <- sprintf('%s_%s_%s', m$position, m$sample, m$tuple_length)
h$id <- sprintf('%s_%s_%s', h$position, h$sample, h$tuple_length)

## ## subset a chromosome (18)
##  print('beware  the head!')
## d <- inner_join(head(m, 1e5)  %>% filter(grepl("^18:", position)),
##                head(h, 1e5) %>% filter(grepl("^18:", position)),
##                by = 'id', suffix = c('.m', '.h'))


## p1 <- ggplot(d, aes(x=category.m, y=value.m, fill = sample.m )) +
##     geom_boxplot(position=position_dodge(1)) +
##     facet_wrap(~tuple_length.m+ category.h) +
##     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


## p1


