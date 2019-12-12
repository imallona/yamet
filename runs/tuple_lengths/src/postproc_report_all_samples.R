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
    library(pheatmap)
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


h <- gather(do.call(rbind.data.frame, h), tuple_length, entropy, CG_2:CG_9, na.rm = TRUE)

table(is.na(h$entropy))


p1 <- ggplot(h, aes(x=tuple_length, y=entropy, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(plot = p1, filename = file.path(opt$output_path,
                                       sprintf('entropy_vs_tuplelength.png')),
       device = 'png')

rm(p1)

p2 <- ggplot(h, aes(x=tuple_length, y=entropy, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 
    

ggsave(plot = p2, filename = file.path(opt$output_path,
                                       sprintf('entropy_vs_tuplelength_log10.png')),
       device = 'png')

rm(p2)

p3 <- ggplot(h, aes(x=sample, y=entropy)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~tuple_length) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 


ggsave(plot = p3, filename = file.path(opt$output_path,
                                       sprintf('entropy_vs_tuplelength_sample_log10.png')),
       device = 'png')

rm(p3)

## similarly with DNA meth data

m_files <- list.files(opt$input_path, sprintf('.*ggpairs_methylation_data.rds'), recursive = TRUE)

m <- list()
m <- lapply(m_files, readRDS)
names(m) <- m_files

for (i in 1:length(m)) {
    m[[i]]$sample <- gsub('_ggpairs_methylation_data.rds', '',  basename(names(m)[i]))
}


m <- gather(do.call(rbind.data.frame, m), tuple_length, methylation, CG_2:CG_9, na.rm = TRUE)

p1 <- ggplot(m, aes(x=tuple_length, y=methylation, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(plot = p1, filename = file.path(opt$output_path,
                                       sprintf('methylation_vs_tuplelength.png')),
       device = 'png')

rm(p1)

p2 <- ggplot(m, aes(x=tuple_length, y=methylation, fill = sample)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 
    

ggsave(plot = p2, filename = file.path(opt$output_path,
                                       sprintf('methylation_vs_tuplelength_log10.png')),
       device = 'png')

rm(p2)

p3 <- ggplot(m, aes(x=sample, y=methylation)) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~tuple_length) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(trans='log10') 


ggsave(plot = p3, filename = file.path(opt$output_path,
                                       sprintf('methylation_vs_tuplelength_sample_log10.png')),
       device = 'png')

rm(p3)


## integrating meth and entropy

m$m_cat <- ifelse(is.na(m$methylation), NA,
           ifelse(m$methylation >= 0.75, 'm>=0.75',
           ifelse(m$methylation >= 0.5, '0.5<m<0.75',
           ifelse(m$methylation >= 0.25, '0.25<m<0.5',
                  'm<0.25'))))

## for entropy
h$h_cat <- ifelse(is.na(h$entropy), NA,
           ifelse(h$entropy >= 1, 'h>=1',
           ifelse(h$entropy >= 0.5, '0.5<h<1',
           ifelse(h$entropy >= 0.1, '0.1<h<0.5',
                  'h<0.1'))))


stopifnot(all(sprintf('%s_%s_%s', m$position, m$sample, m$tuple_length) ==
              sprintf('%s_%s_%s', h$position, h$sample, h$tuple_length)))
## so colpasting is feasible

## downsample

set.seed(667)
sampled <- bind_cols(m, h) %>% sample_frac(0.1, replace = FALSE)
dim(sampled)

p1 <- ggplot(sampled,
             aes(x=m_cat, y=entropy, fill = tuple_length )) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


p1

ggsave(plot = p1, filename = file.path(opt$output_path,
                                       sprintf('integrative_p1.png')),
       device = 'png')

rm(p1)


p2 <- ggplot(sampled,
             aes(x=m_cat, y=entropy, fill = sample )) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~tuple_length) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


ggsave(plot = p2, filename = file.path(opt$output_path,
                                       sprintf('integrative_p2.png')),
       device = 'png', width =300, height = 100)

rm(p2)



p3 <- ggplot(sampled,
             aes(x=h_cat, y=methylation, fill = tuple_length )) +
    geom_boxplot(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


ggsave(plot = p3, filename = file.path(opt$output_path,
                                       sprintf('integrative_p3.png')),
       device = 'png')

rm(p3)



p4 <- ggplot(sampled,
             aes(x=methylation, y=entropy, col = tuple_length)) +
    geom_point(position=position_dodge(1)) +
    facet_wrap(~sample) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

ggsave(plot = p4, filename = file.path(opt$output_path, sprintf('integrative_p4.png')), device = 'png')

rm(p4)

## finally, for each coordinate, plot the correlation between entropy and tuple length,
## and  between meth and tuple length

## we take 10 positions only which are top-represented within a sample

top_intra <- as.data.frame(table(bind_cols(m, h)$position, bind_cols(m, h)$sample))
top_intra <- top_intra[top_intra$Freq == max(top_intra$Freq),]

## intra <- bind_cols(m, h) %>% filter(position %in% head(as.character(top_intra$Var1, 10)))
## intra <- bind_cols(m, h) %>% filter(position == "18:10454771(+)")

set.seed(1)
idx <- m$position %in% sample(as.character(top_intra$Var1), 10)
intra <- bind_cols(m, h)[idx,]

saveRDS(intra, file = file.path(opt$output_path, 'intra_data.rds'))

p5 <- ggplot(intra,
             aes(x = tuple_length, y = entropy, col = sample)) +
    geom_point(position=position_dodge(1)) +
    facet_wrap(~position) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


p6 <- ggplot(intra,
             aes(x = tuple_length, y = methylation, col = sample)) +
    geom_point(position=position_dodge(1)) +
    facet_wrap(~position) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

## similarly with the top 10 positions represented across samples

top_inter <- as.data.frame(table(bind_cols(m, h)$position))
top_inter <- top_inter[top_inter$Freq == max(top_inter$Freq),]

set.seed(3)
idx <- m$position %in% sample(as.character(top_inter$Var1), 10)
inter <- bind_cols(m, h)[idx,]

saveRDS(inter, file = file.path(opt$output_path, 'inter_data.rds'))

p7 <- ggplot(inter,
             aes(x = tuple_length, y = entropy, col = sample)) +
    geom_point(position=position_dodge(1)) +
    facet_wrap(~position) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

p8 <- ggplot(inter,
             aes(x = tuple_length, y = methylation, col = sample)) +
    geom_point(position=position_dodge(1)) +
    facet_wrap(~position) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


# @todo finetune width and height
ggsave(plot = p5, filename = file.path(opt$output_path, sprintf('integrative_p5.png')),
       device = 'png')
ggsave(plot = p6, filename = file.path(opt$output_path, sprintf('integrative_p6.png')),
       device = 'png')
ggsave(plot = p7, filename = file.path(opt$output_path, sprintf('integrative_p7.png')),
       device = 'png')
ggsave(plot = p8, filename = file.path(opt$output_path, sprintf('integrative_p8.png')),
       device = 'png')

rm(top_inter, top_intra, p5, p6, p7, p8)
