#!/usr/bin/env Rscript
##
## Tuple lengths per-sample report
##
## re-written 28 nov 2019
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
})

option_list = list(
    make_option(c("-s", '--sample'),
                type = 'character', default = NULL,
                help = 'sample identifier'),
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

 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat('Running with options:\n')
print(str(opt))

setDTthreads(threads = opt$nthreads)

dir.create(opt$output_path, showWarnings = FALSE)

fns <- list.files(opt$input_path,
                  sprintf('%s.*cov_%s.*entropy_and_meth.bed.gz', opt$sample, opt$coverage))

d <- list()

for (fn in fns) {
    d[[fn]] <- fread(file.path(opt$input_path, fn),
                          header = FALSE)
    colnames(d[[fn]]) <- c('chr', 'start', 'end', 'id', 'entropy', 'strand', 'methylation')
}


## 1 grid of methylation distributions #####################################################

tmp <- lapply(d, function(x)  ggplot(x[,c('methylation', 'entropy')], aes(methylation)) +
                              xlim(c(0,1)) +
                              geom_density())


## dirty way of adding the tuple lengths back
cnames <- gsub('_entropy_and_meth.bed.gz', '', gsub(sprintf('%s_', opt$sample), '', names(tmp)))

for (i in 1:length(tmp)) {
    tmp[[i]] <- tmp[[i]] + labs(title = cnames[[i]])
}

p1 <- grid.arrange(grobs = tmp, ncol = 2,
                   top = opt$sample)

ggsave(plot = p1, filename = file.path(opt$out, sprintf('%s_meth_densities_panel.png',
                                                        opt$sample)),
       device = 'png')

saveRDS(object = tmp, file = file.path(opt$out, sprintf('%s_meth_densities_panel_data.rds',
                                                        opt$sample)))


rm(tmp, cnames, p1)

## 2 grid of entropies distributions ########################################################

tmp <- lapply(d, function(x)  ggplot(x[,c('methylation', 'entropy')], aes(entropy)) +
                              xlim(c(0,2)) + 
                              geom_density())

cnames <- gsub('_entropy_and_meth.bed.gz', '', gsub(sprintf('%s_', opt$sample), '', names(tmp)))

for (i in 1:length(tmp)) {
    tmp[[i]] <- tmp[[i]] + labs(title = cnames[[i]])
}

p2 <- grid.arrange(grobs = tmp, ncol = 2,
                   top = opt$sample)

ggsave(plot = p2, filename = file.path(opt$out, sprintf('%s_entropy_densities_panel.png',
                                                        opt$sample)),
       device = 'png')

saveRDS(object = tmp, file = file.path(opt$out, sprintf('%s_entropy_densities_panel_data.rds',
                                                        opt$sample)))

rm(tmp, cnames, p2)

## grid of meth vs entropy correlations #####################################################

tmp <- lapply(d, function(x)  ggplot(x[,c('methylation', 'entropy')], aes(entropy, methylation)) +
                              xlim(c(0,2)) +
                              geom_point() + 
                              stat_density_2d())

cnames <- gsub('_entropy_and_meth.bed.gz', '', gsub(sprintf('%s_', opt$sample), '', names(tmp)))

for (i in 1:length(tmp)) {
    tmp[[i]] <- tmp[[i]] + labs(title = cnames[[i]])
}

p3 <- grid.arrange(grobs = tmp, ncol = 2,
                   top = opt$sample)

ggsave(plot = p3, filename = file.path(opt$out,
                                       sprintf('%s_meth_vs_entropy_panel.png', opt$sample)),
       device = 'png')

saveRDS(object = tmp, file = file.path(opt$out, sprintf('%s_meth_vs_entropy_data.rds',
                                                        opt$sample)))

rm(tmp, p3, cnames)

## grid of pairwise methylation values correlation ########################################

tmp <- do.call(rbind.data.frame,
               lapply(names(d), function(x)
                   data.frame(methylation = d[[x]]$methylation,
                              position = sprintf('%s:%s(%s)', d[[x]]$chr, d[[x]]$start,
                                                 d[[x]]$strand),
                              sample = x))) %>%
    mutate(sample = gsub('_entropy_and_meth.bed.gz', '', sample),
           sample = gsub(sprintf('%s_', opt$sample), '', sample),
           sample = gsub(sprintf('cov_%s_', opt$coverage), '', sample)) %>%
    group_by(position, sample) %>%
    summarize(Mean = mean(methylation, na.rm=TRUE)) %>%
    spread(key= sample, value = Mean)

p4 <- ggpairs(tmp[,-1], aes(alpha = 0.05),
              title = opt$sample, xlab = 'DNA methylation',
              ylab = 'DNA methylation') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot = p4, filename = file.path(opt$out, sprintf('%s_ggpairs_methylation.png', opt$sample)),
       device = 'png')

saveRDS(object = tmp, file = file.path(opt$out, sprintf('%s_ggpairs_methylation_data.rds',
                                                        opt$sample)))

rm(tmp, p4)

## grid of entropies correlation  ##########################################################

tmp <- do.call(rbind.data.frame,
               lapply(names(d), function(x)
                   data.frame(methylation = d[[x]]$entropy,
                              position = sprintf('%s:%s(%s)', d[[x]]$chr, d[[x]]$start,
                                                 d[[x]]$strand),
                              sample = x))) %>%
    mutate(sample = gsub('_entropy_and_meth.bed.gz', '', sample),
           sample = gsub(sprintf('%s_', opt$sample), '', sample),
           sample = gsub(sprintf('cov_%s_', opt$coverage), '', sample)) %>%
    group_by(position, sample) %>%
    summarize(Mean = mean(methylation, na.rm=TRUE)) %>%
    spread(key= sample, value = Mean)


p5 <- ggpairs(tmp[,-1], aes(alpha = 0.05),
              title = opt$sample, xlab = "Shannon's entropy (H)",
              ylab = "Shannon's entropy (H)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot = p5, filename = file.path(opt$out, sprintf('%s_ggpairs_entropy.png', opt$sample)),
       device = 'png')

saveRDS(object = tmp, file = file.path(opt$out, sprintf('%s_ggpairs_entropy_data.rds',
                                                        opt$sample)))

rm(tmp, p5)
