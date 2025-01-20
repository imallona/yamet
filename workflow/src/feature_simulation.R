#!/usr/bin/env R
##
## Generates basic, simplistic simulated data to test yamet
## Doesn't write expected truths (@todo)
##
## Run: `Rscript feature_simulation.R`
##
## Started 28 Nov 2024

NCPGS=1000
GRAFT= 10 ## grafted feature cpg length

stopifnot (4 * GRAFT < NCPGS, GRAFT < 100)

sims <- list()
for (chr in c('chrM', 'chrU', 'chr1_1','chr3_1', 'chr_random')) {
    fd = data.frame(start = seq(from = 1, to = NCPGS * 2, by = 2),
                    beta = NA,
                    chr = chr)
    sims[[chr]] <- fd
}

## fully methylated
sims[['chrM']]$beta <- 1
sims[['chrM']]$status <- 'M'

## fully unmethylated
sims[['chrU']]$beta <- 0
sims[['chrU']]$status <- 'U'

## alternate fully meth, fully unmeth
sims[['chr1_1']]$beta <- c(0,1)
sims[['chr1_1']]$status <- '1_1'

## repeat three meth, one unmeth
sims[['chr3_1']]$beta <- c(1,1,1,0)
sims[['chr3_1']]$status <- '3_1'

## Bernoulli trials
set.seed(6574)
sims[['chr_random']]$beta <- sample(c(0,1), NCPGS, replace = TRUE)
sims[['chr_random']]$status <- 'random'


## As for the regions, we graft five $GRAFT-CpG-long regions per chromosome:
##  one fully methylated in position 1;
##  one fully unmethylated in position 101;
##  one alternating in position 201;
##  one with three meth, 1 unmeth in position 301;
##  a random one in position 401.

grafted <- sims
for (chr in names(sims)) {
    grafted[[chr]][1: (1 + GRAFT -1), c('beta', 'status')] <- sims[['chrM']][1: (1+ GRAFT -1),
                                                                             c('beta', 'status')]
    grafted[[chr]][101: (101 + GRAFT -1),c('beta', 'status')] <- sims[['chrU']][  101: (101 + GRAFT -1),
                                                                                c('beta', 'status')]
    grafted[[chr]][201: (201 + GRAFT -1), c('beta', 'status')] <- sims[['chr1_1']][201: (201 + GRAFT -1),
                                                                                   c('beta', 'status')]
    grafted[[chr]][301: (301 + GRAFT -1), c('beta', 'status')] <- sims[['chr3_1']][301: (301 + GRAFT -1),
                                                                                   c('beta', 'status')]
    grafted[[chr]][401: (401 + GRAFT -1), c('beta', 'status')] <- sims[['chr_random']][  401: (401 + GRAFT -1),
                                                                                       c('beta', 'status')]
}

lapply(sims, function(x) table(c(x$chr, x$status)))
lapply(grafted, function(x) table(c(x$chr, x$status)))

## to generate the annotated regions/intervals file, with their true status, we do a bit of a convoluted thing

collapsed <- do.call(rbind.data.frame, grafted)

ref <- collapsed[c('chr','start', 'status')]

sel <- list()

last_status <- ref$status[1]
last_chr =  ref$chr[1]
block_start = last_start = ref$start[1]
block_end <- NCPGS

for (i in 2:nrow(ref)) {
    update = FALSE
    curr_status <- ref$status[i]
    curr_chr <- ref$chr[i]

    cat('.')

    if (curr_chr != last_chr){
        cat('+ ', curr_chr, '\n')
        block_end = NCPGS
        update = TRUE
    } else {
        if (curr_status != last_status) {
            cat('-')
            block_end = ref$start[i] -1
            update = TRUE
        }
    }
    if (i == nrow(ref)) {
        cat('+ ending\n')
        block_end = NCPGS
        update = TRUE
    }
    if (update) {
        sel[[as.character(i)]] <- c(last_chr, last_start, block_end, last_status)
        block_start = ref$start[i]
        last_start = block_start
    }

    last_end <- ref$start[i]    
    last_status <- ref$status[i]    
    last_chr <- ref$chr[i]
    
}

regions <- as.data.frame(do.call(rbind, sel))
colnames(regions) <- c('chr', 'start', 'end', 'status')
regions$start <- as.numeric(regions$start)
regions <- regions[order(regions$chr, regions$start),]

## num methylated reads
collapsed$methylated <- collapsed$beta

## total num reads
collapsed$total <- ifelse(collapsed$beta == 0, yes = 2, no = 1)

stopifnot(all(collapsed$beta == collapsed$methylated/collapsed$total))

## and add sparsity by eliminating every nineth element of the simulated data
collapsed <- collapsed[1:nrow(collapsed) %% 9 != 0, c('chr', 'start', 'methylated', 'total', 'beta')]

## write outputs
write.table(regions, file = 'regions.bed', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(collapsed, file = 'simulations.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(ref, file = 'reference.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

print('Done.')
