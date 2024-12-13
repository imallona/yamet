#!/usr/bin/env R
##
## Generates basic, simplistic simulated data to test yamet
## Doesn't write expected truths (@todo)
##
## Run: `Rscript generate_test_inputs.R`
##
## Started 28 Nov 2024

NCPGS=100

sims <- list()
for (chr in c('chrM', 'chrU', 'chr1_1','chr3_1', 'chr_random')) {
    fd = data.frame(start = seq(from = 1, to = NCPGS * 2, by = 2),
                    beta = NA,
                    chr = chr)
    sims[[chr]] <- fd
}

## fully methylated
sims[['chrM']]$beta <- 1

## fully unmethylated
sims[['chrU']]$beta <- 0

## alternate fully meth, fully unmeth
sims[['chr1_1']]$beta <- c(0,1)

## repeat three meth, one unmeth
sims[['chr3_1']]$beta <- c(1,1,1,0)

## Bernoulli trials
set.seed(6574)
sims[['chr_random']]$beta <- sample(c(0,1), NCPGS, replace = TRUE)

## now we generate the reference file, where all simulated CpGs are listed
collapsed <- do.call(rbind.data.frame, sims)
ref <- collapsed[c(3,1)]
ref$end <- ref$start +1

## as well as the regions we'd like to check, just one per chr:
regions <- data.frame(chr = unique(ref$chr),
                      start = min(ref$start),
                      end = max(ref$end) + 1)

## num methylated reads
collapsed$methylated <- collapsed$beta

## total num reads
collapsed$total <- ifelse(collapsed$beta == 0, yes = 2, no = 1)

stopifnot(all(collapsed$beta == collapsed$methylated/collapsed$total))

## and add sparsity by eliminating every nineth element of the simulated data
collapsed <- collapsed[1:nrow(collapsed) %% 9 != 0, c('chr', 'start', 'methylated', 'total', 'beta')]

## write outputs
write.table(regions, file = 'regions.bed', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(collapsed, file = 'simulations.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write.table(ref, file = 'reference.tsv', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

print('Done.')
