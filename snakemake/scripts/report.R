#!/bin/env R

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

