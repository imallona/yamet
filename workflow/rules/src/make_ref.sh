#!/bin/bash
##
# Combines different C loci coordinate/reference files, from different chromosomes, into a genome reference file
##
## Atreya Choudhury, Izaskun Mallona
## Started 26th Feb 2025

## Each per-chromosome .ref file is already sorted by position.
## Merge-sort (-m) combines them without loading all rows into memory at once.
## glob expansion is lexicographic, which matches the sort key (chr1 < chr10 < chr2 ...),
## so no needs to resort
sort -m -k1,1 -k2,2n ${snakemake_params[base]}/*.ref | gzip >${snakemake_output[0]}
