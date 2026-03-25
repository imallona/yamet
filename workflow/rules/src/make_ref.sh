#!/bin/bash
##
# Combines different C loci coordinate/reference files, from different chromosomes, into a genome reference file
##
## Atreya Choudhury, Izaskun Mallona
## Started 26th Feb 2025

## Each per-chromosome .ref file is already sorted by position.
## Merge-sort (-m) combines them without loading all rows into memory at once.
## Use snakemake_input directly so only the chromosomes requested by the rule are merged
## (avoids picking up stale .ref files from previous full-genome runs).
sort -m -k1,1 -k2,2n ${snakemake_input[@]} | gzip >${snakemake_output[0]}
