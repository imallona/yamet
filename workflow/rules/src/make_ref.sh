#!/bin/bash
##
## Combines different reference files, from different chromosomes, into a genome reference file
##
## Atreya Choudhury
## Started 26th Feb 2025

cat /dev/null >ref.tmp

## Combine all the individual chromosome reference files into the final output
cat ${snakemake_params[base]}/*.ref >>ref.tmp

## Sort and compress the combined reference file
sort -k1,1 -k2,2n ref.tmp | gzip >${snakemake_output[0]}

rm ref.tmp
