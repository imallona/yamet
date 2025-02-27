#!/bin/bash
##
## Retrieves assembly for one chromosome
##
## Started 26th Feb 2025

curl "${snakemake_params[base]}${snakemake_params[fa]}".gz --output ${snakemake_params[fa]}.gz

gunzip ${snakemake_params[fa]}.gz
faSize ${snakemake_params[fa]} -detailed >${snakemake_params[fa]}.sizes

bedtools makewindows -g ${snakemake_params[fa]}.sizes -w 3 -s 1 |
    bedtools getfasta -fi ${snakemake_params[fa]} -bed stdin -tab |
    awk '$2 ~ /^CG/ {split($1, a, ":"); split(a[2], b, "-"); print "chr" a[1], b[1], $2}' >${snakemake_output[0]}

rm ${snakemake_params[fa]}*
