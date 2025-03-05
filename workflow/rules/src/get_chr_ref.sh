#!/bin/bash
##
## Retrieves assembly for one chromosome
##
## Started 26th Feb 2025

curl "${snakemake_params[base]}${snakemake_params[fa]}".gz --output ${snakemake_params[fa]}.gz

gunzip ${snakemake_params[fa]}.gz

faSize ${snakemake_params[fa]} -detailed >${snakemake_params[fa]}.sizes

if [[ ${snakemake_wildcards[meth_pat]} == "CG" ]]; then
    regex='^CG'
elif [[ ${snakemake_wildcards[meth_pat]} == "CHH" ]]; then
    regex='^C[ATC][ATC]'
else
    echo "Invalid context. Use 'CG' or 'CHH'."
    exit 1
fi

# make windows and filter for CG sites
bedtools makewindows -g ${snakemake_params[fa]}.sizes -w 3 -s 1 |
    bedtools getfasta -fi ${snakemake_params[fa]} -bed stdin -tab |
    awk -v pattern="$regex" '
        $2 ~ pattern {split($1, a, ":");
        split(a[2], b, "-");
        print "chr" a[1], b[1], $2}' >${snakemake_output[0]}

rm ${snakemake_params[fa]}*
