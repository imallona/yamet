#!/bin/bash
##
## Parses meth cytosine reports from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
## Single-cell Multi-omics Sequencing and Analyses of Human Colorectal Cancer
## for use with yamet
##
## Atreya Choudhury, Izaskun Mallona
## Started 20thDec 2024

gunzip -c ${snakemake_params[raw]}/${snakemake_wildcards[file]} |
    grep "CpG$" |
    awk '
            {OFS="\t";} {
                if ($4 == "-") $2 = $2 - 2;
                else if ($4 == "+") $2 = $2 - 1;
                print $1, $2, $2+1, "", "", $4, $6, $5;
            }
        ' |
    bedtools merge -c 7,8 -o sum |
    awk '
            {OFS=FS="\t";
                bin = ($4/$5 > 0.1) ? 1 : 0;
                print $1,$2,$4,$5,bin}
        ' |
    sort -k1,1 -k2,2n |
    gzip -c >${snakemake_output[0]}
