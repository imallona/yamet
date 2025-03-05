#!/bin/bash
##
## Retrieves only meth cytosine accessors from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
## Single-cell Multi-omics Sequencing and Analyses of Human Colorectal Cancer
##
## Started 16th Feb 2025

esearch -db sra -query PRJNA382695 |
  efetch -format runinfo |
  cut -f11,13,14,15,30 -d"," |
  grep Bis |
  cut -f5 -d"," >"${snakemake_output[gsm]}"
