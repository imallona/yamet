#!/bin/bash
##
## Retrieves the mm10 ensGene genes and sends them to stdout
##
## Atreya Choudhury
## 6th Mar 2024

mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
      'SELECT chrom, min(txStart), max(txEnd), name, strand
       FROM hg19.knownGene g1
       GROUP BY name
       ORDER by chrom, min(txStart);' |
      awk '{OFS=FS="\t"; {print $1, $2, $3, $4, ".", $5}}' |
      sort -k1,1 -k2,2n |
      gzip -c >${snakemake_output[0]}
