#!/bin/bash
##
## Retrieves the mm10 ensGene promoters, as defined by their TSSs, and sends them to stdout
##
## Izaskun Mallona
## 2nd Dec 2024

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
      'SELECT chrom, min(txStart), max(txEnd), name2, strand 
       FROM mm10.wgEncodeGencodeBasicVM25
       GROUP BY name2
       ORDER by chrom, min(txStart);' |
      awk '{OFS=FS="\t"; {print $1, $2+1-2000, $2+2000, $4, ".", $5}}' |
      gzip -c >${snakemake_output[0]}
