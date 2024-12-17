#!/bin/bash
##
## Retrieves the mm10 ensGene genes and sends them to stdout
##
## Izaskun Mallona
## 11th Dec 2024

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
      'SELECT chrom, min(txStart), max(txEnd), name2, strand 
       FROM mm10.wgEncodeGencodeBasicVM25
       GROUP BY name2
       ORDER by chrom, min(txStart);' |  awk  '{OFS=FS="\t"; {print $1, $2, $3, $4, ".", $5}}'
