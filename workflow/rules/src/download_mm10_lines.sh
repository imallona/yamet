#!/bin/bash
##
## Downloads RepeatMasker (rmsk) for mm10 and extracts only SINE elements
##

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
  'SELECT genoName, genoStart, genoEnd, repName, repClass, repFamily
     FROM mm10.rmsk
     WHERE repClass = "LINE";' |
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' |
  bedtools merge -i - |
  gzip -c > ${snakemake_output[0]}
