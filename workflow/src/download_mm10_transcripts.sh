#!/bin/bash
##
## Retrieves the mm10 ensGene transcript and sends them to stdout
## DEPRECATED, use download_mm10_promoters.sh instead
##
## Izaskun Mallona
## 2nd Dec 2024

echo Deprecated

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
      'SELECT chrom, txStart, txEnd, name, alignID,  strand 
      FROM mm10.knownGene;' |  awk  '{OFS=FS="\t"; {print $1, $2, $3, $4, $5, $6}}'
