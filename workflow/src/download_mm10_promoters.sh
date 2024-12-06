#!/bin/bash
##
## Retrieves the mm10 ensGene promoters, as defined by their TSSs, and sends them to stdout
##
## Izaskun Mallona
## 2nd Dec 2024

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
      'SELECT chrom, txStart, txEnd, name, alignID,  strand 
      FROM mm10.knownGene;' |  awk  '{OFS=FS="\t"; {print $1, $2+1-2000, $3+2000, $4, $5, $6}}'
