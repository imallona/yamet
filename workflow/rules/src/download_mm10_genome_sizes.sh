#!/bin/bash
##
## Retrieves mm10's genome sizes and writes them to the stdout
##
## Izaskun Mallona

mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
      'SELECT chrom,size FROM mm10.chromInfo' |
      grep -E -w "chr[0-9XY]{1,2}" |
      sort -k1,1 >${snakemake_output[0]}
