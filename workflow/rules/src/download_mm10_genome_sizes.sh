#!/bin/bash
##
## Retrieves mm10's genome sizes and writes them to the stdout

mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
      'SELECT chrom,size FROM mm10.chromInfo'
