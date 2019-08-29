#!/bin/bash
##
## Fetches hg38 refGene coordinates, merges the overlapping/bookended ones
## and reports (bed) their promoters
##
## 29 aug 2019
## Izaskun Mallona

BEDTOOLS=~/soft/bedtools/bedtools2/bin/bedtools
## ~/soft/bedtools/bedtools2/bin/bedtools --version
## bedtools v2.27.1

cd ~/src/meth_entropies/runs/glios

mkdir data
cd $_

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A \
      -B -N \
      -e "select 
           chrom, 
           txStart,
           txEnd,
           name,
           strand           
          from 
           hg38.refGene;" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,"0",$5}' | $BEDTOOLS sort > refGene_hg38.bed

mysql --user=genome --host=genome-mysql.cse.ucsc.edu \
      -A -e "select chrom, size from hg38.chromInfo" > hg38.genome


## bedtools expanding overlapping transcripts and getting the upstream region,
## taking into account the strand

$BEDTOOLS merge -s -i refGene_hg38.bed -c 4,6 -o distinct |  \
    $BEDTOOLS slop -i - -g hg38.genome -l 1000 -r 500 -s > hg38_nonoverlapping_promoters.bed
