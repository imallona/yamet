#!/bin/bash
##
## Merges transcripts from gencocde 19 (hg19) into nonoverlapping records
##
## 6th Mar 2024

# mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
#       'SELECT chrom, min(txStart), max(txEnd), name, strand
#        FROM hg19.knownGene g1
#        GROUP BY name
#        ORDER by chrom, min(txStart);' |
#       awk '{OFS=FS="\t"; {print $1, $2, $3, $4, ".", $5}}' |
#       sort -k1,1 -k2,2n |
#       bedtools merge -i - |
#       gzip -c >${snakemake_output[0]}

## remove superlong genes / telomeric genes (longer than 100k nt)
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -s -e \
  'SELECT chrom, MIN(txStart), MAX(txEnd), name2, strand 
   FROM hg19.wgEncodeGencodeBasicV19
   GROUP BY name2
   HAVING MAX(txEnd) - MIN(txStart) < 100000
   ORDER BY chrom, MIN(txStart);' |
    awk '{OFS=FS="\t"; {print $1, $2, $3, $4, ".", $5}}' |
    sort -k1,1 -k2,2n |
    bedtools merge -i - |
    gzip -c >${snakemake_output[0]}
