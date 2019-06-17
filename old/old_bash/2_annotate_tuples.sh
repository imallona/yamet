#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 15th november 2018

export HOME=/home/imallona
export TASK="cg_shadows"
export WD="$HOME"/"$TASK"

export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export METILENE="$SOFT"/metilene/metilene_v0.2-7/metilene
export METILENE_INPUT="$SOFT"/metilene/metilene_v0.2-7/metilene_input.pl
export METILENE_OUTPUT="$SOFT"/metilene/metilene_v0.2-7/metilene_output.pl

NTHREADS=12

mkdir -p $WD/data
cd $_


## why some cpgs have a value only for a measurement? @todo check the metilene
# bedtools merge call
# fn=sorjuela_metilene

# rather use the methtuple counts
fn=sorjuela.CG.2.tsv

# the idea is to bedtools annotate these by chromatin colors

# why not imr90? https://www.ncbi.nlm.nih.gov/sra?term=SRX332735

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -e \
      "SELECT chrom, chromStart, chromEnd, name \
    FROM hg19.wgEncodeBroadHmmH1hescHMM" | \
    awk '{ FS="\t"; OFS = "\t";  if (NR> 1) print $0}' > hg19_h1_hmm_colors.bed

$BEDTOOLS sort -i hg19_h1_hmm_colors.bed > foo
mv -f foo hg19_h1_hmm_colors.bed

awk '{OFS=FS="\t"; print "chr"$1,$3,$3+1,$2}' $fn | sed '1d' | $BEDTOOLS sort | \
    $BEDTOOLS closest -a - \
              -b hg19_h1_hmm_colors.bed -t first > left.closest

## check indexing (does it retrieve the cg?)
awk '{OFS=FS="\t"; print "chr"$1,$4-1,$4,$2}' $fn | sed '1d' | $BEDTOOLS sort | \
    $BEDTOOLS closest -a - \
              -b hg19_h1_hmm_colors.bed -t first > right.closest


wc -l $fn *closest

# paste <(sed '1d' $fn | awk '{FS=OFS="\t"; print $1,$2,$3,$4}') <(awk '{print $8}' left.closest ) <(awk '{print $8}' right.closest) > merged.all ## this is wrong

## this is not correct because sorting differs from left.closest, fn nd right closest!!!
