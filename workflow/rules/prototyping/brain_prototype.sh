#!/bin/bash

##
## Processes the brain data - prototype
## https://www.nature.com/articles/s41586-020-03182-8
##
## Started 04 March 2025
##
## Izaskun Mallona
## GPLv3

WD=$HOME/tmp/yamet_brain
REFERENCE=/home/imallona/mm10_cpg.ref.gz


mkdir -p $WD ; cd $WD

## reference generation start ##################################################################

## Kent binaries from https://hgdownload.soe.ucsc.edu/admin/exe/
## parallel no idea where does it come from (!)
export PATH=/home/atchox/ref/kent:/home/atchox/ref/parallel:$PATH

## Chromosomes to process
CHROMOSOMES=($(seq 1 19) X Y)

## Create output directory
mkdir -p out
cd out

## Define the function to process each chromosome
process_chromosome() {
    CHR=$1
    FA="https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/$CHR.fa"

    curl https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/"$FA".gz --output "$FA".gz

    gunzip "$FA".gz

    faSize $FA -detailed >"$FA".sizes

    bedtools makewindows -g "$FA".sizes -w 3 -s 1 |
        bedtools getfasta -fi "$FA" -bed stdin -tab |
        awk '$2 ~ /^CG/ {split($1, a, ":"); split(a[2], b, "-"); print "chr" a[1], b[1], $2}' >"$CHR".mm10.cg.bed

    gzip "$FA"
}

export -f process_chromosome

## Create an empty final output file
touch mm10.ref

## Use parallel to process chromosomes in parallel (one per core)
$SCRIPT_DIR/parallel/src/parallel process_chromosome ::: "${CHROMOSOMES[@]}"

## Combine all the individual chromosome files into the final output
cat *.mm10.cg.bed >>mm10.ref

rm *.cg.bed

## Compress the combined CG sites file
sort -k1,1 -k2,2n mm10.ref | gzip > mm10.ref.gz

rm mm10.ref


# reference generation end ######################################################################
