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
## parallelization procedure largely from Atreya (found the script somewhere his home)

## Kent binaries from https://hgdownload.soe.ucsc.edu/admin/exe/
## parallel no idea where does it come from (!)
export PATH=/home/atchox/ref/kent:/home/atchox/ref/parallel/src:$PATH

## Chromosomes to process
CHROMOSOMES=($(seq 1 19) X Y)

## Define the function to process each chromosome
process_chromosome() {
    CHR=$1
    FA="chr${CHR}.fa"

    curl https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/"$FA".gz --output "$FA".gz

    gunzip "$FA".gz

    faSize $FA -detailed >"$FA".sizes

    bedtools makewindows -g "$FA".sizes -w 3 -s 1 |
        bedtools getfasta -fi "$FA" -bed stdin -tab |
        awk '$2 ~ /^CG/ {split($1, a, ":"); split(a[2], b, "-"); print a[1], b[1], $2}' >"$CHR".mm10.cg.bed

    gzip "$FA"
}

export -f process_chromosome

## Create an empty final output file
touch mm10.ref

## Use parallel to process chromosomes in parallel (one per core)
parallel process_chromosome ::: "${CHROMOSOMES[@]}"

## Combine all the individual chromosome files into the final output
cat *.mm10.cg.bed >> mm10.ref

rm *.cg.bed

## Compress the combined CG sites file
sort -k1,1 -k2,2n mm10.ref | gzip > mm10.ref.gz

rm mm10.ref *fa.gz *fai *sizes


# reference generation end ######################################################################


# data files location / slicing

# this is largely automated by https://github.com/imallona/yamet/blob/8a09723c2e2163aeb643ea16a955d237f201c787/workflow/rules/ecker.smk
# we rely on the `harmonized_ecker_metadata.tsv.gz` output
cd ~/src/yamet/workflow/ # this very commit

snakemake --use-conda --cores 1 ecker_data/harmonized_ecker_metadata.tsv.gz annotation/mm10/done.flag -p

# we skip download and assume CpG reports are at ecker_data/*tsv.tar

# let's inspect the reports

cd $WD
ln -s ~/src/yamet/workflow/ecker_data/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed.tsv.tar cell.tsv.tar
tar xvf cell.tsv.tar 

# (yamet) imallona@mlshainan:~/tmp/yamet_brain$ tar xvf allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed.tsv.tar
# allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed.tsv.gz.idx
# allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H9_AD010_indexed.tsv.gz

zcat allc_*/allc_*_indexed.tsv.gz | grep -w "^19" | grep "CG[ACTG]" | head
# 19      3094546 -       CGG     1       1       1
# 19      3094550 -       CGA     0       1       1
# 19      3102678 -       CGT     1       1       1
# 19      3103765 -       CGT     1       1       1
# 19      3121839 -       CGA     1       1       1
# 19      3123849 -       CGA     1       1       1
# 19      3126120 -       CGG     1       1       1
# 19      3133908 -       CGT     0       1       1
# 19      3140912 +       CGC     2       2       1
# 19      3140959 +       CGC     2       2       1


zcat allc_*/allc_*_indexed.tsv.gz | grep -w "^19" | grep "CG[ACTG]" | grep -E -w "312[0-9]{4}"
# 19      3121839 -       CGA     1       1       1
# 19      3123849 -       CGA     1       1       1
# 19      3126120 -       CGG     1       1       1


# ok, that's strand specific and without `chr` strings, we'll have to harmonize the reference

# let's see how do they encode the cpgs exactly...

zcat allc_*/allc_*_indexed.tsv.gz | grep -w "^19" | grep "CG[ACTG]" | awk '{FS=OFS="\t"; print "chr"$1,$2-2,$2+2,".",".",$3}' | head -10 > theirs.bed
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr19.fa.gz

gunzip chr19.fa.gz
bedtools getfasta -fi chr19.fa -bed theirs.bed -s -tab
# chr19:3094544-3094548(-)        accg
# chr19:3094548-3094552(-)        gtcg
# chr19:3102676-3102680(-)        tacg
# chr19:3103763-3103767(-)        gccg
# chr19:3121837-3121841(-)        accg
# chr19:3123847-3123851(-)        CTCG
# chr19:3126118-3126122(-)        ctcg
# chr19:3133906-3133910(-)        TACG
# chr19:3140910-3140914(+)        ccgc
# chr19:3140957-3140961(+)        gcgc

# so +1,+2 if +; and +2,+3 if -



# what about the annotations?

ln -s ~/src/yamet/workflow/ecker_data/harmonized_ecker_metadata.tsv.gz .

zcat harmonized_ecker_metadata.tsv.gz | cut -f14-16 | head -1000 | sort | uniq -c | sed 's/^ *//g' | sort -V | tail -20

# 9 "Inh" "MGE-Pvalb"     "MGE-Pvalb Gfra2"
# 10 "Exc"        "CT-L6" "CT-L6 Il1rap"
# 11 "Exc"        "CT-L6" "CT-L6 Map4"
# 11 "Inh"        "CGE-Vip"       "CGE-Vip Galnt17"
# 11 "Inh"        "MGE-Sst"       "MGE-Sst Rerg"
# 12 "Exc"        "IT-L6" "IT-L6 Fstl4"
# 12 "Exc"        "IT-L23"        "IT-L23 Foxp1"
# 13 "Exc"        "PT-L5" "PT-L5 Tmtc2"
# 13 "Inh"        "MGE-Pvalb"     "MGE-Pvalb Sema5a"
# 13 "NonN"       "MGC"   "MGC mgc-all"
# 14 "NonN"       "ASC"   "ASC mid"
# 19 "Inh"        "MGE-Sst"       "MGE-Sst Frmd6"
# 21 "Exc"        "IT-L5" "IT-L5 Etv1"
# 23 "Inh"        "MGE-Pvalb"     "MGE-Pvalb Thsd7a"
# 32 NA   NA      NA
# 36 "Inh"        "CGE-Lamp5"     "CGE-Lamp5 Grk5"
# 59 "Exc"        "IT-L5" "IT-L5 Grik3"
# 78 "Exc"        "IT-L5" "IT-L5 Cdh8"
# 125 "Exc"       "IT-L4" "IT-L4 Shc3"
# 284 "Exc"       "IT-L23"        "IT-L23 Cux1"


# looks reasonable
