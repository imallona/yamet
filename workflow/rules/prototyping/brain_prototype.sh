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

# let's run yamet first. On CpN contexts. Disregarding references altogether.

DATA_PATH=/home/imallona/src/yamet/workflow/ecker_data/
mkdir -p sample

# https://lhqing.github.io/ALLCools/start/input_files.html#columns-in-allc-file
# so 5 is methylated and 6 is coverage
# yamet needs where the columns are the 1 chromosome, 
#                                         2 position, 3 number of methylated reads, 
#                                         4 total number of reads and the rate 
#                                         respectively
# caution ceiling any C with a methylated read! and subsetting chr1!

for cell in allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H6_AD008_indexed.tsv.tar \
                allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H4_AD010_indexed.tsv.tar \
                allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H2_AD012_indexed.tsv.tar \
                allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H5_AD010_indexed.tsv.tar
do
    bn=$(basename $cell .tsv.tar)
    echo $bn
    tar xvf ${DATA_PATH}/${cell}
    
    zcat "$bn"/*tsv.gz | grep -w "^1" |\
        awk '{ OFS=FS="\t";
              if($3=="+")
                print "chr"$1, $2-1, $2,    "", "",$5,$6
              else
                print "chr"$1, $2-2, $2-1,  "", "",$5,$6
              }' | \
        bedtools merge -c 6,7 -o sum  | \
        awk 'BEGIN {OFS=FS="\t"} 
           function ceil(v) { 
              return (v==int(v)) ? v : int(v)+1 
           }
           {
              print $1,$2,$4,$5,ceil($4/$5)
           }' | \
               sort -k1,1 -k2,2n -k3,3n | gzip -c > sample/"$bn".chr1.yametized.gz
done

# now faking a reference borrowing loci from a cell

zcat sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H6_AD008_indexed.chr1.yametized.gz | cut -f1,2 > fake.ref

# chr1 dirty segmentation

TAB="$(printf '\t')"
cat > fake_regions.bed << EOF
chr1${TAB}1${TAB}100000
chr1${TAB}100001${TAB}200000
chr1${TAB}400001${TAB}500000
chr1${TAB}3005608${TAB}4005608
chr1${TAB}4005608${TAB}5005608
EOF

yamet --cell sample/*chr1.yametized.gz \
      --reference fake.ref \
      --intervals fake_regions.bed \
      --out test.yamet.out \
      --cores 50 \
      --det-out test.det.yamet.out \
      --print-sampens F 

# Is able to calculate something, even if with an inadequate reference
# (yamet) imallona@mlshainan:~/tmp/yamet_brain$ tail test.det.yamet.out 
# chr     start   end     sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H2_AD012_indexed.chr1.yametized.gz     sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H4_AD010_indexed.chr1.yametized.gz   sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H5_AD010_indexed.chr1.yametized.gz     sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H6_AD008_indexed.chr1.yametized.gz   shannon avg_meth
# chr1    1       100000  -1      -1      -1      -1      -1      -1
# chr1    100001  200000  -1      -1      -1      -1      -1      -1
# chr1    400001  500000  -1      -1      -1      -1      -1      -1
# chr1    3005608 4005608 0.110306        0.104038        0.0924244       0.0692187       0.317843        0.0303077
# chr1    4005608 5005608 0.0796909       0.089736        0.101747        0.0909062       0.388139        0.0382658
# (yamet) imallona@mlshainan:~/tmp/yamet_brain$ tail test.yamet.out 
# file    sampen  avg_meth
# sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H2_AD012_indexed.chr1.yametized.gz     0.0890592       0.0483029
# sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H4_AD010_indexed.chr1.yametized.gz     0.100165        0.0569034
# sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H5_AD010_indexed.chr1.yametized.gz     0.098031        0.0555058
# sample/allc_180627_CEMBA_mm_P56_P63_5D_CEMBA180612_5D_3_CEMBA180612_5D_4_H6_AD008_indexed.chr1.yametized.gz     0.0799147       0.0322619

:<<EOF
so plan:
evaluate CpH in genebodies and look for between-celltypes differences, and run a gene ontology on those
  genes with differences, using genes-with-calculable-entropy as background

if not enough genes covered: let's focus on DG granule cells, since (from Ecker's 2021):

mCH accumulates throughout the genome during postnatal brain development6,8. We reasoned that DG granule cells, which are continuously replenished by ongoing neurogenesis throughout the lifespan, may accumulate mCH during their post-mitotic maturation. If so, global mCH should correlate with the age and maturity of granule cells. To investigate this, we divided DG granule cells into four groups on the basis of their global mCH levels and investigated regions of differential methylation between the groups. We identified 219,498 gradient CG-DMRs between the four groups, among which 139,387 showed a positive correlation with global mCH (+DMR), and 80,111 were negatively correlated (−DMR) (Fig. 5e). Notably, genes overlapping +DMRs or −DMRs have different annotated functions: genes enriched in +DMRs (+DMRgenes, n = 328) were associated with developmental processes, whereas those enriched in −DMRs (−DMRgenes, n = 112) were related to synaptic function (Extended Data Fig. 10a, b).

So we could check methylation entropy differences for these cells, in CpH contexts in particular
EOF


# ok, let's prototype this further: let's get some excitatory neurons:

mkdir -p exc

for exc in "CLA" \
               "CT-L6" \
               "IT-L23" \
               "IT-L4"
do
    echo $exc
    cell_id=$(basename $(zcat harmonized_ecker_metadata.tsv.gz | grep "$exc" | head -1  |cut -f50 | tr -d '"') .tsv.gz)
    cell_tar=$DATA_PATH/"${cell_id}.tsv.tar"
    tar xvf ${cell_tar}
    mv ${cell_id}/* exc/
done

ll -h exc/

# drwxr-xr-x  2 imallona robinsonlab 4.0K Mar  6 14:37 ./
# drwxr-xr-x 11 imallona robinsonlab 4.0K Mar  6 14:37 ../
# -rw-r--r--  1 imallona robinsonlab  85M Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_1_CEMBA180409_2C_2_A10_AD001_indexed.tsv.gz
# -rw-r--r--  1 imallona robinsonlab  269 Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_1_CEMBA180409_2C_2_A10_AD001_indexed.tsv.gz.idx
# -rw-r--r--  1 imallona robinsonlab  95M Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_1_CEMBA180409_2C_2_A10_AD002_indexed.tsv.gz
# -rw-r--r--  1 imallona robinsonlab  263 Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_1_CEMBA180409_2C_2_A10_AD002_indexed.tsv.gz.idx
# -rw-r--r--  1 imallona robinsonlab 147M Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_3_CEMBA180409_2C_4_G1_AD008_indexed.tsv.gz
# -rw-r--r--  1 imallona robinsonlab  265 Oct 20  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180409_2C_3_CEMBA180409_2C_4_G1_AD008_indexed.tsv.gz.idx
# -rw-r--r--  1 imallona robinsonlab 172M Oct 19  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180410_2C_1_CEMBA180410_2C_2_B1_AD006_indexed.tsv.gz
# -rw-r--r--  1 imallona robinsonlab  275 Oct 19  2018 allc_180508_CEMBA_mm_P56_P63_2C_CEMBA180410_2C_1_CEMBA180410_2C_2_B1_AD006_indexed.tsv.gz.idx

# looks feasible. What about regions for ecker's brain?

# just  non-overlapping 10-kb bins on mm10


mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -s -e \
      'SELECT chrom,size FROM mm10.chromInfo' > mm10.sizes

TAB="$(printf '\t')"

grep -E "chr[0-9XY]{1,2}$TAB" mm10.sizes  > mm10.sizes.can
mv mm10.sizes.can mm10.sizes


bedtools makewindows -g mm10.sizes -w 10000 | sort -k1,1 -k2,2n -k3,3n > mm10_windows.bed
