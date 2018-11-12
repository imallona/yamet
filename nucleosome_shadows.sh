#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 29th june 2018

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

# samtools index SRR299053_bwameth_default.bam
# samtools view -h SRR299053_bwameth_default.bam chr14 | samtools view -BS - > chr14.bam

# run methtuple

wget https://github.com/PeteHaitch/methtuple/blob/master/data/pe_directional_1.fq.gz_bismark_bt2_pe.bam?raw=true \
     --output-document pe_bt2.bam

wget https://github.com/PeteHaitch/methtuple/blob/master/data/se_directional.fq.gz_bismark_bt2.bam?raw=true \
     --output-document se_bt2.bam

samtools sort pe_bt2.bam  > pe_sorted.bam
# mv -f sorted.bam pe_bt2.bam
samtools index pe_sorted.bam

samtools view -H pe_bt2.bam  | fgrep PG

ln -s /home/sorjuela/ScienceCloud/test1_pe.bam
samtools view -s 0.05 -b test1_pe.bam > sorjuela.bam

bam=sorjuela.bam

source $VIRTENVS/methtuple/bin/activate

methtuple -m 2 --methylation-type CG $bam

deactivate

## methtuple outputs are
# chr	strand	pos1	pos2	MM	MU	UM	UU
# 1      2       3       4      5       6      7        8

# fn=pe_bt2.CG.2.tsv
fn=sorjuela.CG.2.tsv

sed '1d' $fn | \
    awk 'BEGIN {OFS=FS="\t"}{
      print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,$5/($5+$6+$7+$8)*1000,$2
     }' | "$BEDTOOLS" sort > "$(basename $fn .tsv)"_scored.bed


sed '1d' $fn | \
    awk 'BEGIN {OFS=FS="\t"}{
      print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,((2*$5)+$6+$7)/(2*($5+$6+$7+$8))*1000,$2
     }' | "$BEDTOOLS" sort > "$(basename $fn .tsv)"_meth.bed

## to call metilene the scores should be like betavalues (from 0 to 1000)
## so divide by 1000

## check the cpg offset, whether 0 or 1! setting up the middle point here

# for bed in  "$(basename $fn .tsv)"_scored.bed  "$(basename $fn .tsv)"_meth.bed
# do
#     awk '{OFS=FS="\t"; print $1, int($2+($3-$2)/2), $5/1000}' $bed > \
#         "$(basename $bed .be)"_metilene.input    
# done

## rather, using starts and ends and NOT merging strands @todo
for bed in  "$(basename $fn .tsv)"_scored.bed  "$(basename $fn .tsv)"_meth.bed
do
    awk '{OFS=FS="\t"; print $1,$2,$3,$5/1000}' $bed > \
        "$(basename $bed .be)"_metilene.input    
done


perl ${METILENE_INPUT} -in1 "$(basename $fn .tsv)"_scored.bed_metilene.input \
     -in2  "$(basename $fn .tsv)"_meth.bed_metilene.input \
     -out "$(basename $fn .tsv)"_metilene \
     -h1 agreement -h2 methylation -b $BEDTOOLS


$METILENE -t $NTHREADS -d 0.1 -m 10 -a agreement -b methylation \
           "$(basename $fn .tsv)"_metilene > no_filtering.file



