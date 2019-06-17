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

function get_data {
    wget https://github.com/PeteHaitch/methtuple/blob/master/data/pe_directional_1.fq.gz_bismark_bt2_pe.bam?raw=true \
         --output-document pe_bt2.bam

    wget https://github.com/PeteHaitch/methtuple/blob/master/data/se_directional.fq.gz_bismark_bt2.bam?raw=true \
         --output-document se_bt2.bam

    samtools sort pe_bt2.bam  > pe_sorted.bam
    # mv -f sorted.bam pe_bt2.bam
    samtools index pe_sorted.bam

    samtools view -H pe_bt2.bam  | fgrep PG

    # ln -s /home/sorjuela/ScienceCloud/test1_pe.bam
    # samtools view -s 0.05 -b test1_pe.bam > sorjuela.bam
}


## @param the path to the bismark bamfile
## output writes a methtuple output
## methtuple outputs are
# chr	strand	pos1	pos2	MM	MU	UM	UU
# 1      2       3       4      5       6      7        8
function run_methtuple {
    bam="$1 "
    
    source $VIRTENVS/methtuple/bin/activate

    methtuple -m 2 --methylation-type CG $bam

    deactivate
}

## parses run_methtupleoutputs to get meth agreement scores
## @param the methtuple path
## output writes a bed6 (score column goes to 0 to 1000)

function get_agreement_score {
    fn="$1"
    sed '1d' $fn | \

        awk 'BEGIN {OFS=FS="\t"}{
      print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,$5+$8/($5+$6+$7+$8)*1000,$2
     }' | "$BEDTOOLS" sort > "$(basename $fn .tsv)"_scored.bed
}

## parses run_methtupleoutputs to get meth betavalues
## @param the methtuple path
## output writes a bed6 (score column goes to 0 to 1000)
function get_methylation {
    fn="$1"
    
    sed '1d' $fn | \
        awk 'BEGIN {OFS=FS="\t"}{
      print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,(($5+(0.5*$6)+(0.5*$7))/($5+$6+$7+$8))*1000,$2
     }' | "$BEDTOOLS" sort > "$(basename $fn .tsv)"_meth.bed    
}

function prepare_metilene_input {
    agreement_score_fn="$1"
    background_score_fn="$2"
    sample_name=$3

    ## to call metilene the scores should be like betavalues (from 0 to 1000)
    ## so divide by 1000

    ## check the cpg offset, whether 0 or 1! setting up the middle point here

    # for bed in  "$(basename $fn .tsv)"_scored.bed  "$(basename $fn .tsv)"_meth.bed
    # do
    #     awk '{OFS=FS="\t"; print $1, int($2+($3-$2)/2), $5/1000}' $bed > \
    #         "$(basename $bed .be)"_metilene.input    
    # done

    ## rather, using starts and ends and NOT merging strands @todo
    
    ## rather, using starts and ends and NOT merging strands @todo
    for bed in $agreement_score_fn $background_score_fn
    do
        awk '{OFS=FS="\t"; print $1,$2,$3,$5/1000}' $bed > \
            "$bed"_metilene.input    
    done

    perl ${METILENE_INPUT} -in1 "$agreement_score_fn"_metilene.input \
         -in2  "$background_score_fn"_metilene.input \
         -out "$sample_name"_metilene \
         -h1 agreement -h2 background -b $BEDTOOLS
}


function call_metilene {
    agreement_tag="$1"
    background_tag="$2"
    metilene_input="$3"
    
    $METILENE -t $NTHREADS -d 0.1 -m 10 -a "$agreement_tag" -b "$background_tag" \
              "$metilene_input" > "$metilene_input".output

}

function postprocess_metilene {
    agreement_tag="$1"
    background_tag="$2"
    dmr_input="$3"
    
    perl "${METILENE_OUTPUT}" \
         -q "$dmr_input" \
         -c 5 \
         -d 0.05 \
         -a "$agreement_tag" \
         -b "$background_tag" 
}

sample=sorjuela
bam="$sample".bam

run_methtuple "$bam"

fn="$sample".CG.2.tsv

get_agreement_score $fn

get_methylation $fn

prepare_metilene_input "$sample".CG.2_scored.bed "$sample".CG.2_meth.bed "$sample"

call_metilene agreement background "$sample"_metilene

postprocess_metilene agreement background "$sample"_metilene.output

awk '{OFS=FS="\t"; print "chr"$1,$2,$3,"foo",$4*1000}' metilene_qval.0.05.bedgraph | \
    grep -P "^(chr[0-9]{1,2})\t" | fgrep - > minus.bedgraph

awk '{OFS=FS="\t"; print "chr"$1,$2,$3,"foo",$4*1000}' metilene_qval.0.05.bedgraph | \
    grep -P "^(chr[0-9]{1,2})\t" | fgrep -v - > plus.bedgraph 

echo 'track type=bedGraph name="plus" description="plus" visibility=full color=200,100,0 altColor=0,100,200 priority=20' > foo

cat foo plus.bedgraph > plus2.bedgraph

echo 'track type=bedGraph name="minus" description="minus" visibility=full color=200,100,0 altColor=0,100,200 priority=20' > foo

cat foo minus.bedgraph > minus2.bedgraph


## let's correlate the values with ctcf from encode (assuming it's colon! and not sure the
# bam file is hg19)

# https://www.encodeproject.org/experiments/ENCSR236YGF/ 
cat << EOF >> ctcf_colon.url
https://www.encodeproject.org/metadata/type=Experiment&files.file_type=bed+narrowPeak/metadata.tsv -X GET -H "Accept: text/tsv" -H "Content-Type: application/json" --data '{"elements": ["/experiments/ENCSR236YGF/"]}'
https://www.encodeproject.org/files/ENCFF719EWH/@@download/ENCFF719EWH.bed.gz
https://www.encodeproject.org/files/ENCFF302QQX/@@download/ENCFF302QQX.bed.gz
https://www.encodeproject.org/files/ENCFF559LUB/@@download/ENCFF559LUB.bed.gz
https://www.encodeproject.org/files/ENCFF654IUM/@@download/ENCFF654IUM.bed.gz
EOF

xargs -L 1 curl -O -L < ctcf_colon.url

gunzip *bed.gz

for fn in $(find -name "*ENC*bed")
do
    $BEDTOOLS sort -i $fn > foo
    mv foo $fn
done

$BEDTOOLS sort -i plus.bedgraph > foo
mv foo plus.bedgraph

"$BEDTOOLS" jaccard -a plus.bedgraph \
            -b ENC*.bed

$BEDTOOLS sort -i minus.bedgraph > foo
mv foo minus.bedgraph

"$BEDTOOLS" jaccard -a minus.bedgraph \
            -b ENC*.bed



for fn in $(find -name "*ENC*bed")
do
    "$BEDTOOLS" reldist \
                -a plus.bedgraph \
                -b $fn > "$fn"_plus.reldist
done


for fn in $(find -name "*ENC*bed")
do
    "$BEDTOOLS" reldist \
                -a minus.bedgraph \
                -b $fn > "$fn"_minus.reldist
done


## note that EWH and QQX are peaks and LUB and IUM are inputs


## what if dnase encode data?
cat << EOF >> dnasei_colon.url
https://www.encodeproject.org/metadata/type=Experiment&files.file_type=bed+broadPeak&files.file_type=bed+narrowPeak/metadata.tsv -X GET -H "Accept: text/tsv" -H "Content-Type: application/json" --data '{"elements": ["/experiments/ENCSR462XTM/"]}'
https://www.encodeproject.org/files/ENCFF071XTK/@@download/ENCFF071XTK.bed.gz
https://www.encodeproject.org/files/ENCFF421NEH/@@download/ENCFF421NEH.bed.gz
https://www.encodeproject.org/files/ENCFF330CTT/@@download/ENCFF330CTT.bed.gz
https://www.encodeproject.org/files/ENCFF006EIZ/@@download/ENCFF006EIZ.bed.gz
https://www.encodeproject.org/files/ENCFF722LRQ/@@download/ENCFF722LRQ.bed.gz
https://www.encodeproject.org/files/ENCFF257QKM/@@download/ENCFF257QKM.bed.gz
https://www.encodeproject.org/files/ENCFF756XIR/@@download/ENCFF756XIR.bed.gz
https://www.encodeproject.org/files/ENCFF113TVH/@@download/ENCFF113TVH.bed.gz
EOF



xargs -L 1 curl -O -L < dnasei_colon.url

gunzip *bed.gz

for fn in $(find -name "*ENC*bed")
do
    $BEDTOOLS sort -i $fn > foo
    mv foo $fn
done

$BEDTOOLS sort -i plus.bedgraph > foo
mv foo plus.bedgraph

"$BEDTOOLS" jaccard -a plus.bedgraph \
            -b ENC*.bed

$BEDTOOLS sort -i minus.bedgraph > foo
mv foo minus.bedgraph

"$BEDTOOLS" jaccard -a minus.bedgraph \
            -b ENC*.bed



for fn in $(find -name "*ENC*bed")
do
    "$BEDTOOLS" reldist \
                -a plus.bedgraph \
                -b $fn > "$fn"_plus.reldist
done


for fn in $(find -name "*ENC*bed")
do
    "$BEDTOOLS" reldist \
                -a minus.bedgraph \
                -b $fn > "$fn"_minus.reldist
done
