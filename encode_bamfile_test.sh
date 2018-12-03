#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 29th june 2018

## @todo parallelize by chromosome!

export HOME=/home/imallona
export TASK="cg_shadows"
export WD="$HOME"/mnt/nfs/"$TASK"

export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export METILENE="$SOFT"/metilene/metilene_v0.2-7/metilene
export METILENE_INPUT="$SOFT"/metilene/metilene_v0.2-7/metilene_input.pl
export METILENE_OUTPUT="$SOFT"/metilene/metilene_v0.2-7/metilene_output.pl

export SAMTOOLS=/usr/local/bin/samtools

NTHREADS=20

mkdir -p $WD/data
cd $_


function sort_by_queryname {
    bam="$1"
    nthreads="$2"

    $SAMTOOLS sort -@ "$nthreads" -n $bam -o "$(basename $bam .bam)"_n_sorted.bam
}

## @param the path to the bismark bamfile
## output writes a methtuple output
## methtuple outputs are
# chr	strand	pos1	pos2	MM	MU	UM	UU
# 1      2       3       4      5       6      7        8
function run_methtuple {
    bam="$1"
    
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

function filter_by_coverage {
    fn="$1"
    mincov="$2"

    
awk -v mincov="$mincov" '
NR==1{next}
{
  if (($5+$6+$7+$8) >= mincov) print $0


}' $fn

}

## unnormalized Shannon entropies (could return E/log(4) to be normalized)
function test_entropy_string {
    
 cat << EOF > /tmp/forentropy2
header
1 2 3 4
1 1 1 1
0 0 0 1
1 1 1 1
1 2 3 4
100000000 0 0 1
EOF


awk '
NR==1{next}
{
  H["MM"] = $1
  H["UU"] = $2
  H["UM"] = $3
  H["MU"] = $4
  N  = ($1+$2+$3+$4)
  E = 0
  p = ""
  for (i in H) {
    if (H[i] != 0) {
      p = H[i]/N;
      E -=  p * log(p);
    }
  }
print E


}'  /tmp/forentropy2

rm -f /tmp/forentropy2

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



## parses run_methtupleoutputs to get entropyvalues
## @param the methtuple path
## output writes almost a bed6 (score column is the entropy and DOES NOT go from 0 to 1000)
## chr     strand  pos1    pos2    MM      MU      UM      UU
    
function get_entropy {
    ## entropy
    fn="$1"
awk '
BEGIN {OFS=FS="\t"}
NR==1{next}
{
  H["MM"] = $5
  H["UU"] = $8
  H["UM"] = $7
  H["MU"] = $6
  N  = ($5+$6+$7+$8)
  E = 0
  p = ""
  for (i in H) {
    if (H[i] != 0) {
      p = H[i]/N;
      E -=  p * log(p);
    }
  }
print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,E,$2

}' $fn | "$BEDTOOLS" sort > "$(basename $fn .tsv)"_entropy.bed    
}

function collapse_strands {
    echo 'todo'
}

sample=ENCFF857QML

mkdir -p $sample

cd $_

echo Uncomment to retrieve data
# 78 GB
# wget https://www.encodeproject.org/files/ENCFF857QML/@@download/ENCFF857QML.bam

bam="$sample".bam

sort_by_queryname "$bam" "$NTHREADS"

mv -f "$(basename $bam .bam)"_n_sorted.bam "$bam"

run_methtuple "$bam"

fn="$sample".CG.2.tsv

filter_by_coverage $fn 5 > foo
mv -f foo $fn

get_entropy $fn

get_methylation $fn

## get chromatin colors and look for association

echo 'Retrieving chromatin colors'

mysql --user=genome --host=genome-mysql.cse.ucsc.edu  \
      --skip-column-names \
      -A -e \
      'SELECT chrom, chromStart, chromEnd, name 
     FROM hg19.wgEncodeAwgSegmentationChromhmmH1hesc;' | \
    "$BEDTOOLS" sort  > h1_hmm.bed

# for color in "$(cut -f 4 h1_hmm.bed | sort | uniq)"
# do
#     grep "$color" h1_hmm.bed > curr.bed

#     "$BEDTOOLS" jaccard -a  \
#             -b curr.bed
        
# done

# ## get top quartiles for entropy and dnameth
# ## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003181
## mind the strand collapsing at the groupbystep
$BEDTOOLS intersect -wa -wb \
          -a "$sample".CG.2_entropy.bed \
          -b h1_hmm.bed | \
    $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_entropy_colored.bed


$BEDTOOLS intersect -wa -wb \
          -a "$sample".CG.2_meth.bed \
          -b h1_hmm.bed | \
    $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_meth_colored.bed

