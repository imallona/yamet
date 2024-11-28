#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 17 feb

## @todo parallelize by chromosome!
## Data from https://www.encodeproject.org/experiments/ENCSR617FKV/

export HOME=/home/imallona
export TASK="cg_shadows"
export SUBTASK="encode_bulk"
export WD=/home/Shared_s3it/imallona/"$TASK"/"$SUBTASK"

export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export SAMTOOLS=/usr/local/bin/samtools

export NTHREADS=20
export MIN_COV=10

mkdir -p $WD
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
## ideally should split into two lines, first and scond meth value (do they overlap?)
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
    
function get_entropy_4_states {
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





## parses run_methtupleoutputs to get entropyvalues
## @param the methtuple path
## output writes almost a bed6 (score column is the entropy and DOES NOT go from 0 to 1000)
## chr     strand  pos1    pos2    MM      MU      UM      UU
    
function get_entropy_3_states {
    ## entropy
    fn="$1"
awk '
BEGIN {OFS=FS="\t"}
NR==1{next}
{
  H["MM"] = $5
  H["UU"] = $8
  H["UM"] = $7+$8
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



cd $WD

echo "sample identifiers/manifest for downlading" ##############################################

## id,exp_url,bam_url,assembly,description,replicate,end
cat << EOF > 01_encode.conf
ENCFF112TXF,https://www.encodeproject.org/experiments/ENCSR888FON/,https://www.encodeproject.org/files/ENCFF112TXF/@@download/ENCFF112TXF.bam,GRCh38,IMR90,1,single
ENCFF957OIM,https://www.encodeproject.org/experiments/ENCSR881XOU/,https://www.encodeproject.org/files/ENCFF957OIM/@@download/ENCFF957OIM.bam,GRCh38,HepG2,1,paired
ENCFF572KNK,https://www.encodeproject.org/experiments/ENCSR881XOU/,https://www.encodeproject.org/files/ENCFF572KNK/@@download/ENCFF572KNK.bam,GRCh38,HepG2,2,paired
ENCFF193RVP,https://www.encodeproject.org/experiments/ENCSR550RTN/,https://www.encodeproject.org/files/ENCFF193RVP/@@download/ENCFF193RVP.bam,GRCh38,HeLa-S3,1,paired
ENCFF845VFH,https://www.encodeproject.org/experiments/ENCSR550RTN/,https://www.encodeproject.org/files/ENCFF845VFH/@@download/ENCFF845VFH.bam,GRCh38,HeLa-S3,2,paired
ENCFF079RGH,https://www.encodeproject.org/experiments/ENCSR440MLE/,https://www.encodeproject.org/files/ENCFF079RGH/@@download/ENCFF079RGH.bam,GRCh38,GM23248,1,paired
ENCFF119ELB,https://www.encodeproject.org/experiments/ENCSR440MLE/,https://www.encodeproject.org/files/ENCFF119ELB/@@download/ENCFF119ELB.bam,GRCh38,GM23248,2,paired
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -f1 -d ',')
    bam_url=$(echo $line | cut -f3 -d ',')

    
    echo "$sample" start "$date"
    mkdir $sample; cd $_

    wget --quiet $bam_url

    bam="$sample".bam

    ## sorting by coordinate and filtering chr19 only
    samtools sort -@ $NTHREADS "$bam" > sorted_"$sample".bam
    mv -f sorted_"$sample".bam $bam
    
    samtools index -@ $NTHREADS  "$bam"
    samtools view -b "$bam" chr19 > mini_"$sample".bam
    mv -f mini_"$sample".bam "$bam"
    
    sort_by_queryname "$bam" "$NTHREADS"

    mv -f "$(basename $bam .bam)"_n_sorted.bam "$bam"

    run_methtuple "$bam"

    fn="$sample".CG.2.tsv

    filter_by_coverage $fn $MIN_COV > foo
    mv -f foo $fn

    get_entropy_4_states $fn

    get_methylation $fn

    
    echo "$sample" end "$date"
    cd ..
done < 01_encode.conf
