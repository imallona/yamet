#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 17 feb

echo 'todo filter chr19 only!'

## @todo parallelize by chromosome!
## Data from https://www.encodeproject.org/experiments/ENCSR617FKV/

export HOME=/home/imallona
export TASK="cg_shadows"
export WD="$HOME"/mnt/nfs/"$TASK"

export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

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

cat << EOF >> 01_encode.conf
id,exp_url,bam_url,assembly,description,replicate,end
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
    echo "$sample" start

    echo "$sample" end
    
done < 01_encode.conf

:<<EOF

    
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

get_entropy_4_states $fn

get_methylation $fn

## get chromatin colors and look for association

## to be liftovered!
mysql --user=genome --host=genome-mysql.cse.ucsc.edu  \
      --skip-column-names \
      -A -e \
      'SELECT chrom, chromStart, chromEnd, name 
     FROM hg19.wgEncodeAwgSegmentationChromhmmH1hesc;' | \
    "$BEDTOOLS" sort  > hg19_h1_hmm.bed


wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver

chmod a+x liftOver

./liftOver hg19_h1_hmm.bed hg19ToHg38.over.chain hg38_h1_hmm.bed liftover_h1_unmapped.bed

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
          -b hg38_h1_hmm.bed | \
    $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_entropy_colored.bed


$BEDTOOLS intersect -wa -wb \
          -a "$sample".CG.2_meth.bed \
          -b hg38_h1_hmm.bed | \
    $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_meth_colored.bed


## on ctcf and dnasei associations

wget https://www.encodeproject.org/files/ENCFF368LWM/@@download/ENCFF368LWM.bed.gz -O ENCFF368LWM_ctcf.bed.gz

wget https://www.encodeproject.org/files/ENCFF655STT/@@download/ENCFF655STT.bed.gz -O ENCFF655STT_dnasei_broad.bed.gz

wget https://www.encodeproject.org/files/ENCFF030XPN/@@download/ENCFF030XPN.bed.gz -O ENCFF030XPN_dnasei_narrow.bed.gz

## uncompress and intersect
for encode in ENCFF030XPN ENCFF655STT ENCFF368LWM
do
    fn=$(find . -name "$encode*bed.gz")

    gunzip -f $fn

    bed="$(basename $fn .gz)"

    $BEDTOOLS sort -i "$bed" | cut -f1-3 |  grep -P "^(chr[0-9]{1,2})\t" > sorted
    mv -f sorted "$bed"
    
    echo $fn

    ## only autosomal chromosomes
    
    # awk some column valued 0 or not
    awk '{OFS=FS="\t"; if ($5 == 0) print }' "$sample".CG.2_entropy.bed |  \
        grep -P "^(chr[0-9]{1,2})\t" | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > "$sample"_zero_entropy_vs_"$(basename $fn .bed.gz)".reldist
    
    # awk some column valued 0 or not
    awk '{OFS=FS="\t"; if ($5 > 0.5) print }' "$sample".CG.2_entropy.bed | \
        grep -P "^(chr[0-9]{1,2})\t" | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > "$sample"_high_entropy_vs_"$(basename $fn .bed.gz)".reldist

done


## not really interesting, what about rampage? (
## this is hESC H7 and not H1

wget https://www.encodeproject.org/files/ENCFF788SVK/@@download/ENCFF788SVK.bed.gz -O\
     ENCFF788SVK_rampage.bed.gz


for encode in ENCFF788SVK 
do
    fn=$(find . -name "$encode*bed.gz")

    gunzip -f $fn

    bed="$(basename $fn .gz)"

    $BEDTOOLS sort -i "$bed" | cut -f1-3 |  grep -P "^(chr[0-9]{1,2})\t" > sorted
    mv -f sorted "$bed"
    
    echo $fn

    ## only autosomal chromosomes
    
    # awk some column valued 0 or not
    awk '{OFS=FS="\t"; if ($5 == 0) print }' "$sample".CG.2_entropy.bed |  \
        grep -P "^(chr[0-9]{1,2})\t" | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > "$sample"_zero_entropy_vs_"$(basename $fn .bed.gz)".reldist
    
    # awk some column valued 0 or not
    awk '{OFS=FS="\t"; if ($5 > 0.5) print }' "$sample".CG.2_entropy.bed | \
        grep -P "^(chr[0-9]{1,2})\t" | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > "$sample"_high_entropy_vs_"$(basename $fn .bed.gz)".reldist

done

## what if stratifying by dnameth as well?

entropy_fn="$sample".CG.2_entropy.bed
meth_fn="$sample".CG.2_meth.bed

# unmeth and meth split
# paste $entropy_fn $meth_fn | \
#     column -s $'\t' | \
#     cut -f1,2,3,4,6,11,5 | \




paste $entropy_fn $meth_fn | \
    column -s $'\t' | \
    awk '{FS=OFS="\t"; if ($11 < 200 && $5 < 0.1) print $1,$2,$3,$4,$6,$11,$5}' > \
        "$sample"_lowly_meth_low_entropy.bed

paste $entropy_fn $meth_fn | \
    column -s $'\t' | \
    awk '{FS=OFS="\t"; if ($11 < 200 && $5 >= 0.3) print $1,$2,$3,$4,$6,$11,$5}' > \
        "$sample"_lowly_meth_high_entropy.bed


paste $entropy_fn $meth_fn | \
    column -s $'\t' | \
    awk '{FS=OFS="\t"; if ($11 >= 800 && $5 < 0.1) print $1,$2,$3,$4,$6,$11,$5}' > \
        "$sample"_highly_meth_low_entropy.bed

paste $entropy_fn $meth_fn | \
    column -s $'\t' | \
    awk '{FS=OFS="\t"; if ($11 >= 800 && $5 >= 0.3) print $1,$2,$3,$4,$6,$11,$5}' > \
        "$sample"_highly_meth_high_entropy.bed




wget https://www.encodeproject.org/files/ENCFF562OAN/@@download/ENCFF562OAN.bed.gz
gunzip ENCFF562OAN.bed.gz
mv -f ENCFF562OAN.bed ENCFF562OAN_kdm1a.bed

mkdir -p reldist

for item in "$sample"_lowly_meth_low_entropy.bed \
                     "$sample"_highly_meth_low_entropy.bed \
                     "$sample"_lowly_meth_high_entropy.bed \
                     "$sample"_highly_meth_high_entropy.bed
do
    for bed in ENCFF030XPN_dnasei_narrow.bed \
                   ENCFF368LWM_ctcf.bed \
                   ENCFF655STT_dnasei_broad.bed \
                   ENCFF788SVK_rampage.bed \
                   ENCFF562OAN_kdm1a.bed
    do
        
    grep -P "^(chr[0-9]{1,2})\t" $item | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > reldist/"$(basename $item .bed)"__vs__"$(basename $bed .bed.gz)".reldist
    done
done


## the reldists are biased by gc content, i.e. cpgislands signals

# KDM1A


## remove cpgislands just in case


# hg38    Primary Table: cpgIslandExt
mysql --user=genome --host=genome-mysql.cse.ucsc.edu  \
      --skip-column-names \
      -A -e \
      'SELECT chrom, chromStart, chromEnd, name 
     FROM hg38.cpgIslandExt' | \
    "$BEDTOOLS" sort  > hg38_cpgislandext.bed

mkdir -p reldist_no_cpgi

for item in "$sample"_lowly_meth_low_entropy.bed \
                     "$sample"_highly_meth_low_entropy.bed \
                     "$sample"_lowly_meth_high_entropy.bed \
                     "$sample"_highly_meth_high_entropy.bed
do
    for bed in ENCFF030XPN_dnasei_narrow.bed \
                   ENCFF368LWM_ctcf.bed \
                   ENCFF655STT_dnasei_broad.bed \
                   ENCFF788SVK_rampage.bed \
                   ENCFF562OAN_kdm1a.bed
    do
        
        grep -P "^(chr[0-9]{1,2})\t" $item | \
            $BEDTOOLS intersect -v -a stdin -b hg38_cpgislandext.bed | \
        $BEDTOOLS reldist -a stdin \
                  -b "$bed" > reldist_no_cpgi/"$(basename $item .bed)"__vs__"$(basename $bed .bed.gz)".reldist
    done
done


EOF
