#!/bin/bash
##
## colors by H1 HMM entropy tuples
##
## needs encode_bulk_run.sh to be run first
## 04 march

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

# that might make no sense, using h1 annotations for these elements

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -f1 -d ',')
    bam_url=$(echo $line | cut -f3 -d ',')
    
    echo "$sample" start "$date"

    cd "$WD"/"$sample"
    ## get chromatin colors and look for association

    # ## get top quartiles for entropy and dnameth
    # ## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003181
    ## mind the strand collapsing at the groupbystep
    $BEDTOOLS intersect -wa -wb \
              -a "$sample".CG.2_entropy.bed \
              -b ../hg38_h1_hmm.bed | \
        $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_entropy_colored.bed


    $BEDTOOLS intersect -wa -wb \
              -a "$sample".CG.2_meth.bed \
              -b ../hg38_h1_hmm.bed | \
        $BEDTOOLS groupby -g 1-6 -c 10 -o concat > "$sample".CG.2_meth_colored.bed

    
    echo "$sample" end "$date"
    cd $WD
done < 01_encode.conf
