#!/bin/bash
##
## process encrypted unaligned bam files to get bismark RRBS-based BAM files
## suitabe for the entropy calculator
##
## GPL
##
## 23rd July 2019
## Izaskun Mallona

## not run yet, untested!!

export SUBTASK="glios_mapping"
export TASK="ega_secure"
export WD="$HOME"/mnt/nfs/"$TASK"/"$SUBTASK"


export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export NTHREADS=16
export BISMARK_NTHREADS=4
## guess this one
export MM10=/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/STARIndex/Ensembl_GRCm38.90.dna.primary_assembly_126/
export MM10_GTF=/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf

export STAR=/home/imallona/soft/star/STAR-2.6.0c/source/STAR
export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$SOFT"/sickle/sickle-1.33/sickle
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export FEATURECOUNTS="$SOFT"/subread/subread-1.6.2-source/bin/featureCounts
export BAMTOOLS="/usr/bin/bamtools"
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export BISMARK="$SOFT"/bismark/Bismark_v0.19.1/bismark
export BOWTIE_PATH="$SOFT"/bowtie/bowtie2-2.3.4.1/

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

export EGACLIENT="$SOFT"/egaclient/EGA_download_client_2.2.2/EgaDemoClient.jar

export HG38_GENOME=/home/Shared_taupo/data/annotation/Human/Ensembl_GRCh38.90/genome/

echo 'decrypt'

PASS=""
java -jar "$EGACLIENT" -p izaskun.mallona@dmmd.uzh.ch $PASS -dc *cip \
     -dck first_attempt


echo 'check all files decrypted' ##########################


for cip in $(find . -name "*bam.cip")
do
    bam="$(basename $cip .cip)"
    if test -f "$bam"
    then
        echo ""
    else
        echo "$bam doesn't"
    fi
done

mkdir -p "$WD"/unmapped_bams
mv "$WD"/../*bam "$WD"/unmapped_bams

echo 'get fastq'

mkdir -p $WD

cd $WD

mkdir -p raw_fastq

for bam in $(find $WD/unmapped_bams -name "*bam" | sort)
do
    echo $bam
    tag=$(basename $bam _unmapped.bam)
    
    echo "$tag" start "$(date)"
    
    mkdir -p "$WD"/raw_fastq/"$tag"

    bamToFastq -i "$bam" \
               -fq "$WD"/raw_fastq/"$tag"/"$tag".fastq

    gzip "$WD"/raw_fastq/"$tag"/"$tag".fastq

    echo "$tag" end "$(date)"
    
    cd $WD
done


echo 'cutadapt and sickle'

for bam in $(find $WD/unmapped_bams -name "*bam" | sort)
do
    echo $bam
    tag=$(basename $bam _unmapped.bam)
    
    echo "$tag" start "$(date)"
    
    r="$WD"/raw_fastq/"$tag"/"$tag"
        
    mkdir -p "$WD"/trimmed_fastq/"$tag"

    trimmed_r="$WD"/trimmed_fastq/"$tag"/"$tag"

    source $VIRTENVS/cutadapt/bin/activate
    
    cutadapt \
	-j $NTHREADS \
	-b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
	-o "$trimmed_r"_cutadapt.fastq.gz \
	"$r".fastq.gz &> "$trimmed_r"_cutadapt.log

    deactivate

    rm -f "$r".fastq.gz 
    
    "$SICKLE" se \
	      -f "$trimmed_r"_cutadapt.fastq.gz \
	      -o "$trimmed_r"_sickle.fastq.gz \
	      -t sanger \
	      -g &> "$trimmed_r"_cutadapt_sickle.log

    rm -f "$r"_cutadapt.fastq.gz
    
    echo "$tag" end "$(date)"
    
    cd $WD
done

## till here

echo 'bismark map'

mkdir -p "$WD"/bismark
cd "$WD"

for fastq in $(find $WD -name "*sickle.fastq.gz" | sort)
do
    echo $fastq

    tag=$(basename $fastq _sickle.fastq.gz)

    bismark_bam="$tag"_bismark_bt2.bam

    mkdir -p "$WD"/bismark/"$tag"
    ( "$BISMARK" --bowtie2 \
                 --path_to_bowtie "$BOWTIE_PATH" \
                 --gzip \
                 --output_dir "$WD"/bismark/"$tag" \
                 --parallel "$BISMARK_NTHREADS" \
                 --genome "$HG38_GENOME" \
                 --se  $fastq ) 2>&1 | tee "$WD"/bismark/"$tag"/"$tag"_bismark.log

done
                    
# echo 'dedup'
                    
# "$(dirname $BISMARK)"/deduplicate_bismark -p \
#                      --output_dir $WD \
#                      --bam "$sample"1_cutadapt_sickle_bismark_bt2_se.bam







