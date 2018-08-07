#!/bin/bash
##
## methtuples to look for nucleosome positioning
## 29th june 2018

export HOME=/home/imallona
export TASK="lsh"
export VIRTENVS="$HOME"/virtenvs
export WD="$HOME"/"$TASK"

mkdir -p $WD
cd $_

# get a sampling from stadlers' data

ln -s /home/imallona/mnt/nfs/wgbs_nina/SRR299053/SRR299053_bwameth_default.bam

# samtools index SRR299053_bwameth_default.bam
# samtools view -h SRR299053_bwameth_default.bam chr14 | samtools view -BS - > chr14.bam

# run methtuple

source $VIRTENVS/methtuple/bin/activate

methtuple -h
# oops, does it work with outputs distinct from bismarks?

deactivate
