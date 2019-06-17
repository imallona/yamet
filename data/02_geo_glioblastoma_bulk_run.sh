#!/bin/bash
##
## https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=100351
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100351
##
## Izaskun Mallona
## 14 june 2019

export HOME=/home/imallona
export TASK="cg_shadows"
export SUBTASK="glioblastoma_bulk"
export WD=/home/Shared_s3it/imallona/"$TASK"/"$SUBTASK"

export SOFT="$HOME"/soft
export VIRTENVS="$HOME"/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export SAMTOOLS=/usr/local/bin/samtools

export NTHREADS=20
export MIN_COV=10

mkdir -p $WD
cd $_

cd $WD

## requested access to https://www.ebi.ac.uk/ega/datasets/EGAD00001004074 on june the 14th
