#!/bin/bash
#
# intersects pairs of meth/entropies bedfiles to correlate them
# needs to be snakemaked
# 26 june 2019


intersect () {
    a="$1"
    b="$2"
    intersected="$3"
    bna=$(basename $a)
    bnb=$(basename $b)

    echo -e $bna'_entropy\t'$bna'_meth\t'$bnb'_entropy\t'$bnb'_meth\t' > header
    bedtools intersect -wa -wb -a "$a" \
             -b "$b" | cut -f5,7,12,14 > foo_"$bna"_"$bnb"

    cat header foo_"$bna"_"$bnb" | gzip -c > "$intersected"
    rm header foo_"$bna"_"$bnb"
}

# a=/home/Shared_s3it/imallona/meth_entropies_snakemake/ENCFF079RGH/ENCFF079RGH_cov_20_entropy_and_meth.bed.gz
# b="$a"

# intersect $a $b foo.bed.gz


## finding all files and getting all the pairwise intersects 

WD=/home/Shared_s3it/imallona/meth_entropies_snakemake
cd $WD/sandbox

fns=$(find $WD -name "*entropy_and_meth.bed.gz" | tr -s '\n' ' ')
echo $fns

for fn in $fns
do
    for fn2 in $fns
    do
        echo $fn
        echo $fn2
        intersect  $fn $fn2 intersect_"$(basename $fn)"_"$(basename $fn2)"
    done
done
