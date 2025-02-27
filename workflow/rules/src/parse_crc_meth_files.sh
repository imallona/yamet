#!/bin/bash
##
## Parses meth cytosine reports from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
## Single-cell Multi-omics Sequencing and Analyses of Human Colorectal Cancer
## for use with yamet
##
## Started 20thDec 2024

mkdir -p ${snakemake_params[crc]}

ls -1 ${snakemake_params[raw]}/G*${snakemake_params[patient]}_*.txt.gz | while IFS= read -r file; do
    filename=$(basename "$file")
    found_files=$(ls "${snakemake_params[crc]}" | grep "^$filename" || true)
    if [[ -n "$found_files" ]]; then
        echo "Already processed $file"
    else
        echo "Processing $file"
        zcat $file |
            grep "CpG$" |
            awk '
                $9 ~ /^CG/ {
                    if ($4 == "-") $2 = $2 - 2;
                    else if ($4 == "+") $2 = $2 - 1;
                    bin = ($8 > 0.1) ? 1 : 0;
                    print $1, $2, $6, $5, bin;
                }
            ' |
            sort -k1,1 -k2,2n |
            gzip >${snakemake_params[crc]}/$filename
    fi
done
