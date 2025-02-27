#!/bin/bash
##
## Retrieves only meth cytosine reports from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97693
## Single-cell Multi-omics Sequencing and Analyses of Human Colorectal Cancer
##
## Started 20thDec 2024

mkdir -p ${snakemake_params[raw]}

while read -r gsm; do
    short="$(echo $gsm | cut -c1-7)"
    url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"$short"nnn/"$gsm"/suppl/
    found_files=$(ls "${snakemake_params[raw]}" | grep "^$gsm" || true)
    if [[ -n "$found_files" ]]; then
        echo "Skipping $gsm (already exists: $found_files)"
    else
        echo "Downloading $gsm"
        wget -e robots=off -r -k -A gz -nd -P "${snakemake_params[raw]}" $url
        echo "Done"
    fi
done <"${snakemake_input[gsm]}"
