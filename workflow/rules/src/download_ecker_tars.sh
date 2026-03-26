#!/bin/bash
##
## Retrieves Ecker tar files
##
## Atreya Choudhury
## Started 19thMar 2025

mkdir -p ${snakemake_params[raw]}

while read -r url; do
    filename=$(basename "$url")
    found_files=$(ls "${snakemake_params[raw]}" | grep "^$filename" || true)
    if [[ -n "$found_files" ]]; then
        echo "Skipping $filename as it already exists"
    else
        echo "Downloading $url"
        wget -e robots=off -P "${snakemake_params[raw]}" $url
        echo "Done"
    fi
done <"${snakemake_input[0]}"
