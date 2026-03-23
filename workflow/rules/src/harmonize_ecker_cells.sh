#!/bin/bash
##
## Extract Ecker tar files and convert allc to yamet format (CG context only)
##
## Atreya Choudhury / Izaskun Mallona
## Started March 2025

set -euo pipefail

mkdir -p "${snakemake_params[harmonized]}"

for tar_file in "${snakemake_params[raw]}"/*.tsv.tar; do
    [[ -f "$tar_file" ]] || { echo "No tar files in ${snakemake_params[raw]}"; break; }

    cell_name=$(basename "$tar_file" .tar)
    harmonized_out="${snakemake_params[harmonized]}/${cell_name}.gz"

    if [[ -f "$harmonized_out" ]]; then
        echo "Skipping $cell_name"
        continue
    fi

    tmpdir=$(mktemp -d)
    tar xf "$tar_file" -C "$tmpdir" 2>/dev/null || true
    allc_file=$(find "$tmpdir" -name "*.tsv.gz" | head -1)

    if [[ -z "$allc_file" ]]; then
        echo "Warning: no tsv.gz in $tar_file, skipping"
        rm -rf "$tmpdir"
        continue
    fi

    zcat "$allc_file" \
      | awk -v chr_filter="${snakemake_params[chr10_only]}" '
            BEGIN { OFS="\t" }
            $4 ~ /^CG/ {
                if (chr_filter == "True" && $1 != "chr10") next
                cpg_start = ($3 == "-") ? $2 - 2 : $2 - 1
                print $1, cpg_start, cpg_start + 1, ".", ".", $3, $5, $6
            }' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -c 7,8 -o sum \
      | awk 'BEGIN { OFS="\t" } {
            meth_bin = ($5 > 0 && $4 / $5 > 0.1) ? 1 : 0
            print $1, $2, $4, $5, meth_bin
        }' \
      | gzip -c > "$harmonized_out"

    rm -rf "$tmpdir"
done
