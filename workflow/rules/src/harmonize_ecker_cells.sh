#!/bin/bash
##
## Extract Ecker tar files and convert allc to yamet format (CG context only)
##
## Atreya Choudhury / Izaskun Mallona
## Started March 2025

set -euo pipefail

mkdir -p "${snakemake_params[harmonized]}"

process_one_cell() {
    tar_file="$1"
    harmonized_dir="$2"
    chr_filter="$3"

    cell_name=$(basename "$tar_file" .tar)
    harmonized_out="${harmonized_dir}/${cell_name}.gz"

    if [[ -f "$harmonized_out" ]]; then
        echo "Skipping $cell_name"
        return
    fi

    tmpdir=$(mktemp -d)
    tar xf "$tar_file" -C "$tmpdir" 2>/dev/null || true
    allc_file=$(find "$tmpdir" -name "*.tsv.gz" | head -1)

    if [[ -z "$allc_file" ]]; then
        echo "Warning: no tsv.gz in $tar_file, skipping"
        rm -rf "$tmpdir"
        return
    fi

    zcat "$allc_file" \
      | awk -v chr_filter="$chr_filter" '
            BEGIN { OFS="\t" }
            $7 == 1 {
                if (chr_filter == "True" && $1 != "10") next
                cpg_start = ($3 == "-") ? $2 - 2 : $2 - 1
                print $1, cpg_start, cpg_start + 1, ".", ".", $3, $5, $6
            }' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -c 7,8 -o sum \
      | awk 'BEGIN { OFS="\t" } {
            meth_bin = ($5 > 0 && $4 / $5 > 0) ? 1 : 0
            print $1, $2, $4, $5, meth_bin
        }' \
      | gzip -c > "$harmonized_out"

    rm -rf "$tmpdir"
}

export -f process_one_cell

# check there are tar files to process
shopt -s nullglob
tar_files=("${snakemake_params[raw]}"/*.tsv.tar)
shopt -u nullglob

if [[ ${#tar_files[@]} -eq 0 ]]; then
    echo "No tar files in ${snakemake_params[raw]}"
    exit 0
fi

printf '%s\n' "${tar_files[@]}" \
    | xargs -P "${snakemake_params[threads]}" -I{} \
        bash -c 'process_one_cell "$@"' _ {} \
            "${snakemake_params[harmonized]}" \
            "${snakemake_params[chr10_only]}"
