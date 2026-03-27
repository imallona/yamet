#!/bin/bash
##
## Convert Argelaguet et al. scNMT-seq methylation files to yamet format.
##
## Input format (met/cpg_level/*.tsv.gz):  chr  pos  met_reads  nonmet_reads  rate
## Output format (yamet):                  chr  pos  m  t  meth_bin
##
## where m = met_reads, t = met_reads + nonmet_reads,
## meth_bin = 1 if met_reads > 0 else 0.
##
## CpG positions are already strand-collapsed in the source files,
## so no strand correction is needed.

set -euo pipefail

process_cell() {
    id_met="$1"
    raw_dir="$2"
    harmonized_dir="$3"
    chr10_only="$4"

    src="${raw_dir}/${id_met}.tsv.gz"
    dst="${harmonized_dir}/${id_met}.gz"

    [[ -f "$dst" ]] && return 0
    [[ ! -f "$src" ]] && { echo "Warning: $src not found, skipping"; return 0; }

    zcat "$src" | awk -v chr10="$chr10_only" '
      NR > 1 {
        if (chr10 == "True" && $1 != "10") next
        OFS = "\t"
        print $1, $2, $3, $3 + $4, ($3 > 0 ? 1 : 0)
      }
    ' | sort -k1,1 -k2,2n | gzip -c > "$dst"
}
export -f process_cell

mkdir -p "${snakemake_params[harmonized]}"

id_met_col=$(
    zcat "${snakemake_input[meta]}" |
        head -1 |
        tr '\t' '\n' |
        awk '/^id_met$/ { print NR }'
)

zcat "${snakemake_input[meta]}" |
    awk -F'\t' -v col="$id_met_col" 'NR > 1 && $col != "NA" { print $col }' |
    sort -u |
    xargs -P "${snakemake_params[threads]}" -I{} \
        bash -c 'process_cell "$@"' _ {} \
            "${snakemake_params[raw]}" \
            "${snakemake_params[harmonized]}" \
            "${snakemake_params[chr10_only]}"
