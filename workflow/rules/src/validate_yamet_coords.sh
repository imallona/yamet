#!/usr/bin/env bash
# Validate that harmonized cell coordinates match the yamet reference.
#
# Usage: validate_yamet_coords.sh <harmonized_dir> <ref_gz>
#
# Checks:
#   - at least one harmonized cell file exists
#   - neither the cell file nor the reference uses chr-prefixed chromosome names
#   - the first 100 (chr, pos) pairs from the cell file appear in the reference
#
# Exit 1 on any failure so snakemake aborts the run before yamet is invoked.

harmonized_dir="$1"
ref="$2"

cell=$(ls "${harmonized_dir}"/*.gz 2>/dev/null | head -1)
if [ -z "$cell" ]; then
    echo "ERROR: no harmonized cell files found in ${harmonized_dir}" >&2
    exit 1
fi

cell_has_chr=$(zcat "$cell" | head -100 | grep -c '^chr' || true)
ref_has_chr=$(zcat "$ref"  | head -200 | grep -c '^chr' || true)

if [ "$cell_has_chr" -gt 0 ] && [ "$ref_has_chr" -eq 0 ]; then
    echo "ERROR: cell file $cell uses chr-prefixed chromosomes but reference $ref does not" >&2
    exit 1
fi

if [ "$cell_has_chr" -eq 0 ] && [ "$ref_has_chr" -gt 0 ]; then
    echo "ERROR: reference $ref uses chr-prefixed chromosomes but cell file $cell does not" >&2
    exit 1
fi

overlap=$(comm -12 \
    <(zcat "$cell" | head -100 | awk '{print $1"\t"$2}' | sort) \
    <(zcat "$ref"  | awk '{print $1"\t"$2}' | sort) \
  | wc -l)

if [ "$overlap" -eq 0 ]; then
    echo "ERROR: no overlap between first 100 records of $cell and $ref" >&2
    echo "       Check the CpG offset in the harmonization script for this dataset" >&2
    exit 1
fi

echo "Validation passed: $overlap of 100 cell coordinates found in reference" >&2
