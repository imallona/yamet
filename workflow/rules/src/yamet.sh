#!/bin/bash
##
## Run yamet
##
## Atreya Choudhury
## Started 27thFeb 2025

mkdir -p ${snakemake_params[base]}

${snakemake_input[yamet]} \
    --cell ${snakemake_input[cells]} \
    --reference ${snakemake_input[ref]} \
    --intervals ${snakemake_input[intervals]} \
    --cores ${snakemake[threads]} \
    --print-sampens F \
    --out ${snakemake_output[out]%.gz} \
    --det-out ${snakemake_output[det_out]%.gz} \
    --norm-det-out ${snakemake_output[norm_det_out]%.gz} \
    --meth-out ${snakemake_output[meth_out]%.gz}

gzip ${snakemake_output[out]%.gz}
gzip ${snakemake_output[det_out]%.gz}
gzip ${snakemake_output[norm_det_out]%.gz}
gzip ${snakemake_output[meth_out]%.gz}
