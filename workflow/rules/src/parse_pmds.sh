#!/bin/bash
##
## Parses PMDs
##
## Started 27thFeb 2025

zcat ${snakemake_input[0]} |
    grep "${snakemake_params[filter]}$" >${snakemake_output[0]}
