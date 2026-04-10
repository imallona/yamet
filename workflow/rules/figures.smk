"""
Figure assembly rules for manuscript figures.

Each rule reads pre-computed RDS files saved to the workflow root by
the upstream analysis Rmd scripts and assembles publication-quality
PDF figures. PDFs are saved as side effects inside the Rmd scripts.
Snakemake tracks only the HTML output per rule (Rmd script constraint).
"""

import os.path as op

MANUSCRIPT_FIGURES = op.join("..", "..", "sc_dna_methylation_entropy_2025", "figures")


rule render_fig_ecker:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join(ECKER_BASE, "results", "ecker.html")
    output:
        op.join(ECKER_BASE, "results", "fig_ecker.html")
    params:
        out_dir = MANUSCRIPT_FIGURES
    threads:
        4
    log:
        log = op.join("logs", "fig_ecker.log")
    script:
        "src/fig_ecker.Rmd"


rule render_fig_argelaguet:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join(ARGELAGUET_BASE, "results", "argelaguet.html")
    output:
        op.join(ARGELAGUET_BASE, "results", "fig_argelaguet.html")
    params:
        out_dir = MANUSCRIPT_FIGURES
    threads:
        4
    log:
        log = op.join("logs", "fig_argelaguet.log")
    script:
        "src/fig_argelaguet.Rmd"


rule render_fig_crc:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join(CRC, "results", "crc.html"),
        op.join(CRC, "results", "crc_embeddings_10000.html"),
    output:
        op.join(CRC, "results", "fig_crc.html")
    params:
        out_dir = MANUSCRIPT_FIGURES
    threads:
        4
    log:
        log = op.join("logs", "fig_crc.log")
    script:
        "src/fig_crc.Rmd"


rule render_fig_crc_diffentropy:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        de = op.join(CRC, "results", "de_list_10000.rds"),
        embeddings = op.join(CRC, "results", "crc_embeddings_10000.html"),
    output:
        op.join(CRC, "results", "fig_crc_diffentropy.html")
    params:
        out_dir = MANUSCRIPT_FIGURES,
        corrected_sce = op.join(CRC, "results", "sce_windows_colon_corrected_10000.rds")
    threads:
        4
    log:
        log = op.join("logs", "fig_crc_diffentropy.log")
    script:
        "src/fig_crc_diffentropy.Rmd"
