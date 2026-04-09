"""
Figure assembly rules for manuscript figures.

Each rule reads pre-computed RDS files saved to the workflow
"""

import os.path as op

MANUSCRIPT_FIGURES = op.join("..", "..", "sc_dna_methylation_entropy_2025", "figures")

rule render_fig_ecker:
    conda:
        op.join("..", "envs", "r.yml")
    input:
        op.join(ECKER_BASE, "results", "ecker.html")
    output:
        html    = op.join(ECKER_BASE, "results", "fig_ecker.html"),
        panels  = op.join(MANUSCRIPT_FIGURES, "fig_ecker_panels.pdf"),
        heatmap = op.join(MANUSCRIPT_FIGURES, "fig_ecker_heatmap.pdf"),
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
        html    = op.join(ARGELAGUET_BASE, "results", "fig_argelaguet.html"),
        panels  = op.join(MANUSCRIPT_FIGURES, "fig_argelaguet_panels.pdf"),
        heatmap = op.join(MANUSCRIPT_FIGURES, "fig_argelaguet_heatmap.pdf"),
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
        html    = op.join(CRC, "results", "fig_crc.html"),
        panels  = op.join(MANUSCRIPT_FIGURES, "fig_crc_panels.pdf"),
        heatmap = op.join(MANUSCRIPT_FIGURES, "fig_crc_heatmap.pdf"),
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
        html = op.join(CRC, "results", "fig_crc_diffentropy.html"),
        pdf  = op.join(MANUSCRIPT_FIGURES, "fig_crc_diffentropy.pdf"),
    params:
        out_dir       = MANUSCRIPT_FIGURES,
        corrected_sce = op.join(CRC, "results", "sce_windows_colon_corrected_10000.rds")
    threads:
        4
    log:
        log = op.join("logs", "fig_crc_diffentropy.log")
    script:
        "src/fig_crc_diffentropy.Rmd"
