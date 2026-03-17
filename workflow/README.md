## How to run

From within the `Snakefile`-containing folder:

```bash
snakemake --cores 1 --use-conda -p
```

## Analyses

### CRC (human, hg19)

Single-cell bisulfite data from colorectal cancer biopsies (Bian et al. 2018, GSE97693). Cells are grouped by patient and biopsy site (normal colon, primary tumour, lymph node, metastasis). Methylation entropy is computed with yamet over several genomic compartments: genes, LINEs, SINEs, CpG islands, PMDs/HMDs, chromatin states (Segway), ChIP marks (H3K4me3, H3K9me3, H3K27me3), and LADs. A windows-based analysis (10 kb tiles) is also run. Outputs: `data/crc/results/crc.html` (features) and `data/crc/results/crc_windows_*.html` (windows).

Key files: `rules/crc.smk`, `rules/hg19.smk`, `rules/src/crc.Rmd`.

### Ecker (mouse, mm10)

Single-cell methylation atlas of the mouse primary motor cortex (Liu et al. 2021, Ecker lab). Data are downloaded from NeMO. Cells are grouped by major cell type (Exc, Inh, ASC, ODC, OPC, MGC). Only CpG context methylation is used. Entropy is computed over genes, promoters, LINEs, SINEs, and ENCODE ChIP marks (H3K4me3, H3K9me3, H3K27me3, H3K4me1, H3K27ac). Output: `ecker/results/ecker.html`.

Set `ECKER_CHR10_ONLY = True` in `rules/ecker.smk` to restrict processing to chr10 (fast) or `False` for the full genome.

Key files: `rules/ecker.smk`, `rules/mm10.smk`, `rules/src/ecker.Rmd`.
