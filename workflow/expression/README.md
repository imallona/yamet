# RNA-seq (Bian 2018 scTrio-seq2 CRC, GSE97693)

Standalone Snakefile for per-cell transcriptomes. Independent of the DNA-meth
pipeline in `workflow/`. One cell equals one SRR for the scTrio-seq2 libraries.

## Cell selection

Cells enter the pipeline when their GEO `Sample_title` matches
`SAMPLE_TITLE_PREFIX` (currently `scTrioSeq2Rna_`). Single point of assignment:
rule `parse_sctrioseq2_rna_cells`, which joins the GSE97693 series matrix to
SRA runinfo on GSM and emits `data/rna/meta/sctrioseq2_rna_cells.tsv`
(srr, gsm, title). On GSE97693 this yields 1112 cells.

Downsampling to `N_CELLS` is a seeded random shuffle in rule
`select_rna_subset` (`SUBSET_SEED`, default 42). This avoids the
single-patient / single-metastasis bias of an alphabetical pick. Extending to
the multiplexed protocol later means adding a second title prefix and a
`protocol` column in `parse_sctrioseq2_rna_cells` rather than dispersing the
logic across rules.

## Reference build

GRCh38 primary assembly, GENCODE v47 transcripts, decoy-aware salmon index.
DNA-meth stays on hg19; RNA-to-meth joins happen at the gene level (Ensembl
IDs or symbols), so the build mismatch is a non-issue.

## Analysis plan

Two-stage target:

- Default (`all`): raw FastQC + MultiQC on N_CELLS cells. Inspect before
  deciding what to trim.
- `quantify`: trim, post-trim FastQC, salmon quant per cell, mapping-rate
  gate at `MIN_MAPPING_RATE` (60%), full MultiQC.

Rules defined but not run by default: `trim_adapters_and_primers`,
`fastqc_trimmed`, `download_gencode_references`, `build_salmon_index`,
`salmon_quant`, `check_mapping_rate`, `multiqc_full`.

## Trimming decision

Inspect raw MultiQC and act on these signals:

- Adapter Content non-flat toward 3' ends: keep trim_galore (it auto-detects
  Illumina adapters). If flat, trim_galore is still a cheap safety net.
- Overrepresented Sequences containing `AAGCAGTGGTATCAACGCAGAGT...` (Tang 2009
  pre-amp primer / Smart-seq TSO anchor): enable the cutadapt `-g` step.
  If absent, skip cutadapt.
- Per-base sequence content biased past position ~20: 5' primer or TSO
  carry-over; add a hard clip (`--clip_R1 15 --clip_R2 15` in trim_galore) and
  rerun post-trim FastQC.
- Quick carry-over sniff per FASTQ:
  `zcat R1.fastq.gz | awk 'NR%4==2' | head -1000000 | awk '/^AAGCAGTGGTATCAACGCAGAG/{c++} END{print c}'`
  More than ~1% of the million sampled reads starting with the anchor means
  cutadapt is warranted.

Post-trim verification: adapter content flat, no overrepresented primer/TSO
entries, per-base content flat past position 20.

## Mapping-rate gate

`check_mapping_rate` reads `aux_info/meta_info.json` per cell and fails if any
`percent_mapped` is below `MIN_MAPPING_RATE`. If a cell fails, do not silently
raise the threshold; first check trimming and library-type strandedness
(`libType A` is auto, but a forced `ISR`/`ISF` may be more reliable on
scRNA-seq).

## Disk footprint

A 10-cell default run on GSE97693 pulls ~2.8 GB of .sra and decompresses to
~10-12 GB of FASTQ. salmon index build adds another ~25 GB between the
gentrome download, the uncompressed FASTAs and the k-mer index.

## Running

From `workflow/expression/`:

    source ~/miniconda3/bin/activate
    conda activate snakemake
    snakemake -s Snakefile --use-conda --cores 10
    # review data/rna/qc/multiqc_raw/multiqc_report.html
    snakemake -s Snakefile --use-conda --cores 10 quantify

Conda envs live in `../envs/rnaseq.yml`. When that file changes, drop the
cached env directory before rerunning; the hash-based trigger is not reliable.

## Outputs

- `data/rna/meta/sctrioseq2_rna_cells.tsv` full cell list (srr, gsm, title)
- `data/rna/meta/sctrioseq2_rna_cells_subset.tsv` N_CELLS subset (seeded)
- `data/rna/raw/{srr}_{1,2}.fastq.gz`
- `data/rna/qc/multiqc_raw/multiqc_report.html`
- `data/rna/trimmed/{srr}_{1,2}.trim.fastq.gz`
- `data/rna/qc/multiqc_full/multiqc_report.html`
- `data/rna/quant/{srr}/quant.sf`
- `data/rna/quant/mapping_rate_summary.tsv`
