# RNA-seq (Bian 2018 scTrio-seq2 CRC, GSE97693)

Standalone Snakefile for per-cell transcriptomes. Independent of the DNA-meth pipeline in `workflow/`. One cell equals one SRR for the scTrio-seq2 libraries.

## Cell selection

Cells enter the pipeline when their GEO `Sample_title` matches `SAMPLE_TITLE_PREFIX` (currently `scTrioSeq2Rna_`). Single point of assignment: rule `parse_sctrioseq2_rna_cells`, which joins the GSE97693 series matrix to SRA runinfo on GSM and emits `data/rna/meta/sctrioseq2_rna_cells.tsv` (srr, gsm, title, protocol). On GSE97693 this yields 1112 cells.

Downsampling to `N_CELLS` is a seeded random shuffle in rule `select_rna_subset` (`SUBSET_SEED`, default 42). This avoids the single-patient / single-metastasis bias of an alphabetical pick.

## Two protocols under one GEO prefix

Raw FastQC on the downsampled 10 cells revealed that `scTrioSeq2Rna_*` does not map to a single library chemistry. Two batches coexist:

- GSM3241xxx / SRR7461xxx: Clontech SMARTer / Smart-seq2 chemistry. Evidence: overrepresented sequences in R1 dominated by tandem copies of the SMARTer II A TSO (`AAGCAGTGGTATCAACGCAGAGTACATGGG`, ~32% of reads in SRR7461573_1). R2 carries the complementary Clontech Universal / CDS primer (`CCCATGTACTCTGCGTTGATACCACTGCTT`). Per-base GC and content are oriented: R1 sits near the 5' TSO, R2 near the 3' polyA tail.
- GSM2697xxx / SRR5824xxx: Tang 2009 oligo-dT + restriction-site polylinker cassette. Evidence: zero Clontech TSO hits; instead, poly-T and the BamHI-AscI-SalI linker `GGATCCGGCGCGCCGTCGAC` dominate the overrepresented table, plus standard Illumina TruSeq read-through adapter signal.

Both reads are 150 nt cDNA on both protocols, not cell barcodes. The single-cell identity is encoded at the GSM level (one GSM equals one cell), so no demultiplexing is needed.

The split is encoded numerically via `PROTOCOL_GSM_THRESHOLD` (3000000): GSM >= threshold is labelled `smartseq2`, else `tang`. The `protocol` column is written once in `parse_sctrioseq2_rna_cells`, propagated through `select_rna_subset`, and consumed by `trim_adapters_and_primers` via the `protocol_for_srr()` helper. To extend to a third chemistry, add a new key to `TRIM_CUTADAPT_ARGS` and widen the threshold logic; no other rule changes.

## Reference build

GRCh38 primary assembly, GENCODE v47 transcripts, decoy-aware salmon index. DNA-meth stays on hg19; RNA-to-meth joins happen at the gene level (Ensembl IDs or symbols), so the build mismatch is a non-issue.

## Analysis plan

Two-stage target:

- Default (`all`): raw FastQC + MultiQC on N_CELLS cells. Inspect before deciding what to trim.
- `quantify`: trim, post-trim FastQC, salmon quant per cell, mapping-rate gate at `MIN_MAPPING_RATE` (60%), full MultiQC.

Rules defined but not run by default: `trim_adapters_and_primers`, `fastqc_trimmed`, `download_gencode_references`, `build_salmon_index`, `salmon_quant`, `check_mapping_rate`, `multiqc_full`.

## Trimming

Rule `trim_adapters_and_primers` runs cutadapt alone (no trim_galore) with `threads: workflow.cores`. Per-cell protocol is looked up via `protocol_for_srr()` against the subset TSV, and the corresponding `TRIM_CUTADAPT_ARGS` entry is merged with a common backbone:

- Illumina TruSeq read-through on both mates (`-a AGATCGGAAGAGC -A AGATCGGAAGAGC`)
- Quality trim `-q 20`
- `--minimum-length 20`
- smartseq2 protocol: 5' TSO on R1 and its complement on R2.
- tang protocol: 5' polylinker on both mates, 3' poly-A trimming.

Post-trim FastQC is run on every surviving mate and aggregated in `multiqc_full`. Expected verification: adapter content flat, no overrepresented Clontech or polylinker entries, per-base content flat past position 20.

## Mapping-rate gate

`check_mapping_rate` reads `aux_info/meta_info.json` per cell and fails if any `percent_mapped` is below `MIN_MAPPING_RATE`. If a cell fails, do not silently raise the threshold; first check trimming and library-type strandedness (`libType A` is auto, but a forced `ISR`/`ISF` may be more reliable on scRNA-seq).

## Disk footprint

A 10-cell default run on GSE97693 pulls ~2.8 GB of .sra and decompresses to ~10-12 GB of FASTQ. salmon index build adds another ~25 GB between the gentrome download, the uncompressed FASTAs and the k-mer index.

## Running

From `workflow/expression/`:

    source ~/miniconda3/bin/activate
    conda activate snakemake
    snakemake -s Snakefile --use-conda --cores 10
    # review data/rna/qc/multiqc_raw/multiqc_report.html
    snakemake -s Snakefile --use-conda --cores 10 quantify

Conda envs live in `../envs/rnaseq.yml`.

## Outputs

- `data/rna/meta/sctrioseq2_rna_cells.tsv` full cell list (srr, gsm, title, protocol)
- `data/rna/meta/sctrioseq2_rna_cells_subset.tsv` N_CELLS subset (seeded)
- `data/rna/raw/{srr}_{1,2}.fastq.gz`
- `data/rna/qc/multiqc_raw/multiqc_report.html`
- `data/rna/trimmed/{srr}_{1,2}.trim.fastq.gz`
- `data/rna/qc/multiqc_full/multiqc_report.html`
- `data/rna/quant/{srr}/quant.sf`
- `data/rna/quant/mapping_rate_summary.tsv`
