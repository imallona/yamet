*yamet* is just *y*et *a*nother DNA *m*ethylation *e*ntropy *t*ool. It is aimed to compute DNA methylation heterogeneity scores from bisufilte-sequencing data. It is designed to run standardized snakemake workflows from bismark-mapped BAM files.

This repository is also available at https://bitbucket.org/imallona/yamet/src/master/.

The manuscript is available at https://github.com/imallona/yamet_paper and https://bitbucket.org/imallona/yamet_paper/src/master/ .

## FIXME

the op.join runs, tuple_length, src, entropy_stdin.py does compute normalized entropy and not entropy. This change should be applied everywhere! -- 8th Jan 2020

## Repo organization

- `schemas`, yaml sample and configs, includes `samples.schema.yaml` and `config.schema.yaml`
- `runs`,  snakefiles, data getters, config files and accessors for specific runs, e.g. `encode`, `glios`, `simulations`, `test`

## Current runs

- `glios`, explores entropies from https://www.ebi.ac.uk/ega/datasets/EGAD00001004074, for which clinical data are available
- `encode`, explores the association between entropies and HMM segmentations using data from ENCODE
- `test`, uses the methtuple bismark bam testfiles
- `simulations`
- `tuple_lengths`, on exploring different tuple lengths on encode and/or glios datasets

## To do

Implementation

- Tests
- Config yaml for multisample handling

Metric

- Scores: explore alternatives to methylation-standardized entropy
- Strand specificity, skipped currently

Simulations

- Simulations based en epiallele diversity
- Simulations based in DNA meth
- Simulations based in sequencing depth
- Simulations based in read length

Use cases

- Test dataset provided by `methtuple`
- HMM cell-line specificity, e.g. enhancer detection
- RRBS-derived outcome prediction (glioblastomas)

Use cases challenges

- Confounding with purity/deconvolution of purity
- ASM
- CNV

Further capabilities?

- Differential entropy?
- Regional differential entropy?

Benchmark

- Feinberg's matlab
- Guo et al 2017
- ?
