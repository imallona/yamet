# Current release

Under development! :confetti_ball:

Please also check (probably side branches of):

- https://github.com/emsonder/MethQuant
- https://github.com/emsonder/MethQuant-analysis

# Old release

## Repo organization

- `schemas`, yaml sample and configs, includes `samples.schema.yaml` and `config.schema.yaml`
- `runs`,  data getters, config files and accessors for specific runs, e.g. `encode`, `glios`, `simulations`, `test`
- `snakefiles`, snakefiles for different runs

## Current runs
- `glios`, explores entropies from https://www.ebi.ac.uk/ega/datasets/EGAD00001004074, for which clinical data are available
- `encode`, explores the association between entropies and HMM segmentations using data from ENCODE
- `test`, uses the methtuple bismark bam testfiles


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

Global

- Get a name
