# yamet: Yet Another Methylation Entropy Tool

`yamet` is a command line tool written in C++ for computing entropies from methylation data.

`yamet` takes cell methylation data, reference site definitions and intervals of interest as input, then computes sample entropy and 2-mer shannon entropy. The sample entropy metric is computed for every cell at every interval of interest and also aggregated per cell. It is a measure of the entropy within individual cells. The 2-mer shannon entropy metric is computed for every interval of interest and is a measure of the entropy across cells. Both metrics take into account missing/uncovered positions in the methylation files using a reference file of positions.

Please also check (probably side branches of):

- https://github.com/emsonder/MethQuant
- https://github.com/emsonder/MethQuant-analysis

<!-- prettier-ignore -->
> [!NOTE]
> `yamet` is currently under development! :confetti_ball:

## Installation

### Brew

**yamet** can be installed via `brew` on MacOS and Ubuntu

```bash
brew tap atchox/brew

# stable version
brew install yamet

# or build latest version from source
brew install --HEAD yamet
```

### Compiled binaries

Compiled binaries can be downloaded from the [releases](https://github.com/imallona/yamet/releases) page.

### Build from source

```bash
git clone https://github.com/imallona/yamet.git
cd yamet/method
bash build.sh
./build/yamet --help
```

## Usage

`yamet` processes (covered CpG) DNA methylation report(s), a reference file listing all CpG positions in a genome, and a bedfile specifying the genomic regions to calculate scores for (e.g. promoters, genes, etc). Full CLI args:

```text
yamet

input:
  -c [ --cytosine_report ] arg          per-cell cytosine report file(s) 
                                        (resembling Bismark's for covered 
                                        cytosines). Synonyms: --cell, -c.
                                        These tab-separated files describe 
                                        covered CpGs sorted by chromosome and 
                                        position, with columns:
                                          chr   pos   meth_reads   total_reads 
                                          rate
  -r [ --cytosine_locations ] arg       genomic locations of all cytosines 
                                        (typically CpGs). Used to account for 
                                        missing values in cytosine 
                                        reports.Synonyms: --cytosine_locations,
                                        --reference, -r.
                                        Required to reconstruct contiguous CpG 
                                        sequences. Columns:
                                          chr   pos
  -i [ --regions ] arg                  BED file defining genomic regions where
                                        entropies will be computed (both per 
                                        feature and globally)Synonyms: 
                                        --regions, --features, --target, 
                                        --intervals, -i.
                                        Columns:
                                          chr   start   end
  --skip-header-all [=arg(=1)]          Number of header lines to skip in all 
                                        input files (default: 1). Synonyms: 
                                        --skip-header-all, --skip-header.
  --skip-header-cytosine_report [=arg(=1)]
                                        Number of header lines to skip in 
                                        cytosine_report/cell files (default: 
                                        1).
  --skip-header-cytosine_locations [=arg(=1)]
                                        Number of header lines to skip in 
                                        cytosine_locations/reference file 
                                        (default: 1).
  --skip-header-regions [=arg(=1)]      Number of header lines to skip in 
                                        regions/features/target/intervals file 
                                        (default: 1).

output:
  -d [ --det-out ] arg                  (optional) path to detailed output file
  -n [ --norm-det-out ] arg             (optional) path to detailed normalized 
                                        output file
  -m [ --meth-out ] arg                 (optional) path to average methylation 
                                        output file
  --all-meth [=arg(=true)] (=false)     if true, include all CpGs in 
                                        methylation summaries, including those 
                                        not used for template construction 
                                        (default: false).
  --all-meth [=arg(=true)] (=false)
  -o [ --out ] arg                      (optional) path to simple output file

resource utilisation:
  --cores arg (=0)                      number of cores used for simultaneously
                                        parsing methylation files
  --chunk-size arg (=64K)               size of the buffer (per file) used for 
                                        reading data. Can be specified as a 
                                        positive integer (bytes) or with a 
                                        suffix: B (bytes), K (kilobytes), M 
                                        (megabytes), G (gigabytes). Example: 
                                        4096, 64K, 128M, 2G

verbose:
  --print-intervals                     print parsed intervals file
  --print-reference                     print parsed reference file
  --print-sampens [=arg(=true)] (=true) print computed sample entropies

misc:
  -h [ --help ]                         print help message
  --version                             current version information

```

## Repository

- `.github/workflows` and `test`: testing and build automations
- `method`: yamet code. See [our releases](https://github.com/imallona/yamet/releases) (including binaries)
- `simulations`: yamet simulations
- `workflow`: yamet applications, except simulations
- `yamet-r`: yamet R package

You might want to browse the issues and PRs to explore current developments.


## Running the analysis / reproducible figures

For simulations:

```bash
cd simulations
snakemake --use-conda --cores NUM_CORES
```

For other analysis (mind the slow data download):

```bash
cd workflow
snakemake --use-conda --cores NUM_CORES
```

## License

GPLv3

## Contact

izaskun.mallona at gmail.com

## Acknowledgements

We thank [Mark D. Robinson at UZH](https://robinsonlabuzh.github.io/) for partly funding this project.

We thank the [Swiss National Science Foundation](https://data.snf.ch/grants/grant/190824) and the [Graduate Campus at the University of Zurich](https://www.grc.uzh.ch/en/funding) for their support.

`yamet` builds on great FOSS components, including [snakemake](https://snakemake.readthedocs.io/en/stable/), and uses public data from epigenomic atlases and scientific publications. Thank you!
