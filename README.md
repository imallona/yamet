# We have moved!

`yamet` was superseded by [`amet`](https://github.com/imallona/amet) on May 18th 2026.

This project is not longer maintained.


# Aim
`yamet` was a C++ command line tool for computing entropies from DNA methylation data.

It took per-cell cytosine reports, a reference file of all CpG positions, and a BED file of genomic regions, then computed sample entropy (per cell, per interval) and 2-mer Shannon entropy (across cells, per interval). Missing positions were handled via the reference file.

See also:
- https://github.com/emsonder/MethQuant
- https://github.com/emsonder/MethQuant-analysis

## Installation

### Brew

```bash
brew tap atchox/brew
brew install yamet

# or build from source
brew install --HEAD yamet
```

### Compiled binaries

Compiled binaries are available on the [releases](https://github.com/imallona/yamet/releases) page.

### Build from source

```bash
git clone https://github.com/imallona/yamet.git
cd yamet/method
bash build.sh
./build/yamet --help
```

## Usage

```
Usage: yamet -c <cytosine report>... -r <reference> -i <interval> [OPTIONS]

Input:
  -c,--cytosine-report,--cell           per-cell cytosine report file(s), tab-separated,
                                        sorted by chromosome and position:
                                          chr  pos  meth_reads  total_reads  rate
  --metadata                            metadata file (alternative to -c), tab-separated:
                                          cell_id  cluster  path/to/cytosine_report
  -r,--cytosine-locations,--reference   reference file of all CpG positions (required):
                                          chr  pos
  -i,--regions,--intervals              BED file of regions to compute entropies for (required):
                                          chr  start  end
  --skip-header                         lines to skip in all input files
  --skip-header-cytosine-report         lines to skip in cell files
  --skip-header-cytosine-locations      lines to skip in reference file
  --skip-header-regions                 lines to skip in intervals file
  --skip-header-metadata                lines to skip in metadata file

Output:
  -d,--det-out                          path to detailed output file
  -n,--norm-det-out                     path to normalized detailed output file
  -m,--meth-out                         path to average methylation output file
  -o,--out                              path to simple output file
  --all-meth / --templated-meth         include all CpGs in methylation summaries (default: false)

Resources:
  --cores                               cores for parallel file parsing (default: 0, auto)
  --chunk-size                          read buffer size per file, e.g. 64K, 128M (default: 64K)

Verbose:
  --print-intervals                     print parsed intervals file
  --print-reference                     print parsed reference file
  --print-sampens / --no-print-sampens  print computed sample entropies (default: true)

Misc:
  -h,--help                             show help message
  --version                             show version information
```

## Repository layout

- `.github/workflows` and `test`: CI and build automation
- `method`: C++ source. See [releases](https://github.com/imallona/yamet/releases) for binaries.
- `workflow`: Snakemake workflow for real datasets (CRC, Ecker, Argelaguet) and simulations.

## Reproducing analyses

```bash
cd workflow
snakemake --use-conda --cores NUM_CORES
```

Specific targets:

```bash
# coverage simulations
snakemake --use-conda --cores NUM_CORES simulations/results/simulation_figure2.html

# within/between cell simulations
snakemake --use-conda --cores NUM_CORES simulations/results/simulation_08_combined_figure_adj.html
```

Simulations share the Argelaguet gastrulation download, so the tarball is fetched only once.

## License

GPLv3

## Contact

izaskun.mallona at gmail.com

## Acknowledgements

We thank Mark D. Robinson at UZH for partly funding this project, the Swiss National Science Foundation, and the Graduate Campus at the University of Zurich for their support.

`yamet` builds on great open source components including [snakemake](https://snakemake.readthedocs.io/en/stable/) and uses public data from epigenomic atlases.
