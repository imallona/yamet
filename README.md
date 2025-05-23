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

### Compiled Binaries

Compiled binaries can be downloaded from the [releases](https://github.com/imallona/yamet/releases) page.

### Build from source

```bash
git clone https://github.com/imallona/yamet.git
cd yamet/method
bash build.sh
./build/yamet --help
```

## Usage

`yamet` processes (covered CpG) DNA methylation report(s) (`-cell` argument), a reference file listing all CpG positions in a genome (`--reference`), and a bedfile specifying the genomic regions to calculate scores for (`--intervals`; e.g. promoters, genes, etc). Full CLI args:

```text
Usage:
  yamet (-c <cell>... | <cell>...) \
        -r <reference> \
        -i <intervals> \
        [OPTIONS]

Required inputs:
  -c --cell <cell>...                 One or more tab-separated methylation files,
                                      OR provide them directly as positional arguments.
  <cell>...                           (Positional alternative) Cell files (same format as above).
  -r --reference <reference>          Reference file, tab-separated and sorted by chromosome/position.
  -i --intervals <intervals>          BED file of intervals of interest.

Optional output:
  -d --det-out <file>                 Path to detailed output file.
  -m --meth-out <file>                Path to average methylation output file.
  -o --out <file>                     Path to simple output file.

Resource options:
  --cores <n>                         Number of cores for parallel parsing
                                      [default: 0, implying program decides].
  --chunk-size <size>                 Buffer size per file (e.g., 64K, 128M, 2G) [default: 64K].

Optional input control:
  --skip-header[=<n>]                 Skip <n> lines in all input files [default: 1].
  --skip-header-cell[=<n>]            Skip <n> lines in cell files (overrides --skip-header).
  --skip-header-reference[=<n>]       Skip <n> lines in reference file (overrides --skip-header).
  --skip-header-intervals[=<n>]       Skip <n> lines in intervals file (overrides --skip-header).

Verbose/debugging:
  --print-intervals                   Print parsed intervals file.
  --print-reference                   Print parsed reference file.
  --print-sampens[=<true|false>]      Print computed sample entropies [default: true].

Miscellaneous:
  -h --help                           Show help message.
  --version                           Show version information.
```

## Repository

- `method`: yamet code. See [our releases](https://github.com/imallona/yamet/releases) (including binaries)
- `workflow`: yamet applications, including simulations
- `.github/workflows` and `test`: testing
- `old`: archived codebase from old version to be removed

## License

GPLv3
