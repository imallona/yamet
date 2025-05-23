# Description

`yamet` is `y`et `a`nother `m`ethylation `e`ntropy `t`ool.

`yamet` is under development! :confetti_ball:

Please also check (probably side branches of):

- https://github.com/emsonder/MethQuant
- https://github.com/emsonder/MethQuant-analysis

# Repository

- `method`: yamet code. See [our releases](https://github.com/imallona/yamet/releases) (including binaries)
- `workflow`: yamet applications, including simulations
- `.github/workflows` and `test`: testing
- `old`: archived codebase from old version to be removed

# Installation

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

# Usage (CLI arguments)

`yamet` processes (covered CpG) DNA methylation report(s) (`-cell` argument), a reference file listing all CpG positions in a genome (`--reference`), and a bedfile specifying the genomic regions to calculate scores for (`--intervals`; e.g. promoters, genes, etc). Full CLI args:

```bash
$ yamet --help
 
yamet

input:
  -c [ --cell ] arg                     tab separated files, sorted by 
                                        chromosome and position, for different 
                                        cells in the following format
                                        
                                         chr1    5    0    2    0
                                         chr1    9    1    1    1
                                         chr2    2    3    4    1
                                        
                                        where the columns are the chromosome, 
                                        position, number of methylated reads, 
                                        total number of reads and the rate 
                                        respectively
  -r [ --reference ] arg                tab separated file, sorted by 
                                        chromosome and position, for reference 
                                        sites in the following format
                                        
                                         chr1    5     7
                                         chr1    7     9
                                         chr1    9     11
                                         chr1    11    13
                                         chr2    2     4
                                         chr2    4     6
                                        
                                        where the columns are the chromosome, 
                                        start position and the end position 
                                        respectively
  -i [ --intervals ] arg                bed file, sorted by chromosome and 
                                        start position, for intervals of 
                                        interest in the following format
                                        
                                         chr1    5     7
                                         chr1    10    30
                                         chr2    1     6
                                        
                                        where the columns are the chromosome, 
                                        start position and the end position 
                                        respectively
  --skip-header [=arg(=1)]              integer value indicating number of 
                                        lines to skip in all file inputs(this 
                                        is a default value which can be 
                                        overriden by other 'skip-header-*' 
                                        options)
  --skip-header-cell [=arg(=1)]         integer value indicating number of 
                                        lines to skip in the cell 
                                        files(overrides 'skip-header' if 
                                        provided)
  --skip-header-reference [=arg(=1)]    integer value indicating number of 
                                        lines to skip in the reference 
                                        file(overrides 'skip-header' if 
                                        provided)
  --skip-header-intervals [=arg(=1)]    integer value indicating number of 
                                        lines to skip in the intervals 
                                        file(overrides 'skip-header' if 
                                        provided)

output:
  -d [ --det-out ] arg                  (optional) path to detailed output file
  -o [ --out ] arg                      (optional) path to simple output file

resource utilisation:
  --cores arg (=0)                      number of cores used for simultaneously
                                        parsing methylation files
  --threads-per-core arg (=1)           number of threads per core used for 
                                        simultaneously parsing methylation 
                                        files
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
  -h [ --help ]                         produce help message
  --version                             current version information
```

# License

GPLv3
