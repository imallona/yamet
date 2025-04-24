# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Upcoming

### Added

- [method] Default values for some of the less important arguments in user-exposed functions in the yamet library

### Fixed

- [method] k-mer shannon entropy log base changed from $e$ to $2k$ to squeeze output to $[0,\;1]$

## [v1.1.0-rc.3](https://github.com/imallona/yamet/releases/tag/v1.1.0-rc.3)

### Added

- [method] Added error capturing (incorrect line parsing)
- [method] Added `--skip-header`, `--skip-header-cell`, `--skip-header-reference` and `--skip-header-intervals` flags for skipping headers in files at various levels of granularity
- [method] Added `--meth-out` to produce average methylation files per region and cell
- [workflow] Added nonfunctional/untested workflow steps
- [simulations] Added traditional simulations to evaluate robustness against sparsity and average DNA methylation levels dependencies

## [v1.1.0-rc.2](https://github.com/imallona/yamet/releases/tag/v1.1.0-rc.2)

This is the second official release candidate of **yamet**

### Added

- build tools for generating brew bottles
- precompiled binary for linux-aarch
- separated yamet into a library and an executable that builds on it
- compilation for a shared and static library that one can link to
- cmake config so that the library can be discovered and linked using function like `find_package` and `target_link_libraries`
- `--chunk-size` param to control the size of the buffers into which the files are read
- improved algorithm for handling partial lines in buffers

## [v1.1.0-rc.1](https://github.com/imallona/yamet/releases/tag/v1.1.0-rc.1)

### Overview

This is the first official release candidate of **yamet** 🎉

**yamet** (Yet Another Methylation Entropy Tool) is a powerful and efficient tool written in C++ for computing methylation entropy from genomic data. It aims to provide researchers with fast, accurate, and scalable entropy calculations while being easy to integrate into bioinformatics workflows.

### Added

#### Core Features

- sample entropy within cells - per search interval and also aggregated at a cell level
- average methylation within and across cells - per search interval across all cells and also aggregated at a cell level
- k-mer shannon entropy - per search interval across all cells
- cells files and reference files can be in gzipped format
- multithreaded

#### Primary Inputs

All files must be tab separated. The cell files are the cytosine reports for all covered positions. The reference file is a list of all positions in the genome. The intervals file is list of regions of interest within each chromosome. We do not enforce particular file formats but we do require the data to be presented in the following format.

- cell files with 5 columns in the following order:
  1. chromosome
  2. position
  3. number of methylated reads
  4. total number of reads
  5. rate (beta value)
- reference file with 2 columns in the following order:
  1. chromosome
  2. position
- target intervals with 3 columns in the following order:
  1. chromosome
  2. start position
  3. end position

## [v0.1.0-rc.1](https://github.com/imallona/yamet/releases/tag/v0.1.0-rc.1)

### Capabilities

This is the first release candidate of `yamet` v0.1.0 after refactoring and rewriting in C++.

- Offers a sample entropy within cells, Shannon's entropy across cells, average methylation by/across cells in C++
- Includes simplistic tests on simulated data and valgrind and cppcheck profiling

### Known issues

Requires GCC > 13.

### Installation

```{bash}
cd method
bash build.sh
```

### Usage

- CLI takes a reference file listing cytosine coordinates, as many (covered) cytosine reports as cells, and a bedfile to filter in regions to calculate the metrics from. Metrics are calculated per bedfile interval.
- CLI help:

```
input:
  -t [ --tsv ] arg                      tab separated files for different cells
                                        in the following format

                                         1    5    0    2    0
                                         1    9    1    1    1
                                         2    2    3    4    1

                                        where the columns are the chromosome,
                                        position, number of methylated reads,
                                        total number of reads and the rate
                                        respectively
  -r [ --ref ] arg                      tab separated file for reference sites
                                        in the following format

                                         1    5     7
                                         1    7     9
                                         1    9     11
                                         1    11    13
                                         2    2     4
                                         2    4     6

                                        where the columns are the chromosome,
                                        start position and the end position
                                        respectively
  -b [ --bed ] arg                      path to bed file for regions of
                                        interest in the following format

                                         1    5     7
                                         1    10    30
                                         2    1     6

                                        where the columns are the chromosome,
                                        start position and the end position
                                        respectively

output:
  -d [ --det-out ] arg                  (optional) path to detailed output file
  -o [ --out ] arg                      (optional) path to simple output file

resource utilisation:
  --n-cores arg (=0)                    number of cores used for simultaneously
                                        parsing methylation files
  --n-threads-per-core arg (=1)         number of threads per core used for
                                        simultaneously parsing methylation
                                        files

verbose:
  --print-bed                           print parsed regions file
  --print-ref                           print parsed reference file
  --print-sampens [=arg(=true)] (=true) print computed sample entropies

misc:
  -h [ --help ]                         produce help message
  --version                             current version information

```
