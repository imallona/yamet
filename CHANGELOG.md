# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Overview

This is the first official release of **yamet** 🎉

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
