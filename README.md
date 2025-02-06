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

# How to run the workflow

Lorem ipsum

```
cd workflow
snakemake --use-conda --cores 1
```

# License

GPLv3
