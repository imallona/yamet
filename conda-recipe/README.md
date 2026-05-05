# yamet bioconda recipe

CI workflow `.github/workflows/conda.yml` builds this recipe on `v*` tags and uploads `.conda` artifacts to the GitHub release.

To publish to bioconda, copy `meta.yaml` and `build.sh` into `recipes/yamet/` of a fork of `bioconda/bioconda-recipes`, replace the `sha256` placeholder with the value printed by CI, and open a PR.

Local build:

```
conda install -n base -c conda-forge conda-build
conda build -c conda-forge -c bioconda conda-recipe
```
