# barnlarv

[![Travis-CI Build Status](https://travis-ci.org/jinyung/barnlarv.svg?branch=master)](https://travis-ci.org/jinyung/barnlarv)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1068124.svg)](https://doi.org/10.5281/zenodo.1068124)

This is a companion `R` package to the paper:

> Wong, J. Y., Chan, K. Y. K., & Chan, B. K. K. (2018). Phylogenetic, ecological and biomechanical constraints on larval form: A comparative morphological analysis of barnacle nauplii. *PLoS ONE*, 13(11): e0206973. [doi:10.1371/journal.pone.0206973](https://doi.org/10.1371/journal.pone.0206973)

Findings reported in the paper can be reproduced with the datasets and codes from this package. 

---

## Installation:
- using tarball from [releases](https://github.com/jinyung/barnlarv/releases) page*. In `R`:

```R
install.packages("<path to downloaded zip/tarball>",repos = NULL)
```
- using `devtools`. In `R`:

```R
devtools::install_github("jinyung/barnlarv")
```

\*Note:

  1. binaries for Windows not provided, `devtools` method should work across platforms. 

  2. [releases](https://github.com/jinyung/barnlarv/releases) page also contains package manual.

---

## Using data from the package 

### Outlines and sampled outlines dataset
Interspecific barnacle nauplii outline and sampled outline datasets can be found in [`data`](data) directory and called directly in `R`, e.g.:

```R
data(stage2outline)
data(stage2landmark)
```
see `?outline` or `?semilandmark` for help files of these datasets. These datasets were pre-processed from raw images stored in [`inst/extdata/stage2/larvae_drawing`](inst/extdata/stage2/larvae_drawing) and [`inst/extdata/ontogeny/larvae/drawing`](inst/extdata/stage2/larvae_drawing) directories. The processing can be reproduced following the codes stored in [`data-raw`](data-raw) directory. In [`data-raw`](data-raw), the `.tps` files of the sampled outlines are also available and can be imported into other geometric morphometrics software. 

### Covariates and phylogeny data

Data for covariates used in the main analyses can be found [here](inst/extdata/stage2/frontal-horn-summary.csv) or called in `R` with:

```R
read.csv(system.file("extdata/stage2", "frontal-horn-summary.csv",
         package = "barnlarv"))
```

Data of phylogeny for a subset of species can be found [here](inst/extdata/perez-tree-37-species) or called in `R` with:

```R
ape::read.tree(system.file("extdata/perez-tree-37-species", 
               package = "barnlarv"))
```

---

## Reproducing the results

All tables and figures can be reproduced following the dynamic `R markdown` documents, which can be download [here](https://github.com/jinyung/barnlarv/releases/download/v0.0.3/tables_and_figures.zip).

There are three files inside `.zip`:
1. `figures_and_tables.Rmd`
2. `supplementary_figures_and_tables.Rmd`
3. `preamble.tex`

file 1 will produce figures and tables in the main paper; file 2 will produce the supplementary tables and figures; file 3 is used for formatting in `knitr` compilation of `Rmd` files for pdf outputs. Compiling the `Rmd` files will reproduce the results*.

\***NOTE**: the figures, tables, and supplements were created with `barnlarv v0.0.2`, as such can only be 100% reproduced with the specific `R` session environment used when the paper was written. See [v0.0.2 release notes](https://github.com/jinyung/barnlarv/releases/tag/v0.0.2) for details.
