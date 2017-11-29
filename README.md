[![Travis-CI Build Status](https://travis-ci.org/jinyung/barnlarv.svg?branch=master)](https://travis-ci.org/jinyung/barnlarv)


# barnlarv

## About
This is a companion `R` package to the manuscript:

> Wong, J. Y., Chan, K.Y.K., & Chan, B. K. K. (2017). Phylogenetic, ecological and biomechanical constraints on larval form: A comparative morphological analysis of barnacle nauplii.

Findings reported in the manuscript can be reproduced with the datasets and codes from this package.

## Installation:
- using tarball/zip from [releases](releases) page. In `R`:

```R
install.packages("<path to downloaded zip/tarball>",repos = NULL)
```
- using `devtools`. In `R`:

```R
devtools::install_github("jinyung/barnlarv")
```

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

## Reproducing the results

All tables and figures can be reproduced following the dynamic `R markdown` documents:

- [Main figures and tables](inst/doc/figures_and_tables.Rmd) **
- [Supplementary figures and tables](inst/doc/supplementary_figures_and_tables.Rmd) **

These files will produce the figures and tables used for submission (`knitr` engine required for compilation)

---
** Note: will be available after manuscript review (>`v0.0.1` release)
