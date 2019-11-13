# CCmed: cross-condition mediation analysis

The goal of `CCmed` is to provide computationally efficient tools to conduct cross-condition
mediation analysis. `CCmed` identifies trans-associations mediated by cis-effects (e.g. trans-associations of eQTLs mediated by effects on cis-gene expression levels). 
Mediation analyses can be performed at the gene-level or
used to identify trans-associations of complex trait GWAS SNPs/variants.

## Setup: dependencies not on CRAN

Note that in order to use `CCmed`, the `Primo` package must also be installed. Please follow the following steps
to install `Primo`, if not already installed:

`Primo` package uses functions from the `limma` package, which is downloadable from [Bioconductor](https://www.bioconductor.org), and the `lcmix` package, which is downloadable from [R-Forge](https://r-forge.r-project.org). If you have not yet installed the `limma` or `lcmix` packages, please run the following commands prior to installing `Primo`:

  ```R
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  
  install.packages("MASS","matrixStats","nnls","R.methodsS3")
  install.packages("lcmix",repos="http://r-forge.r-project.org")
  ```

Once you have installed `limma` and `lcmix`, you can install and load functions from `Primo`:

  ```R
  devtools::install_github("kjgleason/Primo")
  library("Primo")
  ```

## Setup

Once you have installed `Primo`, you can install and load functions from `CCmed`:

  ```R
  devtools::install_github("kjgleason/CCmed")
  library("CCmed")
  ```

## Citation

To cite `CCmed` in publications, please use:

Fan Yang, Kevin J. Gleason, Jiebiao Wang, The GTEx consortium, Jubao Duan, Xin He, Brandon L. Pierce and Lin S. Chen. CCmed: cross-condition mediation analysis for identifying robust trans-eQTLs and assessing their effects on human traits. bioRxiv (2019), doi:10.1101/803106.
