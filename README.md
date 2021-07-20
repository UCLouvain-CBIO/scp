[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/scp/coverage.svg?branch=master)](https://codecov.io/github/UCLouvain-CBIO/scp?branch=master)
[![R-CMD-check-bioc](https://github.com/UCLouvain-CBIO/scp/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/UCLouvain-CBIO/scp/actions?query=workflow%3AR-CMD-check-bioc)

# Single cell proteomics data processing

The `scp` package is used to process and analyse mass
spectrometry-based single cell proteomics data.  It relies on the
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/)
package to manage and process
[`SingleCellExperiment`](http://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
objects.

## Installation

To install the **stable version** from Bioconductor:

```
if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
BiocManager::install("scp")
```

To install the **development version** from GitHub, you first need to 
ensure that you are using the `devel` release of Bioconductor and make
sure your installed libraries are valid. Then, you can install `scp` 
from Github.

```
## Get the latests version of Bioconductor
BiocManager::install(version = "devel")
stopifnot(BiocManager::valid())
## Get the latests version of scp from GitHub
BiocManager::install("UCLouvain-CBIO/scp")
```

## Asking for help

Feel free to use [Github
issues](https://github.com/UCLouvain-CBIO/scp/issues) or the
[Bioconductor support site](https://support.bioconductor.org/) to ask
question or report problems with `scp`.
