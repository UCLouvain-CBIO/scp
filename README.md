[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/scp/coverage.svg?branch=master)](https://codecov.io/github/UCLouvain-CBIO/scp?branch=master)


# Single cell proteomics data processing

The `scp` package is used to process and analyse mass
spectrometry-based single cell proteomics data.  It relies on the
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/)
package to manage and process
[`SingleCellExperiment`](http://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
objects.

## Installation

To install the stable version from Bioconductor:

```
if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
BiocManager::install("scp")
```

To install the development version from GitHub:

```
BiocManager::install("UCLouvain-CBIO/scp")
```

## Asking for help

Feel free to use [Github
issues](https://github.com/UCLouvain-CBIO/scp/issues) or the
[Bioconductor support site](https://support.bioconductor.org/) to ask
question or report problems with `scp`.
