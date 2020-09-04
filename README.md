[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/scp/coverage.svg?branch=master)](https://codecov.io/github/UCLouvain-CBIO/scp?branch=master)


# Single cell proteomics data processing

The `scp` package is used to process and analyse mass
spectrometry-based single cell proteomics data.  It relies on the
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/)
package to manage and process
[`SingleCellExperiment`](http://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
objects.

## Installation

```
if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
BiocManager::install("UCLouvain-CBIO/scp")
```

## Asking for help

Feel free to use [Github
issues](https://github.com/UCLouvain-CBIO/scp/issues) to ask question
and report problems with `scp`.
