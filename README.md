[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov.io](https://codecov.io/github/UCLouvain-CBIO/scp/coverage.svg?branch=master)](https://codecov.io/github/UCLouvain-CBIO/scp?branch=master)
[![R-CMD-check-bioc](https://github.com/UCLouvain-CBIO/scp/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/UCLouvain-CBIO/scp/actions?query=workflow%3AR-CMD-check-bioc)
[![Bioconductor-devel Build Status](https://bioconductor.org/shields/build/devel/bioc/scp.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/scp/)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)


# Single cell proteomics data processing

The `scp` package is used to process and analyse mass
spectrometry-based single cell proteomics data.  It relies on the
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/)
package to manage and process
[`SingleCellExperiment`](http://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
objects.

<img
src="https://raw.githubusercontent.com/UCLouvain-CBIO/scp/master/sticker/sticker.png"
height="150">

## Installation

To install the **stable version** from Bioconductor:

```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("scp")
```

To install the **devel version** from GitHub, you first need to
ensure that you are using the `devel` release of Bioconductor and make
sure your installed libraries are valid.

```r
BiocManager::install(version = "devel")
stopifnot(BiocManager::valid())
```

Then, you can install `scp` from Github.

```r
BiocManager::install("UCLouvain-CBIO/scp")
```

## Citation

To cite the `scp` package in publications use:

> Vanderaa, Christophe, and Laurent Gatto. 2023. “Revisiting the
  Thorny Issue of Missing Values in Single-Cell Proteomics.” Journal
  of Proteome Research 22 (9): 2775–84.

> Vanderaa Christophe and Laurent Gatto. The current state of
  single-cell proteomics data analysis. Current Protocols 3 (1):
  e658.; doi: https://doi.org/10.1002/cpz1.658 (2023).

> Vanderaa Christophe and Laurent Gatto. Replication
  of Single-Cell Proteomics Data Reveals Important
  Computational Challenges. Expert Review of
  Proteomics, 1–9 (2021).

## Asking for help

Feel free to use [Github
issues](https://github.com/UCLouvain-CBIO/scp/issues) or the
[Bioconductor support site](https://support.bioconductor.org/) to ask
question or report problems with `scp`.

## License

The `scp` code is provided under a permissive
[Artistic 2.0 license](https://opensource.org/licenses/Artistic-2.0).
The documentation, including the manual pages and the vignettes, are
distributed under a
[CC BY-SA license](https://creativecommons.org/licenses/by-sa/2.0/).

## Contributions

Contributions are welcome, and should ideally be provdied through a
Github pull request. Feel free to discuss any more non-trivial
suggestions or changes first in an issue.

### Contributors

- [Micha Birklbauer](https://github.com/michabirklbauer): [fix grep
  pattern in example for
  readSingleCellExperiment](https://github.com/UCLouvain-CBIO/scp/pull/94).