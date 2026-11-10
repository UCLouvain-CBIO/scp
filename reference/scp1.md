# Single Cell QFeatures data

A small
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object with SCoPE2 data. The object is composed of 5 assays, including 3
PSM-level assays, 1 peptide assay and 1 protein assay.

## Usage

``` r
data("scp1")
```

## Format

An object of class `QFeatures` of length 5.

## Details

The dataset is a subset of the SCoPE2 dataset (version 2, Specht et al.
2019, [BioRXiv](https://www.biorxiv.org/content/10.1101/665307v3)). This
dataset was converted to a
[`QFeatures`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html)
object where each assay in stored as a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object. One assay per chromatographic batch ("LCA9", "LCA10", "LCB3")
was randomly sampled. For each assay, 100 proteins were randomly
sampled. PSMs were then aggregated to peptides and joined in a single
assay. Then peptides were aggregated to proteins.

## Examples

``` r
data("scp1")
scp1
#> An instance of class QFeatures (type: bulk) with 5 sets:
#> 
#>  [1] 190321S_LCA10_X_FP97AG: SummarizedExperiment with 166 rows and 11 columns 
#>  [2] 190222S_LCA9_X_FP94BM: SummarizedExperiment with 176 rows and 11 columns 
#>  [3] 190914S_LCB3_X_16plex_Set_21: SummarizedExperiment with 215 rows and 16 columns 
#>  [4] peptides: SummarizedExperiment with 539 rows and 38 columns 
#>  [5] proteins: SummarizedExperiment with 292 rows and 38 columns 
```
