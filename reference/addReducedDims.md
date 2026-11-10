# Add scplainer Component Analysis Results

The function will add the component results computed by
[`scpComponentAnalysis()`](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-ComponentAnalysis.md)
to a `SingleCellExperiment`'s `reducedDims` slot, to all using the many
`scater` functions, such as
[`scater::plotReducedDim()`](https://rdrr.io/pkg/scater/man/plotReducedDim.html),
[`scater::plotTSNE()`](https://rdrr.io/pkg/scater/man/plot_reddim.html),
...

## Usage

``` r
addReducedDims(sce, x)
```

## Arguments

- sce:

  An instance of class
  [SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- x:

  A `List` of `DataFrames` containing principal components. This list is
  typically the `bySample` element produced by
  [`scpComponentAnalysis()`](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-ComponentAnalysis.md).

## Value

A `SingleCellExperiment` with updated `reducedDims`.

## Author

Laurent Gatto and Christophe Vanderaa

## Examples

``` r

library("scater")
#> Loading required package: SingleCellExperiment
#> Loading required package: scuttle
data("leduc_minimal")
pcs <- scpComponentAnalysis(
   leduc_minimal, method = "ASCA",
   effects = "SampleType")$bySample

reducedDims(leduc_minimal)
#> List of length 0
#> names(0): 
leduc_minimal <- addReducedDims(leduc_minimal, pcs)
reducedDims(leduc_minimal)
#> List of length 3
#> names(3): unmodelled residuals ASCA_SampleType
plotReducedDim(leduc_minimal, dimred = "ASCA_SampleType",
               colour_by = "SampleType")

leduc_minimal <- runTSNE(leduc_minimal, dimred = "ASCA_SampleType")
plotTSNE(leduc_minimal, colour_by = "SampleType")
```
