# Read single-cell proteomics tabular data

Function to import and convert tabular data from a spreadsheet or a
`data.frame` into a `SingleCellExperiment` and `QFeatures` object.

## Usage

``` r
readSCP(..., experimentsAsSce = FALSE)

readSCPfromDIANN(..., experimentsAsSce = FALSE)

readSingleCellExperiment(...)
```

## Arguments

- ...:

  Parameters passed to
  [`readSummarizedExperiment()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html),
  [`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
  or
  [`readQFeaturesFromDIANN()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeaturesFromDIANN.html).
  See these respective manual pages for details.

- experimentsAsSce:

  A `logical(1)` indicating whether the `QFeatures` object should be
  composed of `SingleCellExperiment` objects. By default: `FALSE`, the
  `QFeatures` object will be composed of `SummarizedExperiment` objects.
  Note that using `SingleCellExperiment` can impact the performance.

## Value

An instance of class `SingleCellExperiment` or a `QFeatures`, composed
of `SingleCellExperiment` or `SummarizedExperiment` objects.

## Note

The `SingleCellExperiment` class is built on top of the
`RangedSummarizedExperiment` class. This means that some column names
are forbidden in the `rowData`. Avoid using the following names:
`seqnames`, `ranges`, `strand`, `start`, `end`, `width`, `element`

## See also

- The more general
  [`QFeatures::readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
  function, which this function depends on.

- The more general
  [`QFeatures::readQFeaturesFromDIANN()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeaturesFromDIANN.html)
  function, for details and an example on how to read label-free and
  plexDIA (mTRAQ) data processed with DIA-NN.

- The
  [`QFeatures::readSummarizedExperiment()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.html)
  function, which `readSingleCellExperiment()` depends on.

- The
  [`SingleCellExperiment::SingleCellExperiment()`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  class.

## Examples

``` r

######################################################
## Load a single acquisition as a SingleCellExperiment

## Load a data.frame with PSM-level data
data("mqScpData")

## Create the QFeatures object
quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData))
sce <- readSingleCellExperiment(mqScpData, quantCols)
sce
#> class: SingleCellExperiment 
#> dim: 1361 16 
#> metadata(0):
#> assays(1): ''
#> rownames(1361): 1 2 ... 1360 1361
#> rowData names(133): uid Sequence ... residual participated
#> colnames(16): Reporter.intensity.1 Reporter.intensity.2 ...
#>   Reporter.intensity.15 Reporter.intensity.16
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

######################################################
## Load multiple acquisitions as a QFeatures

## Load an example table containing MaxQuant output
data("mqScpData")

## Load the (user-generated) annotation table
data("sampleAnnotation")

## Format the tables into a QFeatures object
readSCP(assayData = mqScpData,
        colData = sampleAnnotation,
        runCol = "Raw.file")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Splitting data in runs.
#>   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.
#> An instance of class QFeatures (type: scp) with 4 sets:
#> 
#>  [1] 190222S_LCA9_X_FP94BM: SummarizedExperiment with 395 rows and 16 columns 
#>  [2] 190321S_LCA10_X_FP97AG: SummarizedExperiment with 487 rows and 16 columns 
#>  [3] 190321S_LCA10_X_FP97_blank_01: SummarizedExperiment with 109 rows and 16 columns 
#>  [4] 190914S_LCB3_X_16plex_Set_21: SummarizedExperiment with 370 rows and 16 columns 
```
