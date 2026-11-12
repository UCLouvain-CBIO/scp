# Four metrics to report missing values

The function computes four metrics to report missing values in
single-cell proteomics.

## Usage

``` r
reportMissingValues(object, i, by = NULL)
```

## Arguments

- object:

  An object of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html).

- i:

  The index of the assay in `object`. The assay must contain an
  identification matrix, that is a matrix where an entry is `TRUE` if
  the value is observed and `FALSE` is the value is missing (see
  examples). `i` may be numeric, character or logical, but it must
  select only one assay.

- by:

  A vector of length equal to the number of columns in assay `i` that
  defines groups for which the metrics should be computed separately. If
  missing, the metrics are computed for the complete assay.

## Value

A `data.frame` with groups as rows and 5 columns:

- `LocalSensitivityMean`: the average number of features per cell.

- `LocalSensitivitySd`: the standard deviation of the local sensitivity.

- `TotalSensitivity`: the total number of features found in the dataset.

- `Completeness`: the proportion of values that are not missing in the
  data.

- `NumberCells`: the number of cells in the dataset.

## Examples

``` r

data("scp1")

## Define the identification matrix
peps <- scp1[["peptides"]]
assay(peps) <- !is.na(assay(peps))
scp1 <- addAssay(scp1, peps, "id")

## Report metrics
reportMissingValues(scp1, "id")
#>     LocalSensitivityMean LocalSensitivitySd TotalSensitivity Completeness
#> all             187.3684            23.2221              539    0.3476223
#>     NumberCells
#> all          38
## Report metrics by sample type
reportMissingValues(scp1, "id", scp1$SampleType)
#>            LocalSensitivityMean LocalSensitivitySd TotalSensitivity
#> Carrier                183.3333           26.85765              539
#> Reference              183.3333           26.85765              539
#> Unused                 191.0000           26.75818              539
#> Monocyte               185.6000           26.13044              539
#> Blank                  197.3333           28.86751              378
#> Macrophage             186.8000           23.00023              539
#>            Completeness NumberCells
#> Carrier       0.3401361           3
#> Reference     0.3401361           3
#> Unused        0.3543599           4
#> Monocyte      0.3443414           5
#> Blank         0.3661101           3
#> Macrophage    0.3465677          20

data
#> new("nonstandardGenericFunction", .Data = function (x, ...) 
#> {
#>     if (!missing(x)) {
#>         name <- substitute(x)
#>         if (is.name(name) && !exists(name, envir = parent.frame())) 
#>             x <- as.character(name)
#>     }
#>     standardGeneric("data")
#> }, generic = "data", package = "BiocGenerics", group = list(), 
#>     valueClass = character(0), signature = "x", default = NULL, 
#>     skeleton = (function (x, ...) 
#>     stop(gettextf("invalid call in method dispatch to '%s' (no default method)", 
#>         "data"), domain = NA))(x, ...))
#> <bytecode: 0x5c83715571b0>
#> <environment: 0x5c83715477d0>
#> attr(,"generic")
#> [1] "data"
#> attr(,"generic")attr(,"package")
#> [1] "BiocGenerics"
#> attr(,"package")
#> [1] "BiocGenerics"
#> attr(,"group")
#> list()
#> attr(,"valueClass")
#> character(0)
#> attr(,"signature")
#> [1] "x"
#> attr(,"default")
#> `\001NULL\001`
#> attr(,"skeleton")
#> (function (x, ...) 
#> stop(gettextf("invalid call in method dispatch to '%s' (no default method)", 
#>     "data"), domain = NA))(x, ...)
#> attr(,"class")
#> [1] "nonstandardGenericFunction"
#> attr(,"class")attr(,"package")
#> [1] "methods"
```
