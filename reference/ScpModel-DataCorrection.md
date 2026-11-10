# Correct single-cell proteomics data

The function uses the data modelling output to generate corrected data
that can be used for downstream analysis. The input is expected to be a
`SummarizedExperiment` object that contains an estimated `ScpModel`.
There are two approaches:

- `scpKeepEffect()`: keep the effects of interests. The reconstructed
  data is the sum of the effect matrices for the variable of interest
  and the residuals. Note that the intercepts (baseline intensity of
  each feature) are not included by default, but they can be added when
  `intercept = TRUE`.

- `scpRemoveBatchEffect()`: remove any undesired effect. The batch
  corrected data is the input data minus the effect matrices that
  correspond to batch effect variables. Note that the intercepts
  (baseline intensity of each feature) are removed by default, but they
  can be kept when `intercept = FALSE`.

Despite the two approaches are conceptually different, they can lead to
similar results if the effects that are used to reconstruct the data are
the ones that are not removed when performing batch correction (see
examples).

The function returns a new `SummarizedExperiment` that contains an assay
with the batch corrected data. Note that the 'ScpModel\` is erased in
this new object.

## Usage

``` r
scpKeepEffect(object, effects = NULL, intercept = FALSE, name)

scpRemoveBatchEffect(object, effects = NULL, intercept = TRUE, name)
```

## Arguments

- object:

  An object that inherits from the `SummarizedExperiment` class. It must
  contain an estimated `ScpModel` in its metadata

- effects:

  A [`character()`](https://rdrr.io/r/base/character.html) vector. For
  `scpKeepEffect()`, which model variable should be used to reconstruct
  the data. For `scpRemoveBatchEffect()`, which model variable should be
  removed from the data. When `NULL` (default), both functions return
  the model residuals.

- intercept:

  A `logical(1)`. For `scpKeepEffect()`, should the intercepts be
  included when reconstructing the data? Defaults to `FALSE`, hence the
  intercepts are not included. For `scpRemoveBatchEffect()`, should the
  intercepts be removed from the data? Defaults to `TRUE`, hence the
  intercepts are removed from the data.

- name:

  A `character(1)` providing the name to use to retrieve the model
  results. When retrieving a model and `name` is missing, the name of
  the first model found in `object` is used.

## See also

- [ScpModel](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-class.md)
  for functions to extract information from the `ScpModel` object

- [ScpModel-Workflow](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-Workflow.md)
  to run a model on SCP data required for batch correction.

## Author

Christophe Vanderaa, Laurent Gatto

## Examples

``` r
data("leduc_minimal")
scpModelFormula(leduc_minimal)
#> ~1 + Channel + Set + MedianIntensity + SampleType

reconstructed <- scpKeepEffect(leduc_minimal, effects = "SampleType")
batchCorreced <- scpRemoveBatchEffect(
    leduc_minimal, effects = c("Channel", "Set", "MedianIntensity")
)
## The two approaches are identical
identical(reconstructed, batchCorreced)
#> [1] TRUE
```
