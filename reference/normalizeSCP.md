# Normalize single-cell proteomics (SCP) data

This function normalises an assay in a `QFeatures` according to the
supplied method (see Details). The normalized data is added as a new
assay

## Usage

``` r
normalizeSCP(object, i, name = "normAssay", method, ...)
```

## Arguments

- object:

  An object of class `QFeatures`.

- i:

  A numeric vector or a character vector giving the index or the name,
  respectively, of the assay(s) to be processed.

- name:

  A `character(1)` naming the new assay name. Defaults is are
  `normAssay`.

- method:

  `character(1)` defining the normalisation method to apply. See
  Details.\`

- ...:

  Additional parameters passed to
  [`MsCoreUtils::normalizeMethods()`](https://rdrr.io/pkg/MsCoreUtils/man/normalize.html).

## Value

A `QFeatures` object with an additional assay containing the normalized
data.

## Details

The `method` parameter in `normalize` can be one of `"sum"`, `"max"`,
`"center.mean"`, `"center.median"`, `"div.mean"`, `"div.median"`,
`"diff.meda"`, `"quantiles`", `"quantiles.robust`" or `"vsn"`. The
[`MsCoreUtils::normalizeMethods()`](https://rdrr.io/pkg/MsCoreUtils/man/normalize.html)
function returns a vector of available normalisation methods.

- For `"sum"` and `"max"`, each feature's intensity is divided by the
  maximum or the sum of the feature respectively. These two methods are
  applied along the features (rows).

- `"center.mean"` and `"center.median"` center the respective sample
  (column) intensities by subtracting the respective column means or
  medians. `"div.mean"` and `"div.median"` divide by the column means or
  medians. These are equivalent to `sweep`ing the column means (medians)
  along `MARGIN = 2` with `FUN = "-"` (for `"center.*"`) or `FUN = "/"`
  (for `"div.*"`).

- `"diff.median"` centers all samples (columns) so that they all match
  the grand median by subtracting the respective columns medians
  differences to the grand median.

- Using `"quantiles"` or `"quantiles.robust"` applies (robust) quantile
  normalisation, as implemented in
  [`preprocessCore::normalize.quantiles()`](https://rdrr.io/pkg/preprocessCore/man/normalize.quantiles.html)
  and
  [`preprocessCore::normalize.quantiles.robust()`](https://rdrr.io/pkg/preprocessCore/man/normalize.quantiles.robust.html).
  `"vsn"` uses the
  [`vsn::vsn2()`](https://huber-group-embl.github.io/vsn/reference/vsn2.html)
  function. Note that the latter also glog-transforms the intensities.
  See respective manuals for more details and function arguments.

For further details and examples about normalisation, see
[`MsCoreUtils::normalize_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/normalize.html).

## See also

[QFeatures::normalize](https://rdrr.io/pkg/BiocGenerics/man/normalize.html)
for more details about `normalize`

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
normalizeSCP(scp1, i = "proteins", name = "normproteins", 
             method = "center.mean")
#> An instance of class QFeatures (type: bulk) with 6 sets:
#> 
#>  [1] 190321S_LCA10_X_FP97AG: SummarizedExperiment with 166 rows and 11 columns 
#>  [2] 190222S_LCA9_X_FP94BM: SummarizedExperiment with 176 rows and 11 columns 
#>  [3] 190914S_LCB3_X_16plex_Set_21: SummarizedExperiment with 215 rows and 16 columns 
#>  [4] peptides: SummarizedExperiment with 539 rows and 38 columns 
#>  [5] proteins: SummarizedExperiment with 292 rows and 38 columns 
#>  [6] normproteins: SummarizedExperiment with 292 rows and 38 columns 
```
