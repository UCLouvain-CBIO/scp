# Compute the sample over carrier ratio (SCR)

The function computes the ratio of the intensities of sample channels
over the intentisty of the carrier channel for each feature. The ratios
are averaged within the assay.

## Usage

``` r
computeSCR(
  object,
  i,
  colvar,
  samplePattern,
  sampleFUN = "mean",
  carrierPattern,
  carrierFUN = sampleFUN,
  rowDataName = "SCR"
)
```

## Arguments

- object:

  A `QFeatures` object.

- i:

  A [`character()`](https://rdrr.io/r/base/character.html) or
  [`integer()`](https://rdrr.io/r/base/integer.html) indicating for
  which assay(s) the SCR needs to be computed.

- colvar:

  A `character(1)` indicating the variable to take from
  `colData(object)` that gives the sample annotation.

- samplePattern:

  A `character(1)` pattern that matches the sample encoding in `colvar`.

- sampleFUN:

  A `character(1)` or `function` that provides the summarization
  function to use (eg mean, sum, media, max, ...). Only used when the
  pattern matches multiple samples. Default is `mean`. Note for custom
  function, `na.rm = TRUE` is passed to `sampleFUN` to ignore missing
  values, make sure to provide a function that accepts this argument.

- carrierPattern:

  A `character(1)` pattern that matches the carrier encoding in
  `colvar`. Only one match per assay is allowed, otherwise only the
  first match is taken

- carrierFUN:

  A `character(1)` or `function` that provides the summarization
  function to use (eg mean, sum, media, max, ...). Only used when the
  pattern matches multiple carriers. Default is the same function as
  `sampleFUN`. Note for custom function, `na.rm = TRUE` is passed to
  `carrierFUN` to ignore missing values, make sure to provide a function
  that accepts this argument.

- rowDataName:

  A `character(1)` giving the name of the new variable in the `rowData`
  where the computed SCR will be stored. The name cannot already exist
  in any of the assay `rowData`.

## Value

A `QFeatures` object for which the `rowData` of the given assay(s) is
augmented with the mean SCR.

## Examples

``` r
data("scp1")
scp1 <- computeSCR(scp1,
                   i = 1,
                   colvar = "SampleType",
                   carrierPattern = "Carrier",
                   samplePattern = "Blank|Macrophage|Monocyte",
                   sampleFUN = "mean",
                   rowDataName = "MeanSCR")
## Check results
rowData(scp1)[[1]][, "MeanSCR"]
#>   [1] 0.015510472 0.017444270 0.013998780 0.032842883         NaN 0.019183768
#>   [7] 0.035322014 0.011685058 0.018249260 0.027776749 0.017919113 0.011189229
#>  [13] 0.006934865 0.042799880 0.019657274 0.000000000 0.022816835 0.045264655
#>  [19] 0.026838068 0.006634830 0.008038560 0.007882511 0.007558454 0.022777305
#>  [25] 0.021362267 0.942637117 0.000000000 0.047027150 0.110357190 0.011164589
#>  [31] 0.027721278 0.047673647 0.011990475 0.016613643 0.015244373 0.009326254
#>  [37] 0.000000000 0.022829389 0.034530363 0.007661202 0.025641138 0.038001092
#>  [43] 0.018844896 0.008099232 0.097323303 0.007196344 0.067520833 0.007050806
#>  [49] 0.018265529 0.028173254 0.012425863 0.246754068 0.098868861 0.056916766
#>  [55] 0.972946418 0.077164878 0.011155643 0.014888172 0.015396338 0.023861222
#>  [61] 0.016795817 0.013005209 0.014811532 0.019834972 0.023055615 0.033241575
#>  [67] 0.007923177 0.016450749         NaN         NaN 0.154418700 0.166120958
#>  [73] 0.024312229 0.028761939 0.021799714 0.016790123         NaN         NaN
#>  [79] 0.032421331 0.023590284 0.009421714 0.008141877 0.209825122         NaN
#>  [85] 0.056103552 0.016963036 0.007254255 0.015776270 0.026212478 0.008341729
#>  [91] 0.013683271 0.127146548 0.262428332 0.010815081         NaN 0.021939122
#>  [97] 0.010556431 0.011106046         NaN         NaN 1.421039985 0.290833213
#> [103] 0.006530733 0.014220928 0.020974690 0.111953099 0.018738790 0.014850186
#> [109] 0.211281919 0.007306382 0.058905617 0.024971824 0.102530538 0.020856375
#> [115] 0.042998307 0.020087350 0.051781427 0.011624164 0.041935595 0.011965542
#> [121] 0.044674680 0.017633856 0.054716089 0.447498978 0.032647408 2.315250048
#> [127] 1.186082092 0.391160056 0.061278850         NaN 0.021773007 0.020502310
#> [133]         NaN 0.090177289 0.056202829         NaN         NaN         NaN
#> [139] 0.019264467 0.011998890         NaN 0.016139295 0.012450886 0.029591602
#> [145] 0.013236440 0.038606232 0.025146765         NaN 0.010298813 0.076093736
#> [151] 0.025677348 0.021617243 0.597881028         NaN 0.930486819 0.030299867
#> [157] 0.002346206 0.024708650 0.012180817 0.019880893 0.011504703 0.022299321
#> [163] 0.030614294 0.006889869 0.120408305 0.037997229
```
