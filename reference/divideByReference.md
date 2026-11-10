# Divide assay columns by a reference column

The function divides the sample columns by a reference column. The
sample and reference columns are defined based on the provided `colvar`
variable and on regular expression matching.

## Usage

``` r
divideByReference(object, i, colvar, samplePattern = ".", refPattern)
```

## Arguments

- object:

  A `QFeatures` object

- i:

  A [`numeric()`](https://rdrr.io/r/base/numeric.html) or
  [`character()`](https://rdrr.io/r/base/character.html) vector
  indicating from which assays the `rowData` should be taken.

- colvar:

  A `character(1)` indicating the variable to take from
  `colData(object)` that gives the sample annotation.

- samplePattern:

  A `character(1)` pattern that matches the sample encoding in `colvar`.
  By default all samples are devided (using the regex wildcard `.`).

- refPattern:

  A `character(1)` pattern that matches the carrier encoding in
  `colvar`. Only one match per assay is allowed, otherwise only the
  first match is taken

## Value

A `QFeatures` object

## Details

The supplied assay(s) are replaced with the values computed after
reference division.

## Examples

``` r
data("scp1")
scp1 <- divideByReference(scp1, 
                          i = 1, 
                          colvar = "SampleType",
                          samplePattern = "Macrophage",
                          refPattern = "Ref")
                          
```
