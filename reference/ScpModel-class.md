# Class to store the results of single-cell proteomics modelling

An `ScpModel` object must be always stored in the `metadata()` of an
object that inherits from the `SummarizedExperiment` class. The
`ScpModel` object should **never be accessed directly** by the user.
Instead, we provide several setter function to retrieve information that
may be useful to the user.The `ScpModel` class contains several slots:

- `scpModelFormula`: a `formula` object controlling which variables are
  to be modelled.

- `scpModelInputIndex`: a `numeric(1)`, selecting the assay to use in
  the `SummarizedExperiment` object as input matrix. Note that this slot
  serves as a pointer, meaning that the quantitative data is not
  duplicated. Any change to the assay in the `SummarizedExperiment` will
  impact the estimation of the `ScpModel` object.

- `scpModelFilterThreshold`: A `numeric(1)` indicating the minimal n/p
  ratio required for a feature to be included in further model
  exploration. n/p is the number of measured values for a features
  divided by the number of coefficients to estimate. n/p cannot be
  smaller than 1 because this would lead to over-specified models.

- `scpModelFitList`: A `List` that contains the model results for each
  feature. Each element is a `ScpModelFit` object (see
  [`ScpModelFit`](https://UCLouvain-CBIO.github.io/scp/reference/ScpModelFit-class.md))

## Usage

``` r
scpModelFormula(object, name)

scpModelInput(object, name, filtered = TRUE)

scpModelFilterThreshold(object, name)

scpModelFilterNPRatio(object, name, filtered = TRUE)

scpModelResiduals(object, name, join = TRUE, filtered = TRUE)

scpModelEffects(object, name, join = TRUE, filtered = TRUE)

scpModelNames(object)

scpModelFilterThreshold(object, name) <- value
```

## Arguments

- object:

  An object that inherits from the `SummarizedExperiment` class.

- name:

  A `character(1)` providing the name to use to store or retrieve the
  modelling results. When retrieving a model and `name` is missing, the
  name of the first model found in `object` is used.

- filtered:

  A `logical(1)` indicating whether the output should return all
  features (`FALSE`) or the features that comply to the n/p ratio
  threshold (`TRUE`).

- join:

  A `logical(1)` indicating whether the output should be combined in a
  single matrix (`TRUE`) or it should be returned as a list with one
  element for each feature (`FALSE`). When `TRUE`, any gaps across
  features will be filled with NA's.

- value:

  An `numeric(1)`, the new value for the n/p ratio threshold

## Getters

Each slot has a getter function associated:

- `scpModelNames()`: returns a vector of names of `ScpModel` objects
  stored in the `SummarizedExperiment` object.

- `scpModelFormula()`: returns the `formula` slot of the `ScpModel`
  within an object that inherits from the `SummarizedExperiment` class.

- `scpModelFilterThreshold()`: returns the n/p ration threshold used for
  feature filtering.

- `scpModelInput()`: returns a `matrix` with the quantitative values
  used as input of the model. Hence, the matrix contains the data before
  modelling. If `filtered = TRUE`, the feature of the matrix are
  restricted to the features that satisfy the n/p ratio threshold.

- `scpModelFilterNPRatio()`: returns the computed n/p ratio for each
  feature. If `filtered = TRUE`, the function returns only the n/p of
  the features that satisfy the n/p ratio threshold.

- `scpModelResiduals()`: when `join = FALSE`, the function returns a
  list where each element corresponds to a feature and contains the
  estimated residuals. When `join = TRUE` (default), the function
  combines the list into a matrix with features in rows and cells in
  columns, and filling the gaps with `NA`. If `filtered = TRUE`, the
  feature of the matrix are restricted to the features that satisfy the
  n/p ratio threshold.

- `scpModelEffects()`: when `join = FALSE`, the function return a list
  where each element of the list corresponds to a feature. Each element
  contains another list with as many elements as variable in the model
  and each element contains the data effect vector for that vector. When
  `join = TRUE` (default), each element of the list is a matrix with
  features in rows and cells in columns where gaps are filled with `NA`.
  If `filtered = TRUE`, the feature of the matrix are restricted to the
  features that satisfy the n/p ratio threshold.

Setter:

- `scpModelFilterThreshold<-()`: the function changes the n/p ratio
  threshold used for filtering features.

## See also

- [ScpModelFit](https://UCLouvain-CBIO.github.io/scp/reference/ScpModelFit-class.md)
  for a description of the class that store modelling results

- [ScpModel-Workflow](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-Workflow.md)
  that uses the class to store the estimated model.

## Author

Christophe Vanderaa, Laurent Gatto

## Examples

``` r

data("leduc_minimal")

####---- Getters ----####

scpModelNames(leduc_minimal)
#> [1] "model"

scpModelFormula(leduc_minimal)
#> ~1 + Channel + Set + MedianIntensity + SampleType

dim(leduc_minimal)
#> [1] 200  73
dim(scpModelInput(leduc_minimal))
#> [1] 140  73
dim(scpModelInput(leduc_minimal, filtered = FALSE))
#> [1] 200  73

head(scpModelFilterNPRatio(leduc_minimal))
#>      SAVEDEGLK       APNVVVTR      IVVVTAGVR GFQEVVTPNIFNSR ENAYDLEANLAVLK 
#>       1.277778       3.476190       1.411765       2.333333       1.235294 
#>      IGPLGLSPK 
#>       1.888889 

dim(scpModelResiduals(leduc_minimal))
#> [1] 140  73
dim(scpModelResiduals(leduc_minimal, filtered = FALSE))
#> [1] 200  73
scpModelResiduals(leduc_minimal, join = FALSE)
#> NumericList of length 140
#> [["SAVEDEGLK"]] wAL00284RI5=0.699508268833165 ...
#> [["APNVVVTR"]] eAL00219RI5=-0.211268540934462 ...
#> [["IVVVTAGVR"]] wAL00284RI5=0.103674283809887 ...
#> [["GFQEVVTPNIFNSR"]] eAL00219RI5=0.637577960221826 ...
#> [["ENAYDLEANLAVLK"]] eAL00220RI5=0.340761402427571 ...
#> [["IGPLGLSPK"]] eAL00220RI5=-0.111917031828062 ...
#> [["VAETANEEEVKK"]] eAL00219RI5=-0.0862785009338332 ...
#> [["TADDPSLSLIK"]] eAL00220RI5=0.0148848162385562 ...
#> [["GPPPSYGGSSR"]] eAL00219RI5=-0.117374478914073 ...
#> [["IDATSASVLASR"]] eAL00219RI5=-0.145390183837019 ...
#> ...
#> <130 more elements>

scpModelEffects(leduc_minimal)
#> List of length 4
#> names(4): Channel Set MedianIntensity SampleType
dim(scpModelEffects(leduc_minimal)$Set)
#> [1] 140  73
dim(scpModelEffects(leduc_minimal, filtered = FALSE)$Set)
#> [1] 200  73
scpModelEffects(leduc_minimal, join = FALSE)[[1]]
#> NumericList of length 4
#> [["Channel"]] wAL00284RI5=-0.0371719605222012 ...
#> [["Set"]] wAL00284RI5=0.0994019567725992 ... wAL00286RI17=0.249998609880628
#> [["MedianIntensity"]] wAL00284RI5=0.119366566488542 ...
#> [["SampleType"]] wAL00284RI5=-0.382965748230049 ...

scpModelFilterThreshold(leduc_minimal)
#> [1] 1

####---- Setter ----####

scpModelFilterThreshold(leduc_minimal) <- 2
scpModelFilterThreshold(leduc_minimal)
#> [1] 2
```
