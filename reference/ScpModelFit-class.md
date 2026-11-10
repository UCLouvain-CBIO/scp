# Class to store the components of an estimated model for a feature

An `ScpModelFit` object is expected to be stored as a list element in
the `scpModelFitList` of an `ScpModel` object. The `ScpModelFit` object
should **never be accessed directly** by the user. Refer to the
[ScpModel](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-class.md)
for a list of function to access the information in an `ScpModelFit`.
The `ScpModelFit` class contains several slots that contain the model
output for a feature:

- `n`: an `integer`, the number of observations for the feature

- `p`: an `integer`, the number of coefficient to estimate

- `coefficients`: a `numeric` vector with the estimated coefficients

- `residuals`: a `numeric` vector with the estimated residuals

- `effects`: a `List` with the

- `df`: an `integer` providing the number of degrees of freedom of the
  model estimation

- `var`: a `numeric` vector with the residual variance of the model
  estimation

- `uvcov`: the unscaled variance covariance `matrix`

- `levels`: a named `List` where each elements corresponds to a
  categorical model variable and contains a vector with the possible
  categories.

## See also

[ScpModel](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-class.md)
for a description of the class that relies on `ScpModelFit`

## Author

Christophe Vanderaa, Laurent Gatto

## Examples

``` r
new("ScpModelFit") ## this should never be used by the user
#> An object of class "ScpModelFit"
#> Slot "n":
#> integer(0)
#> 
#> Slot "p":
#> integer(0)
#> 
#> Slot "coefficients":
#> numeric(0)
#> 
#> Slot "residuals":
#> numeric(0)
#> 
#> Slot "effects":
#> <S4 Type Object>
#> attr(,"elementType")
#> [1] "ANY"
#> attr(,"elementMetadata")
#> `\001NULL\001`
#> attr(,"metadata")
#> list()
#> 
#> Slot "df":
#> numeric(0)
#> 
#> Slot "var":
#> numeric(0)
#> 
#> Slot "uvcov":
#> <0 x 0 matrix>
#> 
#> Slot "levels":
#> <S4 Type Object>
#> attr(,"elementType")
#> [1] "ANY"
#> attr(,"elementMetadata")
#> `\001NULL\001`
#> attr(,"metadata")
#> list()
#> 
```
