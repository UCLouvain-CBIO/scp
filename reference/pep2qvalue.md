# Compute q-values

This function computes q-values from the posterior error probabilities
(PEPs). The functions takes the PEPs from the given assay's `rowData`
and adds a new variable to it that contains the computed q-values.

## Usage

``` r
pep2qvalue(object, i, groupBy, PEP, rowDataName = "qvalue")
```

## Arguments

- object:

  A `QFeatures` object

- i:

  A [`numeric()`](https://rdrr.io/r/base/numeric.html) or
  [`character()`](https://rdrr.io/r/base/character.html) vector
  indicating from which assays the `rowData` should be taken.

- groupBy:

  A `character(1)` indicating the variable name in the `rowData` that
  contains the grouping variable, for instance to compute protein FDR.
  When `groupBy` is not missing, the best feature approach is used to
  compute the PEP per group, meaning that the smallest PEP is taken as
  the PEP of the group.

- PEP:

  A `character(1)` indicating the variable names in the `rowData` that
  contains the PEPs. Since, PEPs are probabilities, the variable must be
  contained in (0, 1).

- rowDataName:

  A `character(1)` giving the name of the new variable in the `rowData`
  where the computed FDRs will be stored. The name cannot already exist
  in any of the assay `rowData`.

## Value

A `QFeatures` object.

## Details

The q-value of a feature (PSM, peptide, protein) is the minimum FDR at
which that feature will be selected upon filtering (Savitski et al.). On
the other hand, the feature PEP is the probability that the feature is
wrongly matched and hence can be seen as a local FDR (Kall et al.).
While filtering on PEP is guaranteed to control for FDR, it is usually
too conservative. Therefore, we provide this function to convert PEP to
q-values.

We compute the q-value of a feature as the average of the PEPs
associated to PSMs that have equal or greater identification confidence
(so smaller PEP). See Kall et al. for a visual interpretation.

We also allow inference of q-values at higher level, for instance
computing the protein q-values from PSM PEP. This can be performed by
supplying the `groupBy` argument. In this case, we adopt the best
feature strategy that will take the best (smallest) PEP for each group
(Savitski et al.).

## References

Käll, Lukas, John D. Storey, Michael J. MacCoss, and William Stafford
Noble. 2008. “Posterior Error Probabilities and False Discovery Rates:
Two Sides of the Same Coin.” Journal of Proteome Research 7 (1): 40–44.

Savitski, Mikhail M., Mathias Wilhelm, Hannes Hahne, Bernhard Kuster,
and Marcus Bantscheff. 2015. “A Scalable Approach for Protein False
Discovery Rate Estimation in Large Proteomic Data Sets.” Molecular &
Cellular Proteomics: MCP 14 (9): 2394–2404.

## Examples

``` r
data("scp1")
scp1 <- pep2qvalue(scp1,
                   i = 1,
                   groupBy = "protein",
                   PEP = "dart_PEP",
                   rowDataName = "qvalue_protein")
## Check results
rowData(scp1)[[1]][, c("dart_PEP", "qvalue_protein")]
#> DataFrame with 166 rows and 2 columns
#>              dart_PEP qvalue_protein
#>             <numeric>      <numeric>
#> PSM3773   7.78683e-05    1.24198e-16
#> PSM9078   3.51840e-01    7.55061e-02
#> PSM9858   1.54857e-04    7.83117e-07
#> PSM11744  1.00000e+00    2.28967e-01
#> PSM21752  3.42850e-01    5.52325e-02
#> ...               ...            ...
#> PSM732069 4.89848e-05    4.41181e-09
#> PSM735396 3.89829e-07    1.24198e-16
#> PSM744756 1.54479e-06    7.82514e-08
#> PSM745037 4.44010e-01    1.14849e-01
#> PSM745130 2.07957e-04    1.83409e-10
```
