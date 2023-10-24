data("leduc_minimal")

####---- Run analysis of variance ----####

(var <- scpVarianceAnalysis(leduc_minimal))

####---- Annotate results ----####

## Add peptide annotations available from the rowData
var <- scpAnnotateResults(
    var, rowData(leduc_minimal), by = "feature", by2 = "Sequence"
)

####---- Plot results ----####

## Plot the analysis of variance through the whole data
scpVariancePlot(var)

## Plot the analysis of variance for the top 20 peptides with highest
## percentage of variance explained by the cell type
scpVariancePlot(
    var, effect = "SampleType", top = 20, combined = FALSE
)

## Same but grouped by protein
scpVariancePlot(
    var, effect = "SampleType", top = 20, combined = FALSE, fcol = "gene"
)

####---- Aggregate results ----####

## Aggregate to protein-level results
varProtein <- scpVarianceAggregate(var, fcol = "gene")
scpVariancePlot(
    varProtein, effect = "SampleType", top = 20, combined = FALSE
)
