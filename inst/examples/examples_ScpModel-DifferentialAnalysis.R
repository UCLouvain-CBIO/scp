library("patchwork")
library("ggplot2")
data("leduc_minimal")
## Add n/p ratio information in rowData
rowData(leduc_minimal)$npRatio <- 
    scpModelFilterNPRatio(leduc_minimal, filtered = FALSE)

####---- Run differential abundance analysis ----####

(res <- scpDifferentialAnalysis(
    leduc_minimal, coefficients =  "MedianIntensity", 
    contrasts = list(c("SampleType", "Melanoma", "Monocyte"))
))
## IHW return a message because of the example data set has only few
## peptides, real dataset should not have that problem.

####---- Annotate results ----####

## Add peptide annotations available from the rowData
res <- scpAnnotateResults(
    res, rowData(leduc_minimal), 
    by = "feature", by2 = "Sequence"
)

####---- Plot results ----####

scpVolcanoPlot(res, textBy = "gene") |>
    wrap_plots(guides = "collect")

## Modify point and label aesthetics
scpVolcanoPlot(
    res, textBy = "gene", top = 20,
    pointParams = list(aes(colour = npRatio), alpha = 0.5),
    labelParams = list(size = 2, max.overlaps = 20)) |>
    wrap_plots(guides = "collect")

####---- Aggregate results ----####

## Aggregate to protein-level results
byProteinDA <- scpDifferentialAggregate(
    res, fcol = "Leading.razor.protein.id"
)
scpVolcanoPlot(byProteinDA) |>
    wrap_plots(guides = "collect")
