library("patchwork")
library("ggplot2")
data("leduc_minimal")
leduc_minimal$cell <- rownames(colData(leduc_minimal))

####---- Run component analysis ----####

(pcs <- scpComponentAnalysis(
    leduc_minimal, method = "ASCA", effects = "SampleType",
    pcaFUN = "auto", residuals = FALSE, unmodelled = FALSE
))

####---- Annotate results ----####

## Add cell annotation available from the colData
bySamplePCs <- scpAnnotateResults(
    pcs$bySample, colData(leduc_minimal), by = "cell"
)

## Add peptide annotations available from the rowData
byFeaturePCs <- scpAnnotateResults(
    pcs$byFeature, rowData(leduc_minimal), 
    by = "feature", by2 = "Sequence"
)

####---- Plot results ----####

## Plot result in cell-space, ie each dot is a cell
scpComponentPlot(
    bySamplePCs, 
    pointParams = list( ## ggplot arguments
        aes(colour = SampleType, shape = lcbatch), 
        alpha = 0.6
    )
) |>
    wrap_plots(guides = "collect")

## Plot result in peptide-space, ie each dot is a peptide
scpComponentPlot(
    byFeaturePCs, 
    pointParams = list(colour = "dodgerblue", alpha = 0.6)
) |>
    wrap_plots(guides = "collect")

## Plot both
scpComponentBiplot(
    bySamplePCs, byFeaturePCs, 
    pointParams = list(aes(colour = SampleType), alpha = 0.6),
    labelParams = list(max.overlaps = 20),
    textBy = "gene"
) |>
    wrap_plots(guides = "collect")

####---- Aggregate results ----####

## Aggregate to protein-level results
byProteinPCs <- scpComponentAggregate(
    byFeaturePCs, fcol = "Leading.razor.protein.id"
)

## Plot result in protein-space, ie each dot is a protein
scpComponentPlot(
    byProteinPCs, 
    pointParams = list(colour = "firebrick", alpha = 0.6)
) |>
    wrap_plots(guides = "collect")
