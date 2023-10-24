data("leduc_minimal")
var <- scpVarianceAnalysis(leduc_minimal)
colnames(var$Residuals)
## Add peptide annotations available from the rowData
var <- scpAnnotateResults(
    var, rowData(leduc_minimal), by = "feature", by2 = "Sequence"
)
colnames(var$Residuals)
