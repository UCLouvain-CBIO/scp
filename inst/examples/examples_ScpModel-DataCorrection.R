data("leduc_minimal")
scpModelFormula(leduc_minimal)

reconstructed <- scpKeepEffect(leduc_minimal, effects = "SampleType")
batchCorreced <- scpRemoveBatchEffect(
    leduc_minimal, effects = c("Channel", "Set", "MedianIntensity")
)
## The two approaches are identical
identical(reconstructed, batchCorreced)
