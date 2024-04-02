
data("leduc_minimal")
leduc_minimal
## Overview of available cell annotations
colData(leduc_minimal)

####---- Model data ----####

f <- ~ 1 + ## intercept
    Channel + Set + ## batch variables
    MedianIntensity +## normalization
    SampleType ## biological variable
leduc_minimal <- scpModelWorkflow(leduc_minimal, formula = f)

####---- n/p feature filtering ----####

## Get n/p ratios
head(scpModelFilterNPRatio(leduc_minimal))

## Plot n/p ratios
scpModelFilterPlot(leduc_minimal)

## Change n/p ratio threshold
scpModelFilterThreshold(leduc_minimal) <- 2
scpModelFilterPlot(leduc_minimal)
