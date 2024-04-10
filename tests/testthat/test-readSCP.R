data("mqScpData")
data("sampleAnnotation")

test_that("readSCP", {
    ## Most unit testing is performed in readQFeatures()
    quantCols <- paste0("Reporter.intensity.", 1:16)
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   quantCols = quantCols)
    expect_true(all(sapply(experiments(scp), inherits, "SingleCellExperiment")))
})

test_that("readSCPfromDIANN", {
    diannData <- read.delim(MsDataHub::Report.Derks2022.plexDIA.tsv())
    diannData$FileFile.Name <- diannData$Run
    scp <- readQFeaturesFromDIANN(diannData, multiplexing = "mTRAQ")
    expect_true(all(sapply(experiments(scp), inherits, "SingleCellExperiment")))
})

test_that("readSingleCellExperiment", {
    ## Most unit testing is performed in readQFeatures()
    quantCols <- paste0("Reporter.intensity.", 1:16)
    sce <- readSingleCellExperiment(mqScpData,
                                    ecol = grep("Reporter.intensity.[0-9]*$",
                                                colnames(mqScpData)))
    expect_true(inherits(sce, "SingleCellExperiment"))
})
