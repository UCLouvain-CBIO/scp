data("mqScpData")
data("sampleAnnotation")

test_that("readSCP", {
    ## Most unit testing is performed in readQFeatures()
    quantCols <- paste0("Reporter.intensity.", 1:16)
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   quantCols = quantCols)
    expect_identical(scp@metadata$._type, "scp")
    expect_true(all(sapply(experiments(scp), inherits, "SummarizedExperiment")))
})

test_that("readSCP with sce", {
    ## Most unit testing is performed in readQFeatures()
    quantCols <- paste0("Reporter.intensity.", 1:16)
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   quantCols = quantCols,
                   experimentsAsSce = TRUE)
    expect_identical(scp@metadata$._type, "scp")
    expect_true(all(sapply(experiments(scp), inherits, "SingleCellExperiment")))
})

test_that("readSCPfromDIANN", {
    diannData <- read.delim(MsDataHub::Report.Derks2022.plexDIA.tsv())
    diannData$FileFile.Name <- diannData$Run
    scp <- readSCPfromDIANN(diannData, multiplexing = "mTRAQ")
    expect_identical(scp@metadata$._type, "scp")
    expect_true(all(sapply(experiments(scp), inherits, "SummarizedExperiment")))
})

test_that("readSCPfromDIANN with sce", {
    diannData <- read.delim(MsDataHub::Report.Derks2022.plexDIA.tsv())
    diannData$FileFile.Name <- diannData$Run
    scp <- readSCPfromDIANN(diannData,
                            multiplexing = "mTRAQ",
                            experimentsAsSce = TRUE)
    expect_identical(scp@metadata$._type, "scp")
    expect_true(all(sapply(experiments(scp), inherits, "SingleCellExperiment")))
})

test_that("readSingleCellExperiment", {
    ## Most unit testing is performed in readQFeatures()
    quantCols <- paste0("Reporter.intensity.", 1:16)
    sce <- readSingleCellExperiment(mqScpData,
                                    quantCols = grep("Reporter.intensity.[0-9]*$",
                                                     colnames(mqScpData)))
    expect_true(inherits(sce, "SingleCellExperiment"))
})
