data("mqScpData")
data("sampleAnnotation")

test_that("readSCP: correct use", {
    ## Multiple batches
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   channelCol = "Channel")
    expect_identical(sort(names(scp)), sort(unique(mqScpData$Raw.file)))
    expect_true(all(dims(scp)[2, ] == 16L))
    expect_true(sum(dims(scp)[1, ]) == nrow(mqScpData))
    ## Make sure all rownames start with "PSM"
    expect_true(all(grepl("^PSM", unlist(rownames(scp)))))
    ## Make sure the column names are as expected
    expectedCols <- paste0(rep(unique(mqScpData$Raw.file), 16),
                           rep(paste0("Reporter.intensity.", 1:16), each = 4))
    expect_identical(character(0), setdiff(unlist(colnames(scp)), as.character(expectedCols)))
    ## Single batch
    onebatch <- mqScpData %>%
        dplyr::filter(Raw.file == "190222S_LCA9_X_FP94BM")
    scp <- readSCP(onebatch,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   channelCol = "Channel")
    expect_identical(dims(scp)[1, ],
                     c("190222S_LCA9_X_FP94BM" = nrow(onebatch)))
    expect_identical(dims(scp)[2, ],
                     c("190222S_LCA9_X_FP94BM" = 16L))
    ## Test remove empty columns
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   channelCol = "Channel",
                   removeEmptyCols = TRUE)
    expect_identical(sort(names(scp)), sort(unique(mqScpData$Raw.file)))
    expect_true(all(dims(scp)[2, ] == c(rep(11, 3), 16)))
    expect_true(sum(dims(scp)[1, ]) == nrow(mqScpData))
    ## Test suffix
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   channelCol = "Channel",
                   suffix = paste0("_TMT", 1:16))
    expectedCols <- paste0(rep(unique(mqScpData$Raw.file), 16),
                           rep(paste0("_TMT", 1:16), each = 4))
    expectedCd <- DataFrame(sampleAnnotation)
    rownames(expectedCd) <- paste0(sampleAnnotation$Raw.file,
                                   sub("Repo.*ty.", "_TMT", sampleAnnotation$Channel))
    expect_true(all(unlist(colnames(scp)) %in% expectedCols))
    expect_identical(colData(scp)[sort(rownames(colData(scp))), ],
                     expectedCd[sort(rownames(expectedCd)), ])
    ## Test sep
    scp <- readSCP(mqScpData,
                   sampleAnnotation,
                   runCol = "Raw.file",
                   channelCol = "Channel",
                   suffix = paste0("TMT", 1:16),
                   sep = "_")
    expectedCols <- paste0(rep(unique(mqScpData$Raw.file), 16),
                           rep(paste0("_TMT", 1:16), each = 4))
    expectedCd <- DataFrame(sampleAnnotation)
    rownames(expectedCd) <- paste0(sampleAnnotation$Raw.file,
                                   sub("Repo.*ty.", "_TMT", sampleAnnotation$Channel))
    expect_true(all(unlist(colnames(scp)) %in% expectedCols))
    expect_identical(colData(scp)[sort(rownames(colData(scp))), ],
                     expectedCd[sort(rownames(expectedCd)), ])

})

test_that("readSCP: warnings", {
    ## Missing batch in metadata
    expect_warning(scp <- readSCP(mqScpData,
                                  dplyr::filter(sampleAnnotation,
                                                Raw.file == "190222S_LCA9_X_FP94BM"),
                                  runCol = "Raw.file",
                                  channelCol = "Channel"),
                   regexp = "Missing metadata. The features are removed")
    expect_identical(names(scp), "190222S_LCA9_X_FP94BM")
    expect_identical(dims(scp)[2, ], c("190222S_LCA9_X_FP94BM" = 16L))
    expect_identical(dims(scp)[1, ], c("190222S_LCA9_X_FP94BM" = sum(mqScpData$Raw.file == "190222S_LCA9_X_FP94BM")))

})

test_that("readSCP: error", {
    ## Suffix has not correct size
    expect_error(scp <- readSCP(mqScpData,
                                sampleAnnotation,
                                runCol = "Raw.file",
                                channelCol = "Channel",
                                suffix = (1:2)),
                 regex = "invalid rownames length")
    ## Suffix is not unique
    expect_error(expect_warning(
        scp <- readSCP(mqScpData,
                                sampleAnnotation,
                                runCol = "Raw.file",
                                channelCol = "Channel",
                                suffix = rep(1, 16)),
        regexp = "non-unique values"),
        regex = "duplicate 'row.names'")
})

test_that("readSingleCellExperiment: correct use", {
    sce <- readSingleCellExperiment(mqScpData,
                                    ecol = grep("Reporter.intensity.[0-9]*$",
                                                colnames(mqScpData)))
    ## Make sure dimensions are correct
    expect_identical(dim(sce),
                     c(nrow(mqScpData), 16L))
    ## Make sure class is correct
    expect_true(inherits(sce, "SingleCellExperiment"))
})

test_that("readSingleCellExperiment: error", {
    ## Some names in the rowData are not allowed by RangedSummarizedExperiment
    badData <- mqScpData
    badData$seqnames <- 1
    sce <- readSingleCellExperiment(badData,
                                    ecol = grep("Reporter.intensity.[0-9]*$",
                                                colnames(badData)))
    expect_error(sce[1, ], regexp = "cannot have columns named \"seqnames\"")
})
