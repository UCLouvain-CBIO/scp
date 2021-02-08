data("mqScpData")
data("sampleAnnotation")


test_that("readSCP: correct use", {
    ## Multiple batches
    scp <- readSCP(mqScpData, 
                   sampleAnnotation, 
                   batchCol = "Set", 
                   channelCol = "Channel")
    expect_identical(dims(scp)[, c("190222S_LCA9_X_FP94BM", "190321S_LCA10_X_FP97_blank_01", "190321S_LCA10_X_FP97AG", "190914S_LCB3_X_16plex_Set_21")],
                     matrix(c(290L, 109L, 375L, 323L, rep(16L, 4)), byrow = TRUE, nrow = 2, 
                            dimnames = list(NULL, c("190222S_LCA9_X_FP94BM", "190321S_LCA10_X_FP97_blank_01", "190321S_LCA10_X_FP97AG", "190914S_LCB3_X_16plex_Set_21"))))
    ## Note: there is something weird in this test... When executing the code locally, 
    ## the assays in `scp` are ordered alphabetically, but when runned by 
    ##  `devtools::test()`, the assays are ordered as they appear in `mqScpData`
    ## hence the need for `dims(scp)[, unique(mqScpData$Set)]` to avoid an error
    
    ## Make sure all rownames start with "PSM"
    expect_true(all(grepl("^PSM", unlist(rownames(scp)))))
    expect_true(all(grepl("^PSM", unlist(rownames(scp@ExperimentList@listData$`190222S_LCA9_X_FP94BM`)))))
    ## Make sur the column names are as expected
    expectedCols <- paste0(rep(unique(mqScpData$Set), 16), "_",
                           rep(paste0("RI", 1:16), each = 4))
    expect_true(all(unlist(colnames(scp)) %in% as.character(expectedCols)))
    
    ## Single batch
    scp <- readSCP(mqScpData %>% 
                     dplyr::filter(Set == "190222S_LCA9_X_FP94BM"), 
                   sampleAnnotation, batchCol = "Set", 
                   channelCol = "Channel")
    expect_identical(dims(scp),
                     matrix(c(290L, 16L), nrow = 2, 
                            dimnames = list(NULL, "190222S_LCA9_X_FP94BM")))
    ## Test remove empty columns
    scp <- readSCP(mqScpData, 
                   sampleAnnotation, 
                   batchCol = "Set", 
                   channelCol = "Channel",
                   removeEmptyCols = TRUE)
    expect_identical(dims(scp)[, c("190222S_LCA9_X_FP94BM", "190321S_LCA10_X_FP97_blank_01", "190321S_LCA10_X_FP97AG", "190914S_LCB3_X_16plex_Set_21")],
                     matrix(c(290L, 109L, 375L, 323L, rep(11L, 3), 16L),
                            byrow = TRUE, nrow = 2, 
                            dimnames = list(NULL, c("190222S_LCA9_X_FP94BM", "190321S_LCA10_X_FP97_blank_01", "190321S_LCA10_X_FP97AG", "190914S_LCB3_X_16plex_Set_21"))))
    
})

test_that("readSCP: warnings", {
  ## Missing batch in metadata 
  expect_warning(scp <- readSCP(mqScpData, 
                                sampleAnnotation  %>% 
                                  dplyr::filter(Set == "190222S_LCA9_X_FP94BM"), 
                                batchCol = "Set", 
                                channelCol = "Channel"),
                 regexp = "Missing metadata. The features are removed")
  expect_identical(dims(scp),
                   matrix(c(290L, 16L), nrow = 2, 
                          dimnames = list(NULL, "190222S_LCA9_X_FP94BM")))
  
})

test_that("readSingleCellExperiment: correct use", {
  sce <- readSingleCellExperiment(mqScpData, 
                                  ecol = grep("RI[0-9]*$",
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
                                  ecol = grep("RI[0-9]*$",
                                              colnames(badData)))
  expect_error(sce[1, ], regexp = "cannot have columns named \"seqnames\"")
})
