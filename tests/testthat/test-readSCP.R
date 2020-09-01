data("mqScpData")
data("sampleAnnotation")


test_that("readSCP: correct use", {
    ## Multiple batches 
    scp <- readSCP(mqScpData, sampleAnnotation, batchCol = "Set", 
                   channelCol = "Channel")
    expect_identical(dims(scp),
                     matrix(c(334L, 433L, 321L, rep(16L, 3)), byrow = TRUE, nrow = 2, 
                            dimnames = list(NULL, c("190222S_LCA9_X_FP94BM", "190321S_LCA10_X_FP97AG", "190914S_LCB3_X_16plex_Set_21"))))
    ## Single batch
    scp <- readSCP(mqScpData %>% 
                     dplyr::filter(Set == "190222S_LCA9_X_FP94BM"), 
                   sampleAnnotation, batchCol = "Set", 
                   channelCol = "Channel")
    expect_identical(dims(scp),
                     matrix(c(334L, 16L), nrow = 2, 
                            dimnames = list(NULL, "190222S_LCA9_X_FP94BM")))
    
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
                   matrix(c(334L, 16L), nrow = 2, 
                          dimnames = list(NULL, "190222S_LCA9_X_FP94BM")))
  
})
