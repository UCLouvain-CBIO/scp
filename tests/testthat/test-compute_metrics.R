data("scp1")

test_that("computeSCR", {
  ## Single assay
  test <- computeSCR(scp1, i = 1, colDataCol = "SampleType", 
                     samplePattern = "Reference", carrierPattern = "Carrier")
  expect_identical(rowData(test[[1]])$meanSCR,
                   assay(scp1[[1]])[, 2] / assay(scp1[[1]])[, 1])
  ## Multiple assays
  test <- computeSCR(scp1, i = 1:2, colDataCol = "SampleType", 
                     samplePattern = "Reference", carrierPattern = "Carrier")
  expect_identical(rowData(test[[2]])$meanSCR,
                   assay(scp1[[2]])[, 2] / assay(scp1[[2]])[, 1])
  ## Warning: multiple match for carrier
  expect_warning(test2 <- computeSCR(scp1, i = 1:2, colDataCol = "SampleType", 
                                     samplePattern = "Reference", carrierPattern = "Carrier|Blank"),
                 regexp = "Multiple carriers found")
  expect_identical(test2, test)
  ## Error: colDataCol not fount
  expect_error(computeSCR(scp1, i = 1:2, colDataCol = "foo", 
                            samplePattern = "Reference", carrierPattern = "Carrier|Blank"),
                 regexp = "invalid names")
  ## Error: pattern not fount
  expect_error(computeSCR(scp1, i = 1:2, colDataCol = "SampleType", 
                          samplePattern = "foo", carrierPattern = "Carrier"),
                 regexp = "Pattern did not match")
})

test_that("computeFDR", {
  ## Correct use
  ## Compute FDR with scp function, then isolte the result for a single peptide
  ## and compute manually the FDR
  fdrFromPEP <- function(x) ## this is calc_fdr from SCoPE2
    return((cumsum(x[order(x)]) / seq_along(x))[order(order(x))])
  computeFDR(scp1, i = 1:3, groupCol = "peptide", pepCol = "dart_PEP") %>%
    rowDataToDF(1:3, vars = c("peptide", "dart_PEP", ".FDR")) %>%
    data.frame() %>%
    filter(peptide == "_AQLGGPEAAK_2") ->
    test
  expect_identical(test$.FDR, fdrFromPEP(test$dart_PEP))
  ## Error: rowData variable not fount
  expect_error(computeFDR(scp1, i = 1, groupCol = "foo", pepCol = "dart_PEP"),
               regexp = paste0("not found in:\n", names(scp1)[1]))
  expect_error(computeFDR(scp1, i = 2, groupCol = "peptide", pepCol = "foo"),
               regexp = paste0("not found in:\n", names(scp1)[2]))
  ## Error: PEP must be a numeric between 0 and 1
  expect_error(computeFDR(scp1, i = 1, groupCol = "peptide", pepCol = "Length"),
               regexp = "must link to a probability")
})

test_that("computeMedianCV", {
  warning("TODO")
})