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
               regexp = "is not a probability")
})

test_that("computeMedianCV", {
  ## TODO improve this test after improving the CV computation algorithm
  expect_equal(computeMedianCV(scp1, i = "peptides", peptideCol = "peptide", 
                                   proteinCol = "protein", 
                                   batchCol = "Set")[["peptides"]]$MedianCV,
                   c(0.645293098683498, 0.773549942815159, 0.663774239146459, 
                     0.72830710905132, 0.778688251750191, 1.03956742480134, 
                     0.526933183452957, 0.447329625349529, 0.346740321446223, 
                     0.441428083348158, 0.774968411920989, 0.66756262970139,
                     1.0564230469663, 0.714686435220506, 0.977975517652186, 
                     1.28733510164692, 0.948691932585763, NA, 0.908961764969587, 
                     1.07509134173885, 0.71998288423204, 0.947756691696601,
                     0.835782730630332, 0.929717820468142, 1.0011663228458, NA, 
                     0.808254681838227, 1.09581518003356, 1.1330059279254,
                     1.12084088292072, 0.806942496031636, 1.19054653629395, 
                     1.1396089190922, NA, NA, NA, NA, NA))
})