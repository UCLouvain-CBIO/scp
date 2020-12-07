data("scp1")

## Create a matrix with the following cases for testing CV computation 
## - 1 NA -> CV = NA
## - >1 NA -> CV = NA
## - 1 NA and 1 numeric -> CV = NA
## - 1 NA and >1 numeric -> CV = numeric
## - 1 numeric -> CV = NA
## - >1 numeric -> CV = numeric
group <- factor(c(2, rep(1, 2), rep(3, 3)))
m <- matrix(c(1, 1, NA, 1, 2, 3, rep(NA, 4), 2, 3), ncol = 2,
            dimnames = list(paste0(group, letters[1:6]), 
                            LETTERS[1:2]))


test_that(".rowCV", {
    ## Take missingness into account
    cvs <- .rowCV(m, group, na.rm = TRUE)
    expect_true(all(is.na(cvs[1:2, 1:2])))
    ## The lines below tests equality instead of identity because of 
    ## precision errors
    expect_equal(cvs[3, "A"], sd(m[4:6, "A"]) / mean(m[4:6, "A"]))
    expect_equal(cvs[3, "B"], sd(m[4:6, "B"], na.rm = TRUE) / 
                     mean(m[4:6, "B"], na.rm = TRUE))
    ## Test no ordering
    expect_identical(dimnames(cvs), 
                     list(as.character(unique(group)), colnames(m)))
    ## Test ordering
    expect_identical(dimnames(.rowCV(m, group, reorder = TRUE)), 
                     list(as.character(sort(unique(group))), colnames(m)))
    ## Test when na.rm = FALSE
    ## The output should be the same as above but the last element is NA
    expect_identical(c(as.vector(cvs[1:5]), NA), 
                     as.vector(.rowCV(m, group, na.rm = FALSE)))
    ## When nobs is too large, we expect only NAs
    expect_true(all(is.na(.rowCV(m, group, nobs = 5))))
    ## Test error: nobs is too small
    expect_error(.rowCV(m, group, nobs = 1), 
                 regexp = "No sd can be computed")
})
    

test_that("computeSCR", {
  ## Single assay
  test <- computeSCR(scp1, i = 1, colDataCol = "SampleType", 
                     samplePattern = "Reference", carrierPattern = "Carrier")
  expect_identical(rowData(test[[1]])$.meanSCR,
                   assay(scp1[[1]])[, 2] / assay(scp1[[1]])[, 1])
  ## Multiple assays
  test <- computeSCR(scp1, i = 1:2, colDataCol = "SampleType", 
                     samplePattern = "Reference", carrierPattern = "Carrier")
  expect_identical(rowData(test[[2]])$.meanSCR,
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
  ## Message: the PEP contains missing values
  rowData(scp1[[1]])$dart_PEP[1] <- NA
  expect_message(tmp <- computeFDR(scp1, i = 1:3, groupCol = "peptide", 
                                   pepCol = "dart_PEP"),
                 regexp = "missing values")
  expect_true(is.na(rowData(tmp[[1]])$dart_PEP[1]))
  ## Error: rowData variable not fount
  expect_error(computeFDR(scp1, i = 1, groupCol = "foo", pepCol = "dart_PEP"),
               regexp = paste0("not found in:\n", names(scp1)[1]))
  expect_error(computeFDR(scp1, i = 2, groupCol = "peptide", pepCol = "foo"),
               regexp = paste0("not found in:\n", names(scp1)[2]))
  ## Error: PEP must be a numeric between 0 and 1
  expect_error(computeFDR(scp1, i = 1, groupCol = "peptide", pepCol = "Length"),
               regexp = "is not a probability")
})

test_that("featureCV", {
    ## Create a SingleCellExperiment
    sce <- SingleCellExperiment(m, rowData = DataFrame(group = group))
    ## No normalization
    expect_identical(.rowCV(m, group = group, reorder = TRUE, 
                            na.rm = TRUE),
                     featureCV(x = sce, group = rowData(sce)$group, 
                               norm = "none", na.rm = TRUE))
    ## With normalization: divide columns by median
    expect_identical(.rowCV(sweep(m, 2, colMedians(m, na.rm = TRUE), "/"),
                            group = group, reorder = TRUE, na.rm = TRUE),
                     featureCV(x = sce, group = rowData(sce)$group, 
                               norm = "div.median", na.rm = TRUE))
    ## Double normalization: divide columns by median and divide rows
    ## by sum
    mproc <- sweep(m, 2, colMedians(m, na.rm = TRUE), "/")
    mproc <- sweep(mproc, 1, rowSums(mproc, na.rm = TRUE), "/")
    expect_identical(.rowCV(mproc, group = group, reorder = TRUE, 
                            na.rm = TRUE),
                     featureCV(x = sce, group = rowData(sce)$group, 
                               norm = c("div.median", "sum"), na.rm = TRUE))
})

test_that("medianCVperCell", {
    ## Check for single assay 
    scpfilt <- filterFeatures(scp1, ~ !is.na(Proteins))
    scp2 <- medianCVperCell(scpfilt, i = 1, groupBy = "Proteins")
    cvs <- colMedians(featureCV(scpfilt[[1]], group = rowData(scpfilt[[1]])$Proteins, 
                                nobs = 5, na.rm = TRUE), na.rm = TRUE)
    expect_identical(colData(scp2)[colnames(scp2)[[1]], "MedianCV"], cvs)
    expect_true(all(is.na(colData(scp2)[unlist(colnames(scp2)[2:3]), "MedianCV"])))
    ## Same for 2 assays
    scp2 <- medianCVperCell(scpfilt, i = 1:2, groupBy = "Proteins")
    cvs1 <- colMedians(featureCV(scpfilt[[1]], group = rowData(scpfilt[[1]])$Proteins, 
                                 nobs = 5, na.rm = TRUE), na.rm = TRUE)
    cvs2 <- colMedians(featureCV(scpfilt[[2]], group = rowData(scpfilt[[2]])$Proteins, 
                                 nobs = 5, na.rm = TRUE), na.rm = TRUE)
    expect_identical(colData(scp2)[unlist(colnames(scp2)[1:2]), "MedianCV"], 
                     c(cvs1, cvs2))
    expect_true(all(is.na(colData(scp2)[colnames(scp2)[[3]], "MedianCV"])))
    ## Warning: all computed median CV are NA (nobs is too high)
    expect_warning(medianCVperCell(scpfilt, i = 1, groupBy = "Proteins", nobs = 100),
                   regexp = "The median CV is NA for at least one column")
    ## Error: the colData name already exists
    expect_error(medianCVperCell(scpfilt, i = 1:5, groupBy = "Proteins",
                                 colData = "SampleType"),
                 regexp = "The colData name 'SampleType' already exists")
    ## Error: the assays contain duplicated samples
    expect_error(medianCVperCell(scp1, i = 1:5, groupBy = "Proteins"),
                 regexp = "Duplicated samples")
    
})

test_that("computeMedianCV_SCoPE2", {
    expect_warning(
        scp2 <- computeMedianCV_SCoPE2(scp1, i = "peptides", 
                                       peptideCol = "peptide", 
                                       proteinCol = "protein", 
                                       batchCol = "Set"),
        regexp = "deprecated")
    expect_equal(scp2[["peptides"]]$.MedianCV,
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

