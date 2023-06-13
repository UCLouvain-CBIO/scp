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
    test <- computeSCR(scp1, i = 1, colvar = "SampleType",
                       samplePattern = "Reference", carrierPattern = "Carrier")
    expect_identical(rowData(test[[1]])$SCR,
                     unname(assay(scp1[[1]])[, 2] / assay(scp1[[1]])[, 1]))
    ## Multiple assays
    test <- computeSCR(scp1, i = 1:2, colvar = "SampleType",
                       samplePattern = "Reference", carrierPattern = "Carrier")
    expect_identical(rowData(test[[2]])$SCR,
                     unname(assay(scp1[[2]])[, 2] / assay(scp1[[2]])[, 1]))
    ## Multiple match for carrier: compute mean carrier
    test <- computeSCR(scp1, i = 1:2, colvar = "SampleType",
                        samplePattern = "Reference", carrierPattern = "Carrier|Blank",
                        carrierFUN = "mean")
    expect_identical(rowData(test[[2]])$SCR,
                     unname(assay(scp1[[2]])[, 2] / rowMeans(assay(scp1[[2]])[, c(1, 5)], na.rm = TRUE)))
    ## Multiple match for sample: compute median instead of mean, and
    ## test rowDataName
    test <- computeSCR(scp1, i = 1:2, colvar = "SampleType",
                       samplePattern = "Macrophage", carrierPattern = "Carrier|Blank",
                       carrierFUN = "mean", sampleFUN = "median",
                       rowDataName = "medianSCR")
    expect_identical(rowData(test[[2]])$medianSCR,
                     unname(rowMedians(assay(scp1[[2]])[, 7:11], na.rm = TRUE) /
                                rowMeans(assay(scp1[[2]])[, c(1, 5)], na.rm = TRUE)))
    ## Error: colvar not fount
    expect_error(computeSCR(scp1, i = 1:2, colvar = "foo",
                            samplePattern = "Reference", carrierPattern = "Carrier|Blank"),
                 regexp = "invalid names")
    ## Error: sample pattern not found
    expect_error(computeSCR(scp1, i = 1:2, colvar = "SampleType",
                            samplePattern = "foo", carrierPattern = "Carrier"),
                 regexp = "No match.*samplePattern")
    ## Error: carrier pattern not found
    expect_error(computeSCR(scp1, i = 1:2, colvar = "SampleType",
                            samplePattern = "Ref", carrierPattern = "foo"),
                 regexp = "No match.*carrierPattern")
    ## Error: the new rowData variable already exists
    expect_error(computeSCR(scp1, i = 1:2, colvar = "SampleType",
                            samplePattern = "Reference", carrierPattern = "Carrier",
                            rowDataName = "dart_PEP"),
                 regexp = "already exists")
})

test_that("pep2qvalue", {
    ## Correct use
    ## Compute q-values with scp function, then compare to the known result
    ## Test with groupBy
    test <- rowData(pep2qvalue(scp1,
                               i = 1:3,
                               groupBy = "protein",
                               PEP = "dart_PEP"))
    expect_equal(unique(test[[1]][test[[1]]$protein == "P61981", "qvalue"]),
                 3.104949e-17)
    ## Test missing groupBy
    expect_identical(rbindRowData(pep2qvalue(scp1, i = 1:3, PEP = "dart_PEP"), 1:3)$qvalue,
                     .pep2qvalue(rbindRowData(scp1, 1:3)$dart_PEP))
    ## Warning: the PEP contains missing values
    rowData(scp1[[1]])$dart_PEP[1] <- NA
    expect_warning(tmp <- pep2qvalue(scp1, i = 1:3, groupBy = "peptide",
                                     PEP = "dart_PEP"),
                   regexp = "no non-missing arguments to min")
    expect_true(is.na(rowData(tmp[[1]])$dart_PEP[1]))
    expect_true(all(!(is.na(rowData(tmp[[1]])$dart_PEP[-1]))))
    ## Error: rowData variable not found
    expect_error(pep2qvalue(scp1, i = 1, groupBy = "foo", PEP = "dart_PEP"),
                 regexp = "not found")
    expect_error(pep2qvalue(scp1, i = 2, groupBy = "peptide", PEP = "foo"),
                 regexp = "not found")
    ## Error: PEP must be a numeric between 0 and 1
    expect_error(pep2qvalue(scp1, i = 1, groupBy = "peptide", PEP = "Length"),
                 regexp = "is not a probability")
    ## Error: the new rowData variable already exists
    expect_error(pep2qvalue(scp1, i = 1, groupBy = "peptide", PEP = "Length",
                            rowDataName = "dart_PEP"),
                 regexp = "already exists")
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
    ## SCoPE2 normalization: the normalization factor is computed by
    ## dividing columns by the median and then dividing rows by the mean.
    nf <- sweep(m, 2, colMedians(m, na.rm = TRUE), "/")
    nf <- rowMeans(nf, na.rm = TRUE)
    expect_identical(.rowCV(sweep(m, 1, nf, "/"), group = group,
                            reorder = TRUE, na.rm = TRUE),
                     featureCV(x = sce, group = rowData(sce)$group,
                               norm = "SCoPE2", na.rm = TRUE))
})

test_that("medianCVperCell", {
    ## Check for single assay
    scpfilt <- filterFeatures(scp1, ~ !is.na(Proteins))
    scp2 <- medianCVperCell(scpfilt, i = 1, groupBy = "Proteins")
    cvs <- colMedians(featureCV(scpfilt[[1]], group = rowData(scpfilt[[1]])$Proteins,
                                nobs = 5, na.rm = TRUE),
                      na.rm = TRUE, useNames = FALSE)
    expect_identical(colData(scp2)[colnames(scp2)[[1]], "MedianCV"], cvs)
    expect_true(all(is.na(colData(scp2)[unlist(colnames(scp2)[2:3]), "MedianCV"])))
    ## Same for 2 assays
    scp2 <- medianCVperCell(scpfilt, i = 1:2, groupBy = "Proteins")
    cvs1 <- colMedians(featureCV(scpfilt[[1]], group = rowData(scpfilt[[1]])$Proteins,
                                 nobs = 5, na.rm = TRUE),
                       na.rm = TRUE, useNames = FALSE)
    cvs2 <- colMedians(featureCV(scpfilt[[2]], group = rowData(scpfilt[[2]])$Proteins,
                                 nobs = 5, na.rm = TRUE),
                       na.rm = TRUE, useNames = FALSE)
    expect_identical(colData(scp2)[unlist(colnames(scp2)[1:2]), "MedianCV"],
                     c(cvs1, cvs2))
    expect_true(all(is.na(colData(scp2)[colnames(scp2)[[3]], "MedianCV"])))
    ## Warning: all computed median CV are NA (nobs is too high)
    expect_warning(medianCVperCell(scpfilt, i = 1, groupBy = "Proteins", nobs = 100),
                   regexp = "The median CV could not be computed for one or more")
    ## Error: the colData name already exists
    expect_error(medianCVperCell(scpfilt, i = 1:5, groupBy = "Proteins",
                                 colData = "SampleType"),
                 regexp = "The colData name 'SampleType' already exists")
    ## Error: the assays contain duplicated samples
    expect_error(medianCVperCell(scp1, i = 1:5, groupBy = "Proteins"),
                 regexp = "Duplicated samples")

})
