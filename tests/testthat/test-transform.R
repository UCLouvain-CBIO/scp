data("scp1")

test_that(".replaceAssay", {
    ## Correct use
    mat <- mat2 <- matrix(1, nrow = nrow(scp1[[1]]), ncol = ncol(scp1[[1]]))
    dimnames(mat2) <- dimnames(scp1[[1]])
    se <- SummarizedExperiment(mat2)
    test <- .replaceAssay(scp1, i = 1, y = se)
    expect_identical(assay(test[[1]]), mat2)
    ## Error: cannot replace multiple assays
    expect_error(.replaceAssay(scp1, i = 1:2, y = list(se, se)),
                 regexp = "Only 1 assay can be replaced at a time.")
    ## Error: assay not a Summarized experiment
    expect_error(.replaceAssay(scp1, i = 1, y = mat),
                 regexp = "must inherits from a 'SummarizedExperiment' object")
    ## Error: dimnames do not match
    expect_error(.replaceAssay(scp1, i = 1, y = SummarizedExperiment(mat)),
                 regexp = "Colnames of old and new assays must match")
})


test_that("divideByReference", {
    ## Correct use
    ## Single assay
    test <- divideByReference(scp1, i = "190321S_LCA10_X_FP97AG", 
                              colDataCol = "SampleType", samplePattern = ".",
                              refPattern = "Reference") 
    expect_identical(assay(test[[1]]), assay(scp1[[1]])/assay(scp1[[1]])[, 2])
    ## Multiple assays 
    test <- divideByReference(scp1, i = 1:3, colDataCol = "SampleType", 
                              samplePattern = ".", refPattern = "Reference") 
    expect_identical(assay(test[[2]]), assay(scp1[[2]])/assay(scp1[[2]])[, 2])
    ## Error: reference not found
    expect_error(divideByReference(scp1, i = 1:3, colDataCol = "SampleType", 
                                   samplePattern = ".", refPattern = "foo"),
                 regexp = "The reference pattern 'foo' did not match")
    ## Error: sample not found
    expect_error(divideByReference(scp1, i = 1:3, colDataCol = "SampleType", 
                                   samplePattern = "foo", refPattern = "Reference"),
                 regexp = "The sample pattern 'foo' did not match")
    ## Warning: multiple references found
    expect_warning(divideByReference(scp1, i = 1:3, colDataCol = "SampleType", 
                                     samplePattern = ".", refPattern = "."),
                   regexp = "Only the first match will be used")
})

