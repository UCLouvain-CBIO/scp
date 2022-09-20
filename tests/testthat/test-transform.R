data("scp1")
scp1 <- updateObject(scp1)

test_that("divideByReference", {
    ## Correct use
    ## Single assay
    test <- divideByReference(scp1, i = "190321S_LCA10_X_FP97AG", 
                              colvar = "SampleType", samplePattern = ".",
                              refPattern = "Reference") 
    expect_identical(assay(test[[1]]), assay(scp1[[1]])/assay(scp1[[1]])[, 2])
    ## Multiple assays 
    test <- divideByReference(scp1, i = 1:3, colvar = "SampleType", 
                              samplePattern = ".", refPattern = "Reference") 
    expect_identical(assay(test[[2]]), assay(scp1[[2]])/assay(scp1[[2]])[, 2])
    ## Error: reference not found
    expect_error(divideByReference(scp1, i = 1:3, colvar = "SampleType", 
                                   samplePattern = ".", refPattern = "foo"),
                 regexp = "The reference pattern 'foo' did not match")
    ## Error: sample not found
    expect_error(divideByReference(scp1, i = 1:3, colvar = "SampleType", 
                                   samplePattern = "foo", refPattern = "Reference"),
                 regexp = "The sample pattern 'foo' did not match")
    ## Warning: multiple references found
    expect_warning(divideByReference(scp1, i = 1:3, colvar = "SampleType", 
                                     samplePattern = ".", refPattern = "."),
                   regexp = "Only the first match will be used")
})

test_that("function: normalizeSCP", {
  ## Test .normalizeSCP
  sce <- scp1[[1]]
  sce_norm <- .normalizeSCP(sce, method = "max")
  expect_identical(dim(sce), dim(sce_norm))
  e <- assay(sce) / rowMax(assay(sce))
  expect_identical(assay(sce_norm), e)
  ## Test normalizeSCP
  scp_norm <- normalizeSCP(scp1, 1, method = "max")
  expect_identical(scp_norm[["normAssay"]], sce_norm)
  ## Check the one-to-one link
  expect_identical(assayLink(scp_norm, "normAssay")@hits@elementMetadata$names_from,
                   assayLink(scp_norm, "normAssay")@hits@elementMetadata$names_to,
                   rownames(sce))
})

test_that("function: all normalize methods", {
    sce <- scp1[[1]]
    for (.method in MsCoreUtils::normalizeMethods()) {
        sce_norm <- .normalizeSCP(sce, method = .method)
        scp_norm <- normalizeSCP(scp1, 1, method = .method)
        expect_identical(sce_norm, scp_norm[["normAssay"]])
    }
})
