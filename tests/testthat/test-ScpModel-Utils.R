
## Internal function that creates a minimal SE object as expected by
## scplainer for unit testing ScpModel class methods
## @param nr Number of rows
## @param nc Number of columns
.createMinimalData <- function(nr = 10, nc = 5) {
    require("SummarizedExperiment")
    a <- matrix(1, nr, nc)
    rownames(a) <- letters[1:nr]
    colnames(a) <- LETTERS[1:nc]
    SummarizedExperiment(assays = List(assay = a))
}

test_that("scpAnnotateResults", {
    library("S4Vectors")
    ## tableList is not a list = error
    expect_error(
        scpAnnotateResults(
            DataFrame(foo = "bar"), DataFrame(foo = "bar"), by = "foo"
        ),
        "'tableList' must be a list of 'DFrame' or 'data.frame'."
    )
    ## tableList elements are not tables = error
    expect_error(
        scpAnnotateResults(
            List(foo = "bar"), DataFrame(foo = "bar"), by = "foo"
        ),
        "'tableList' must be a list of 'DFrame' or 'data.frame'."
    )
    ## tableList is empty = error
    expect_error(
        scpAnnotateResults(
            List(), DataFrame(foo = "bar"), by = "foo"
        ),
        "'tableList' must be a list of 'DFrame' or 'data.frame'."
    )
    ## by not found in table list = error
    d1 <- exp <- DataFrame(id = 1:10, var1 = paste("value1", 1:10))
    d2 <- DataFrame(id = 1:10, var2 = paste("value2", 1:10))
    expect_error(
        scpAnnotateResults(List(test = d1), d2, by = "foo"),
        "'foo' not found in the columns of 'tableList' elements."
    )
    ## by2 not found in annotations = error
    expect_error(
        scpAnnotateResults(List(test = d1), d2, by = "id", by2 = "foo"),
        "'foo' not found in 'annotations'."
    )
    ## New table column = added
    exp$var2 <- d2$var2
    expect_identical(
        scpAnnotateResults(List(test = d1), d2, by = "id"),
        List(test = exp)
    )
    ## Existing table column = overwrite
    d1$var2 <- paste("otherValue", 1:10)
    expect_identical(
        scpAnnotateResults(List(test = d1), d2, by = "id"),
        List(test = exp)
    )
    ## Works for multiple tables
    expect_identical(
        scpAnnotateResults(List(test1 = d1, test2 = d1), d2, by = "id"),
        List(test1 = exp, test2 = exp)
    )
    ## Different merging columns need 'by' and 'by2'
    d2 <- DataFrame(id2 = 1:10, var2 = paste("value2", 1:10))
    expect_identical(
        scpAnnotateResults(List(test1 = d1, test2 = d1), d2, by = "id", by2 = "id2"),
        List(test1 = exp, test2 = exp)
    )
    ## Test different ordering
    d2 <- d2[10:1, ]
    expect_identical(
        scpAnnotateResults(List(test1 = d1, test2 = d1), d2, by = "id", by2 = "id2"),
        List(test1 = exp, test2 = exp)
    )
    ## Test subset of annotations
    exp$var2[1:5] <- NA
    expect_identical(
        scpAnnotateResults(
            List(test1 = d1, test2 = d1), d2[1:5, ], by = "id", by2 = "id2"
        ),
        List(test1 = exp, test2 = exp)
    )
    ## Test subset of table list
    exp <- DataFrame(id = 1:5, var1 = paste("value1", 1:5), var2 = paste("value2", 1:5))
    expect_identical(
        scpAnnotateResults(
            List(test1 = d1[1:5, ], test2 = d1[1:5, ]), d2, by = "id", by2 = "id2"
        ),
        List(test1 = exp, test2 = exp)
    )
})

test_that("scpKeepEffect", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    ## Effect not modelled = error
    expect_error(
        scpKeepEffect(se, effects = "foo"),
        "effects %in% scpModelEffectNames.object, name. is not TRUE"
    )
    expect_error(
        scpKeepEffect(se, effects = c("foo", "condition")),
        "effects %in% scpModelEffectNames.object, name. are not all TRUE"
    )
    expect_error(
        scpKeepEffect(se, effects = c("foo", "foo2")),
        "effects %in% scpModelEffectNames.object, name. are not all TRUE"
    )
    ## No effect = return residuals
    exp <- se[c("a", "c", "e", "h", "i"), ]
    assays(exp) <- list(scpModelResiduals(se))
    metadata(exp) <- list()
    expect_identical(
        scpKeepEffect(se),
        exp
    )
    ## Keep one effect
    assay(exp) <- assay(exp) + scpModelEffects(se)[["condition"]]
    expect_identical(
        scpKeepEffect(se, effects = "condition"),
        exp
    )
    ## Keep more effects
    assay(exp) <- assay(exp) + scpModelEffects(se)[["batch"]]
    expect_identical(
        scpKeepEffect(se, effects = c("condition", "batch")),
        exp
    )
    ## keep intercept
    assay(exp) <- assay(exp) + scpModelIntercept(se)
    expect_identical(
        scpKeepEffect(se, effects = c("condition", "batch"), intercept = TRUE),
        exp
    )
    ## batch correct from second model
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch, name = "model2")
    expect_identical(
        scpKeepEffect(se, effects = c("condition", "batch"),
                      intercept = TRUE, name = "model2"),
        exp
    )
    ## only scp models are removed from metadata
    metadata(se)$foo <- metadata(exp)$foo <- "bar"
    expect_identical(
        metadata(scpKeepEffect(se, effects = c("condition", "batch"),
                               intercept = TRUE, name = "model2")),
        list(foo = "bar")
    )
    ## when SCE, empty the reducedDims slot
    require(SingleCellExperiment)
    sce <- as(se, "SingleCellExperiment")
    reducedDim(sce, "PCA") <- t(assay(sce))
    exp <- List()
    names(exp) <- character()
    expect_identical(
        reducedDims(scpKeepEffect(sce, effects = c("condition", "batch"),
                                  intercept = TRUE, name = "model2")),
        exp
    )
})

test_that("scpRemoveBatchEffect", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    ## Effect not modelled = error
    expect_error(
        scpRemoveBatchEffect(se, effects = "foo"),
        "effects %in% scpModelEffectNames.object, name. is not TRUE"
    )
    expect_error(
        scpKeepEffect(se, effects = c("foo", "condition")),
        "effects %in% scpModelEffectNames.object, name. are not all TRUE"
    )
    expect_error(
        scpKeepEffect(se, effects = c("foo", "foo2")),
        "effects %in% scpModelEffectNames.object, name. are not all TRUE"
    )
    ## No effect = return residuals
    expect_identical(
        scpRemoveBatchEffect(se),
        scpKeepEffect(se)
    )
    ## Remove one effect
    expect_identical(
        scpRemoveBatchEffect(se, effects = "batch"),
        scpKeepEffect(se, effects = "condition")
    )
    ## Remove more effects
    expect_identical(
        scpRemoveBatchEffect(se, effects = c("batch", "condition")),
        scpKeepEffect(se, effects = NULL),
    )
    ## keep intercept
    expect_identical(
        scpRemoveBatchEffect(se, effects = "batch", intercept = FALSE),
        scpKeepEffect(se, effects = "condition", intercept = TRUE)
    )
    ## batch correct from second model
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch, name = "model2")
    expect_identical(
        scpRemoveBatchEffect(se, effects = "batch", intercept = FALSE,
                             name = "model2"),
        scpKeepEffect(se, effects = "condition", intercept = TRUE,
                      name = "model2")
    )
})

test_that("addReducedDims", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    require("nipals")
    compRes <- scpComponentAnalysis(se)
    ## sce is not a SCE = error
    expect_error(
        addReducedDims(se, compRes),
        "'sce' must be a SingleCellExperiment object. Transform your data using 'as.sce, .SingleCellExperiment..'."
    )
    ## empty list = no effect
    require("SingleCellExperiment")
    sce <- as(se, "SingleCellExperiment")
    expect_identical(
        addReducedDims(sce, List()),
        sce
    )
    ## add one table
    exp <- as.matrix(compRes$bySample$unmodelled[, 1:2])
    attr(exp, "proportionVariance") <- metadata(compRes$bySample$unmodelled)$proportionVariance
    expect_identical(
        reducedDims(addReducedDims(sce, compRes$bySample["unmodelled"])),
        List(unmodelled = exp)
    )
    ## add more table
    exp <- lapply(compRes$bySample, function(x) {
        out <- as.matrix(x[, 1:2])
        attr(out, "proportionVariance") <- metadata(x)$proportionVariance
        out
    })
    expect_identical(
        reducedDims(addReducedDims(sce, compRes$bySample)),
        List(exp)
    )
})

test_that(".getPCs", {
    ## x does not contain columns starting with PC = error
    expect_error(
        .getPCs(data.frame(foo = "bar")),
        "Invalid table.*provided the 'bySample'.*scpComponentAnalysis"
    )
    ## x is not a DataFrame = error
    expect_error(
        .getPCs(data.frame(PC = "bar")),
        "Invalid table.*provided the 'bySample'.*scpComponentAnalysis"
    )
    ## x has no metadata = error
    expect_error(
        .getPCs(DataFrame(PC = "bar")),
        "Invalid table.*provided the 'bySample'.*scpComponentAnalysis"
    )
    ## x has no proportionVariance slot in metadata = error
    x <- DataFrame(PC = "bar")
    metadata(x)$foo <- "bar"
    expect_error(
        .getPCs(x),
        "Invalid table.*provided the 'bySample'.*scpComponentAnalysis"
    )
    ## Works for a table with a single PC
    x <- DataFrame(PC = 1:10)
    metadata(x)$proportionVariance <- c(PC = 1)
    exp <- matrix(1:10, ncol = 1, dimnames = list(NULL, "PC"))
    attr(exp, "proportionVariance") <- c(PC = 1)
    expect_identical(
        .getPCs(x),
        exp
    )
    ## Works for multiple PC columns
    x <- DataFrame(PC1 = 1:10, PC2 = 1:10)
    metadata(x)$proportionVariance <- c(PC1 = 1, PC2 = 2)
    exp <- matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("PC1", "PC2")))
    attr(exp, "proportionVariance") <- c(PC1 = 1, PC2 = 2)
    expect_identical(
        .getPCs(x),
        exp
    )
    ## Works for multiple PC columns with annotations, also testing
    ## case sensitivity
    x <- DataFrame(PC = 1:10, PC2 = 1:10,
                   cell = letters[1:10], pc = letters[1:10])
    metadata(x)$proportionVariance <- c(PC = 1, PC2 = 2)
    exp <- matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("PC", "PC2")))
    attr(exp, "proportionVariance") <- c(PC = 1, PC2 = 2)
    expect_identical(
        .getPCs(x),
        exp
    )
})
