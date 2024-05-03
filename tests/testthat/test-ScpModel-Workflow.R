## ---- ScpModel workflow ----

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

test_that("scpModelWorkflow", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 2

    #### Test standard case
    ## Run workflow
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    ## workflow adds a ScpModel object
    expect_true(inherits(metadata(se)$model, "ScpModel"))
    ## ScpModel contains expected formula
    expect_identical(
        scpModel(se)@scpModelFormula,
        ~ 1 + condition
    )
    ## ScpModel contains expected input index
    expect_identical(
        scpModel(se)@scpModelInputIndex,
        1
    )
    ## ScpModel contains expected filter threshold
    expect_identical(
        scpModel(se)@scpModelFilterThreshold,
        1
    )
    ## ScpModel contains expected List of ScpModelFit objects
    expect_true(inherits(scpModel(se)@scpModelFitList, "List"))
    expect_true(all(sapply(scpModel(se)@scpModelFitList, inherits, "ScpModelFit")))
    for (i in 1:nrow(se)) {
        ## ScpModelFit contains expected coefficients
        expect_equal(
            scpModelCoefficients(se)[[i]],
            c("(Intercept)" = 1.5, condition1 = -0.5),
            tolerance = 1E-3
        )
        ## ScpModelFit contains expected residuals
        expect_equal(
            scpModelResiduals(se, join = FALSE)[[i]],
            structure(rep(0, ncol(se)), .Names = colnames(se)),
            tolerance = 1E-3
        )
        ## ScpModelFit contains expected effects
        expect_equal(
            scpModelEffects(se, join = FALSE)[[i]],
            List(condition = structure(
                as.numeric(se$condition) - 1.5,
                .Names = colnames(se)
            )),
            tolerance = 1E-3
        )
        ## ScpModelFit contains expected df
        expect_identical(scpModelDf(se)[[i]], ncol(se) - 2)
        ## ScpModelFit contains expected var
        expect_true(scpModelVar(se)[[i]] < 1E-6)
        ## ScpModelFit contains expected unscaled covariance
        expect_equal(
            scpModelUvcov(se)[[i]],
            matrix(
                c(0.2082, -0.0416, -0.0416, 0.2082), 2, 2,
                dimnames = list(
                    c("(Intercept)", "condition1"),
                    c("(Intercept)", "condition1")
                )
            ),
            tolerance = 1E-4
        )
        ## ScpModelFit contains expected variable levels
        expect_identical(
            scpModelFitLevels(scpModelFitList(se)[[i]]),
            List(condition = as.character(1:2))
        )
    }

    #### Test presence of missing values
    ## Randomly inject missing values
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    ## Model already exists = warning
    expect_warning(
        se <- scpModelWorkflow(se, formula = ~ 1 + condition),
        "An element called 'model' is already present in the metadata. The associated content will be overwritten."
    )
    ## Some features are no longer estimable
    estimable <- c("a", "c", "e", "h", "i")
    ## ScpModelFit contains expected residuals (joined)
    res <- matrix(0, nrow(se), ncol(se), dimnames = dimnames(se))
    res[is.na(assay(se))] <- NA
    expect_equal(
        scpModelResiduals(se),
        res[estimable, ],
        tolerance = 1E-3
    )
    for (i in estimable) {
        ## ScpModelFit contains expected coefficients. The residuals
        ## should be equal (within a small tolerance) as above
        expect_equal(
            scpModelCoefficients(se)[[i]],
            c("(Intercept)" = 1.5, condition1 = -0.5),
            tolerance = 1E-3
        )
        ## ScpModelFit contains expected coefficients. The effect
        ## should be equal (within a small tolerance) as above, subset
        ## for only observed values.
        effects <- structure(as.numeric(se$condition) - 1.5, .Names = colnames(se))
        effects <- effects[!is.na(assay(se)[i, ])]
        expect_equal(
            scpModelEffects(se, join = FALSE)[[i]],
            List(condition = effects),
            tolerance = 1E-3
        )
        ## ScpModelFit contains expected degrees of freedom.
        expect_identical(scpModelDf(se)[[i]], sum(!is.na(assay(se)[i, ])) - 2)
        ## ScpModelFit contains expected variance, similarly to above
        ## but with a higher tolerance since less data points.
        ## Note we don't test Uvcov since it possibly changes for each
        ## feature
        expect_true(scpModelVar(se)[[i]] < 1E-5)
        ## ScpModelFit contains expected levels. The levels
        ## should be equal as above
        expect_identical(
            scpModelFitLevels(scpModelFitList(se)[[i]]),
            List(condition = as.character(1:2))
        )
    }
    ## Features that are not estimated contain an empty ScpModelFit
    for (i in rownames(se)[!rownames(se) %in% estimable]) {
        expect_identical(
            scpModelFitList(se)[[i]],
            ScpModelFit(as.integer(sum(!is.na(assay(se)[i, ]))), 0L)
        )
    }

    ## metadata element (not an ScpModel) already exist = warning
    metadata(se)$foo <- "bar"
    expect_warning(
        se <- scpModelWorkflow(se, formula = ~ 1 + condition, name = "foo"),
        "An element called 'foo' is already present in the metadata. The associated content will be overwritten."
    )
    expect_identical(
        scpModel(se, "foo"),
        scpModel(se, "model")
    )
    ## Test i as numeric, and name as character
    se <- scpModelWorkflow(
        se, formula = ~ 1 + condition, name = "model2", i = 1
    )
    expect_identical(
        scpModel(se, "model"),
        scpModel(se, "model2")
    )
    ## Test verbose
    expect_silent(tmp <- scpModelWorkflow(
        se, formula = ~ 1 + condition, name = "model3", i = 1, verbose = FALSE
    ))
})


test_that(".fitScpModel", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    exp <- scpModel(se)@scpModelFitList

    ## Test standard case
    expect_identical(
        exp,
        .fitScpModel(se, "model", verbose = TRUE)
    )
    ## Test verbose
    expect_silent(.fitScpModel(se, "model", verbose = FALSE))
})

test_that(".checkAnnotations", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 2
    model <- ScpModel()
    model@scpModelFormula <- ~ 1 + condition
    metadata(se)$model <- model
    ## colData has no rownames = error
    colnames(se) <- NULL
    expect_error(
        .checkAnnotations(se),
        "colData.object.' must have row names."
    )
    ## colData has missing values = error
    colnames(se) <- LETTERS[1:5]
    se$condition <- NA
    expect_error(
        .checkAnnotations(se),
        "Sample annotations .colData. cannot contain missing values."
    )
    ## singular designs = error
    se$condition <- 1
    expect_error(
        .checkAnnotations(se),
        "The design matrix is .near. singular."
    )
    ##
})