## ----Constructor ----

test_that("ScpModel", {
    x <- ScpModel()
    expect_true(inherits(x, "ScpModel"))
    expect_true(inherits(x@scpModelFormula, "formula"))
    expect_identical(x@scpModelInputIndex, numeric())
    expect_identical(x@scpModelFilterThreshold, numeric())
    expect_identical(x@scpModelFitList, List())
    x2 <- new("ScpModel")
    x2@scpModelFitList <- List()
    expect_identical(x, x2)
})

## ----Test exported getters ----

## Internal function that creates a minimal SE object as expected by
## scplainer for unit testing ScpModel class methods
## @param nr Number of rows
## @param nc Number of columns
.createMinimalData <- function(nr = 10, nc = 5) {
    require("SummarizedExperiment")
    a <- matrix(1, nr, nc)
    rownames(a) <- letters[1:nr]
    colnames(a) <- LETTERS[1:nc]
    se <- SummarizedExperiment(assays = List(assay = a))
    list(se = se, a = a)
}

## Internal function that creates a mock List of ScpModelFit objects
## for unit testing ScpModel class methods
## @param model An ScpModel object
## @param features A character() with the names of the features for
##     which to create mock ScpModelFit objects
.addScpModelFitList <- function(model, features) {
    fitList <- as(lapply(1L:length(features), function(i) {
        ScpModelFit(n = i * i, p = i)
    }), "List")
    names(fitList) <- features
    model@scpModelFitList <- fitList
    model
}

test_that("scpModelFormula", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model = error
    expect_error(
        scpModelFormula(se),
        regexp = "No 'ScpModel'"
    )
    ## When no formula = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFormula(se),
        regexp = "scpModelFormula.*test1.*scpModelWorkflow"
    )
    ## Retrieve formula
    model@scpModelFormula <- ~ var1 + var2
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFormula(se),
        ~ var1 + var2
    )
})

test_that("scpModelInput", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model = error
    expect_error(
        scpModelInput(se),
        regexp = "No 'ScpModel'.*scpModelWorkflow"
    )
    ## No input = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelInput(se),
        regexp = "scpModelFilterThreshold.*test1.*scpModelWorkflow"
    )

    ## Retrieve model input
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model@scpModelInputIndex <- 1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelInput(se, filtered = FALSE), a)
    ## Test the 'filtered' argument
    model <- .addScpModelFitList(model, rownames(se))
    ## Filter = 5 => remove half of the featutres
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelInput(se, filtered = TRUE), a[5:nrow(a), ])
    ## Same but with 1 row (test drop = FALSE)
    model@scpModelFilterThreshold <- nrow(a)
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelInput(se, filtered = TRUE), a[nrow(a), , drop = FALSE])
    ## Test when filtering is disabled
    expect_identical(scpModelInput(se, filtered = FALSE), a)
})

test_that("scpModelFilterThreshold", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## No model = error
    expect_error(
        scpModelFilterThreshold(se),
        regexp = "scpModelFilterThreshold.*test1.*scpModelWorkflow"
    )
    ## No threshold = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFilterThreshold(se),
        regexp = "scpModelFilterThreshold.*test1.*scpModelWorkflow"
    )
    ## Retrieve threshold
    model@scpModelFilterThreshold <- 3
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelFilterThreshold(se), 3)
})

test_that("scpModelFilterNPRatio", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## When no fit object = error
    expect_error(
        scpModelFilterNPRatio(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## Retrieve NP ratio
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    exp <- as.numeric(1:nrow(a))
    names(exp) <- rownames(a)
    ## No filtering (threshold = 0), filtered = FALSE
    model@scpModelFilterThreshold <- 0
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelFilterNPRatio(se, filtered = FALSE), exp)
    ## No filtering (threshold = 0), filtered = TRUE
    expect_identical(scpModelFilterNPRatio(se, filtered = TRUE),  exp)
    ## With filtering, filtered = TRUE
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelFilterNPRatio(se, filtered = TRUE),
                     exp[5:nrow(a)])
    ## With filtering, filtered = FALSE
    expect_identical(scpModelFilterNPRatio(se, filtered = FALSE), exp)
})

test_that("scpModelResiduals", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelResiduals(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No residuals = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelResiduals(se),
        regexp = "Residuals.*test1.*scpModelWorkflow"
    )
    ## Retrieve residuals
    resids <- lapply(seq_len(nrow(se)), function(x) {
        structure(rep(0, ncol(se)), .Names = colnames(se))
    })
    names(resids) <- rownames(se)
    resids <- as(resids, "List")
    model@scpModelFitList <- mendoapply(function(fl, res) {
        names(res) <- colnames(se)
        fl@residuals <- res
        fl
        }, model@scpModelFitList, resids)
    ## No filtering, no joining
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelResiduals(se, join = FALSE, filtered = FALSE),
        resids
    )
    ## No filtering, with joining
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = FALSE),
        do.call(rbind, resids)
    )
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelResiduals(se, join = FALSE, filtered = TRUE),
        resids[5:nrow(se)]
    )
    ## With filtering, with joining
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = TRUE),
        do.call(rbind, resids[5:nrow(se)])
    )
    ## Test drop = FALSE
    model@scpModelFilterThreshold <- 10
    metadata(se)[["test1"]] <- model
    exp <- t(resids[[10]])
    rownames(exp) <- rownames(se)[10]
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = TRUE),
        exp
    )
})

test_that("scpModelEffects", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelEffects(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No effects in model assays = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelEffects(se),
        regexp = "Effect.*test1.*scpModelWorkflow"
    )
    ## Retrieve effects
    effects <- lapply(seq_len(nrow(se)), function(i) {
        out <- lapply(c("Var1", "Var2"), function(j) {
            structure(rep(0, ncol(se)), .Names = colnames(se))
        })
        names(out) <- c("Var1", "Var2")
        as(out, "List")
    })
    names(effects) <- rownames(se)
    effects <- as(effects, "List")
    model@scpModelFitList <- mendoapply(function(fl, eff) {
        fl@effects <- eff
        fl
    }, model@scpModelFitList, effects)
    ## No filtering, no joining
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelEffects(se, join = FALSE, filtered = FALSE),
        effects
    )
    ## No filtering, with joining
    model@scpModelFormula <- ~ 1 + Var1 + Var2
    metadata(se)[["test1"]] <- model
    eff_mat <- matrix(
        0, ncol = ncol(se), nrow = nrow(se), dimnames = dimnames(se)
    )
    expect_identical(
        scpModelEffects(se, join = TRUE, filtered = FALSE),
        List(Var1 = eff_mat, Var2 = eff_mat)
    )
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelEffects(se, join = FALSE, filtered = TRUE),
        effects[5:nrow(se)]
    )
    ## With filtering, with joining
    expect_identical(
        scpModelEffects(se, join = TRUE, filtered = TRUE),
        List(Var1 = eff_mat[5:nrow(se), ], Var2 = eff_mat[5:nrow(se), ])
    )
    ## Test drop = FALSE
    model@scpModelFilterThreshold <- 10
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelEffects(se, join = TRUE, filtered = TRUE),
        List(Var1 = eff_mat[5:nrow(se), , drop = FALSE],
             Var2 = eff_mat[5:nrow(se), , drop = FALSE])
    )
})

test_that("scpModelNames", {
    require(SummarizedExperiment)
    ## SE metadata is empty = error
    se <- SummarizedExperiment()
    expect_error(
        scpModelNames(se),
        regexp = "No 'ScpModel' found in object.*scpModelWorkflow"
    )
    ## SE metadata does not contain an ScpModel = error
    metadata(se)$foo <- "bar"
    expect_error(
        scpModelNames(se),
        regexp = "No 'ScpModel' found in object.*scpModelWorkflow"
    )
    ## SE metadata contains empty model
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["model1"]] <- model
    expect_identical(scpModelNames(se), "model1")
    ## Additional element in SE metadata does not change output
    metadata(se)[["foo"]] <- "bar"
    expect_identical(scpModelNames(se), "model1")
    ## Works with multiple models
    metadata(se)[["model2"]] <- model
    expect_identical(scpModelNames(se), c("model1", "model2"))
})

## ---- Test internal getters ----

test_that("scpModel", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## There is no metadata = error
    expect_error(
        scpModel(se),
        regexp = "No 'ScpModel' found in object"
    )
    ## There is no ScpModel in metadata = error
    metadata(se)$foo <- "bar"
    expect_error(
        scpModel(se),
        regexp = "No 'ScpModel' found in object"
    )
    ## Retrieving more than one model = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    model@scpModelInputIndex <- 1
    metadata(se)[["test2"]] <- model
    expect_error(
        scpModel(se, c("test1", "test2")),
        regexp = "length.*1.*not TRUE"
    )
    ## Retrieving an assay that does not exist = error
    expect_error(
        scpModel(se, "test3"),
        regexp = "test3.*not found.*scpModelWorkflow"
    )
    ## Retrieving an assay without a name default to first assay
    expect_identical(scpModel(se), scpModel(se, "test1"))
    expect_identical(scpModel(se)@scpModelInputIndex, numeric())
    ## Retrieving an assay with a string
    expect_identical(scpModel(se, "test2")@scpModelInputIndex, 1)
})

test_that("scpModelInputIndex", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelEffects(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## When no input index = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelInputIndex(se),
        regexp = "scpModelInputIndex.*test1.*scpModelWorkflow"
    )
    ## Retrieve model input
    model@scpModelInputIndex <- 1
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelInputIndex(se),
        1
    )
})

test_that("scpModelFitList", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model = error
    expect_error(
        scpModelFitList(se),
        regexp = "ScpModel.*object.*scpModelWorkflow"
    )
    ## When no scpModelFitList index = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFitList(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## Retrieve scpModelFitList
    ## No filtering
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    fl <- model@scpModelFitList
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFitList(se, filtered = FALSE),
        fl
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFitList(se, filtered = TRUE),
        fl[(5:nrow(se))]
    )
})

test_that("scpModelFitElement", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model = error
    expect_error(
        scpModelFitElement(se),
        regexp = "ScpModel.*object.*scpModelWorkflow"
    )
    ## When no scpModelFitList  = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFitElement(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## Unknown element = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFitElement(se, what = "foo"),
        regexp = "foo.*not a slot of an ScpModelFit"
    )
    ## Empty element = error
    expect_error(
        scpModelFitElement(se, what = "Residuals"),
        regexp = "Residuals.*ScpModelFit.*model 'test1'[.]"
    )
    ## Test improving error message
    expect_error(
        scpModelFitElement(se, what = "Residuals", helpMessage = "foo!"),
        regexp = "Residuals.*ScpModelFit.*model 'test1'[.] foo!"
    )
    ## Retrieve element
    resids <- lapply(seq_len(nrow(se)), function(x) {
        structure(rep(0, ncol(se)), .Names = colnames(se))
    })
    names(resids) <- rownames(se)
    resids <- as(resids, "List")
    model@scpModelFitList <- mendoapply(function(fl, res) {
        names(res) <- colnames(se)
        fl@residuals <- res
        fl
    }, model@scpModelFitList, resids)
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFitElement(se, what = "Residuals", filtered = FALSE),
        resids
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFitElement(se, what = "Residuals", filtered = TRUE),
        resids[5:nrow(se)]
    )
})

test_that("scpModelN", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelN(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## note the 'n' slot can never be missing, so no error possible
    ## Retrieve N
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    n <- as.integer(seq_len(nrow(se))^2)
    names(n) <- rownames(se)
    ## No filtering
    expect_identical(scpModelN(se, filtered = FALSE), n)
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelN(se, filtered = TRUE), n[5:nrow(se)])
})

test_that("scpModelP", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelP(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## note the 'n' slot can never be missing, so no error possible
    ## Retrieve N
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    p <- as.integer(seq_len(nrow(se)))
    names(p) <- rownames(se)
    ## No filtering
    expect_identical(scpModelP(se, filtered = FALSE), p)
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelP(se, filtered = TRUE), p[5:nrow(se)])
})

test_that("scpModelCoefficients", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelCoefficients(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No coefficients = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelCoefficients(se),
        regexp = "Coefficients.*test1.*scpModelWorkflow"
    )
    ## Retrieve coefficients
    coefs <- lapply(seq_len(nrow(se)), function(x) {
        structure(rep(0, 3), .Names = paste0("param", 1:3))
    })
    names(coefs) <- rownames(se)
    coefs <- as(coefs, "List")
    model@scpModelFitList <- mendoapply(function(fl, coef) {
        fl@coefficients <- coef
        fl
    }, model@scpModelFitList, coefs)
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelCoefficients(se, filtered = FALSE),
        coefs
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelCoefficients(se, filtered = TRUE),
        coefs[5:nrow(se)]
    )
})

test_that("scpModelDf", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelDf(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No df = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelDf(se),
        regexp = "Df.*test1.*scpModelWorkflow"
    )
    ## Retrieve df
    for (i in seq_len(nrow(se))) {
        model@scpModelFitList[[i]]@df <- i
    }
    df <- structure(1:nrow(se), .Names = rownames(se))
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelDf(se, filtered = FALSE), df)
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelDf(se, filtered = TRUE), df[5:nrow(se)])
})

test_that("scpModelVar", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelVar(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No var = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelVar(se),
        regexp = "Var.*test1.*scpModelWorkflow"
    )
    ## Retrieve var
    for (i in seq_len(nrow(se))) {
        model@scpModelFitList[[i]]@var <- i
    }
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelVar(se, filtered = FALSE),
        structure(1:nrow(se), .Names = rownames(se))
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelVar(se, filtered = TRUE),
        structure(1:nrow(se), .Names = rownames(se))[5:nrow(se)]
    )
})

test_that("scpModelUvcov", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelUvcov(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No Uvcov = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelUvcov(se),
        regexp = "Uvcov.*test1.*scpModelWorkflow"
    )
    ## Retrieve Uvcov
    p <- 3
    uvcov <- lapply(seq_len(nrow(se)), function(x) {
        matrix(0, p, p, dimnames = list(paste0("param", 1:p), paste0("param", 1:p)))
    })
    names(uvcov) <- rownames(se)
    uvcov <- as(uvcov, "List")
    model@scpModelFitList <- mendoapply(function(fl, u) {
        fl@uvcov <- u
        fl
    }, model@scpModelFitList, uvcov)
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelUvcov(se, filtered = FALSE),
        uvcov
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelUvcov(se, filtered = TRUE),
        uvcov[5:nrow(se)]
    )
})

test_that("scpModelVcov", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelVcov(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No Uvcov = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelVcov(se),
        regexp = "Vcov.*test1.*scpModelWorkflow"
    )
    ## Retrieve Vcov
    p <- 3
    uvcov <- lapply(seq_len(nrow(se)), function(x) {
        out <- matrix(0, p, p, dimnames = list(paste0("param", 1:p), paste0("param", 1:p)))
        diag(out) <- 1
        out
    })
    names(uvcov) <- rownames(se)
    uvcov <- as(uvcov, "List")
    var <- 2
    model@scpModelFitList <- mendoapply(function(fl, u) {
        fl@uvcov <- u
        fl@var <- var
        fl
    }, model@scpModelFitList, uvcov)
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelVcov(se, filtered = FALSE),
        endoapply(uvcov, function(x) x * var)
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelVcov(se, filtered = TRUE),
        endoapply(uvcov[5:10], function(x) x * var)
    )
})

test_that("scpModelIntercept", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelIntercept(se),
        regexp = "scpModelFitList.*test1.*scpModelWorkflow"
    )
    ## No coefficients = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelIntercept(se),
        regexp = "Coefficients.*test1.*scpModelWorkflow"
    )
    ## Retrieve coefficients
    coefs <- lapply(seq_len(nrow(se)), function(x) {
        structure(c(1, 0, 0), .Names = c("(Intercept)", paste0("param", 2:3)))
    })
    names(coefs) <- rownames(se)
    coefs <- as(coefs, "List")
    model@scpModelFitList <- mendoapply(function(fl, coef) {
        fl@coefficients <- coef
        fl
    }, model@scpModelFitList, coefs)
    ## No filtering
    metadata(se)[["test1"]] <- model
    exp <- structure(rep(1, nrow(se)), .Names = rownames(se))
    expect_identical(
        scpModelIntercept(se, filtered = FALSE),
        exp
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelIntercept(se, filtered = TRUE),
        exp[5:nrow(se)]
    )
})

test_that("scpModelFeatureNames", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelFeatureNames(se),
        regexp = "scpModelFilterThreshold.*test1.*scpModelWorkflow"
    )
    ## Retrieve feature names
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    model@scpModelFilterThreshold <- 0
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFeatureNames(se),
        rownames(se)
    )
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFeatureNames(se),
        rownames(se)[5:10]
    )
})

test_that("scpModelEffectNames", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model = error
    expect_error(
        scpModelEffectNames(se),
        regexp = "No 'ScpModel'"
    )
    ## When no formula = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelEffectNames(se),
        regexp = "scpModelFormula.*test1.*scpModelWorkflow"
    )
    ## Retrieve effect names: one var
    model@scpModelFormula <- ~ var1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), "var1")
    ## Retrieve effect names: multiple var
    model@scpModelFormula <- ~ var1 + var2 + var3
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), c("var1", "var2", "var3"))
    ## Retrieve effect names: no var
    model@scpModelFormula <- ~ 1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), character())
    ## Retrieve effect names: interaction
    model@scpModelFormula <- ~ var1 * var2
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), c("var1", "var2", "var1:var2"))
})

## ---- Test exported setters ----

test_that("scpModelFilterThreshold<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = List(a = matrix(1, 5, 5)))
    model <- ScpModel()
    metadata(se)[["model"]] <- model
    ## Value is not of class numeric = error
    expect_error(
        scpModelFilterThreshold(se) <- "foo",
        regexp = "scpModelFilterThreshold.*numeric.*not TRUE"
    )
    ## Value has length > 1 = error
    expect_error(
        scpModelFilterThreshold(se) <- 1:3,
        regexp = "length.value. == 1 is not TRUE"
    )
    ## Value is < 1 = error
    expect_error(
        scpModelFilterThreshold(se) <- 0,
        regexp = "value >= 1 is not TRUE"
    )
    ## Value is NULL = length 0 = error
    expect_error(
        scpModelFilterThreshold(se) <- NULL,
        regexp = "length.value. == 1 is not TRUE"
    )
    ## Correct case
    scpModelFilterThreshold(se) <- 1
    expect_identical(metadata(se)$model@scpModelFilterThreshold, 1)
})

## ---- Test internal setters ----

test_that("scpModel<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    ## 'value' is not an ScpModel object
    expect_error(
        scpModel(se) <- matrix(),
        regexp = "ScpModel.*is not TRUE"
    )
    ## Add model
    ## name is not provided and there is no model to replace = error
    expect_error(
        scpModel(se) <- model,
        regexp = "No 'ScpModel'"
    )
    ## name is provided and model not initialised = create
    scpModel(se, "test") <- model
    expect_identical(
        metadata(se)[["test"]],
        model
    )
    ## Replace model
    ## name is not provided and model initialised = replace
    model2 <- model
    model2@scpModelInputIndex <- 1
    scpModel(se) <- model
    expect_identical(
        metadata(se)[["test"]],
        model
    )
    ## name is provided and model initialised = replace
    model2@scpModelInputIndex <- 2
    scpModel(se, "test") <- model2
    expect_identical(
        metadata(se)[["test"]],
        model2
    )
    ## Remove model (ie value is NULL)
    ## name is provided and model not initialised = nothing happens
    se <- SummarizedExperiment()
    se2 <- se
    scpModel(se2, "test") <- NULL
    expect_identical(se, se2)
    ## name is provided and model initialised = remove
    scpModel(se2, "test") <- model ## add a model to immediately remove
    scpModel(se2, "test") <- NULL
    expect_identical(se, se2)
    ## name is not provided and model initialised = remove the default model
    se <- SummarizedExperiment()
    scpModel(se, "test1") <- model ## this is the default model
    scpModel(se, "test2") <- model
    expect_identical(
        names(metadata(se)),
        c("test1", "test2")
    )
    scpModel(se) <- NULL
    expect_identical(metadata(se), list(test2 = model))
    ## name is not provided and model not initialised = error
    se <- SummarizedExperiment()
    expect_error(
        scpModel(se) <- NULL,
        regexp = "No 'ScpModel'"
    )
})

test_that("scpModelInputIndex<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(
        assay1 = matrix(1, 2, 2),
        assay2 = matrix(1, 2, 2)
    ))
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    ## Object has no dimension names
    expect_error(
        scpModelInputIndex(se) <- NA,
        regexp = "!is.null.colnames.object.. is not TRUE"
    )
    colnames(se) <- LETTERS[1:2]
    expect_error(
        scpModelInputIndex(se) <- NA,
        regexp = "!is.null.rownames.object.. is not TRUE"
    )
    rownames(se) <- letters[1:2]
    ## Value has wrong type = error
    expect_error(
        scpModelInputIndex(se) <- NA,
        regexp = "must be a character, numeric or logical"
    )
    expect_error(
        scpModelInputIndex(se) <- factor(1),
        regexp = "must be a character, numeric or logical"
    )
        expect_error(
        scpModelInputIndex(se) <- data.frame(assay = 1),
        regexp = "must be a character, numeric or logical"
    )
    expect_error(
        scpModelInputIndex(se) <- matrix(1),
        regexp = "must be a character, numeric or logical"
    )
    ## Value points to multiple assays = error
    expect_error(
        scpModelInputIndex(se) <- c(TRUE, TRUE),
        regexp = "'i' points to multiple input assays."
    )
    expect_error(
        scpModelInputIndex(se) <- 1:2,
        regexp = "'i' points to multiple input assays."
    )
    expect_error(
        scpModelInputIndex(se) <- c("assay1", "assay2"),
        regexp = "'i' points to multiple input assays."
    )
    ## Value points to out of bound index = error
    expect_error(
        scpModelInputIndex(se) <- c(FALSE, FALSE, TRUE),
        regexp = "out of bounds"
    )
    expect_error(
        scpModelInputIndex(se) <- 3,
        regexp = "out of bounds"
    )
    expect_error(
        scpModelInputIndex(se) <- "assay3",
        regexp = "'assay3' not found"
    )
    ## Replace input assay for missing model = error
    expect_error(
        scpModelInputIndex(se, "missingModel") <- 1,
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace input assay for empty model = add
    scpModelInputIndex(se) <- 1
    expect_identical(
        metadata(se)[["test"]]@scpModelInputIndex,
        1
    )
    ## Replace input assay for non-empty model = replace
    scpModelInputIndex(se) <- 2
    expect_identical(
        metadata(se)[["test"]]@scpModelInputIndex,
        2
    )
    ## Make sure it still works when assays in SE are not named
    se <- SummarizedExperiment(assays = list(
        matrix(1, 2, 2),
        matrix(1, 2, 2)
    ))
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    scpModelInputIndex(se) <- 1
    expect_identical(
        metadata(se)[["test"]]@scpModelInputIndex,
        1
    )
})

test_that("scpModelFitList<-", {
    ## TODO
    ## Test error when fitList has wrong length
    ## Test error when fitList has wrong names
})


test_that("scpModelFormula<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = List(a = matrix(1, 5, 5)))
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    ## Value is not of class formula = error
    expect_error(
        scpModelFormula(se) <- matrix(),
        regexp = "formula.*not TRUE"
    )
    ## Variable in model formula are absent from colData = error
    expect_error(
        scpModelFormula(se, "test1") <- ~ 1 + var1,
        regexp = "empty"
    )
    se$var1 <- 1
    se$var2 <- 2
    ## Replace model formula for missing model = error
    expect_error(
        scpModelFormula(se, "missingModel") <- ~ 1 + var1,
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Variable in model is not in colData = error
    expect_error(
        scpModelFormula(se, "test1") <- ~ 1 + var3,
        regexp = "missing.*var3"
    )
    ## The formula has no intercept = error
    expect_error(
        scpModelFormula(se, "test1") <- ~ 0 + var1,
        regexp = "The formula must contain an intercept"
    )
    ## Replace model formula for existing empty model = add
    scpModelFormula(se) <- ~ 1 + var1
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1
    )
    ## Replace model formula for existing non-empty model = replace
    scpModelFormula(se) <- ~ 1 + var1 + var2
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1 + var2
    )
    ## If model contains response variable = remove response + warning
    expect_warning(
        scpModelFormula(se) <- y + x ~ 1 + var1,
        regexp = "Response variables are ignored"
    )
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1
    )
})


test_that("scpModelDesign<-", {
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    se <- SummarizedExperiment(assays = list(
        assay1 = matrix(1, m, n),
        assay2 = matrix(1, m, n)
    ))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    dm <- matrix(1, n, 3, dimnames = list(colnames(se), letters[1:3]))
    ## Value is not of class matrix = error
    expect_error(
        scpModelDesign(se) <- array(),
        regexp = "matrix.*not TRUE"
    )
    ## Replace model design for missing model = error
    expect_error(
        scpModelDesign(se, "missingModel") <- dm,
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace model design with missing dimension names = error
    expect_error(
        scpModelDesign(se, "test") <- matrix(),
        regexp = "'scpModelDesign' is missing sample names."
    )
    ## Replace model design with wrong dimensions (rows) = error
    expect_error(
        scpModelDesign(se, "test") <- dm[1, , drop = FALSE],
        regexp = paste("Expected", n, "but found 1")
    )
    ## Replace model design with wrong dimension names = error
    se2 <- se
    dm2 <- dm
    colnames(se2) <- letters[1:n]
    rownames(dm2) <- rev(colnames(se2))
    expect_error(
        scpModelDesign(se2, "test") <- dm2,
        regexp = "scpModelDesign.*do not match"
    )
    ## Replace model design for existing empty model = add
    scpModelDesign(se) <- dm
    expect_identical(
        metadata(se)[["test"]]@scpModelDesign,
        dm
    )
    ## Replace model design for existing non-empty model = replace
    scpModelDesign(se) <- dm * 2
    expect_identical(
        metadata(se)[["test"]]@scpModelDesign,
        dm * 2
    )
    ## No dimnames = error
    dimnames(dm) <- NULL
    expect_error(
        scpModelDesign(se) <- dm,
        regexp = "'scpModelDesign' is missing sample names."
    )
})

test_that("scpModelCoefficients<-", {
    n <- 5
    m <- 10
    p <- 2
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(
        assay1 = matrix(1, m, n),
        assay2 = matrix(1, m, n)
    ))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    md <- matrix(1, n, p)
    mc <- matrix(1, m, p)
    dimnames(md) <- list(letters[1:n], letters[1:p])
    dimnames(mc) <- list(letters[1:m], letters[1:p])
    model@scpModelDesign <- md
    metadata(se)[["test"]] <- model
    ## Value is not of class matrix = error
    expect_error(
        scpModelCoefficients(se) <- c(),
        regexp = "matrix.*not TRUE"
    )
    ## Replace model coefficients without dimension names = error
    expect_error(
        scpModelCoefficients(se) <- matrix(1, m, p),
        regexp = "'scpModelCoefficients' is missing feature names."
    )
    ## Replace model coefficients with wrong dimensions = error
    expect_error(
        scpModelCoefficients(se) <- mc[1, , drop = FALSE],
        regexp = paste("Expected", m, "but found 1")
    )
    expect_error(
        scpModelCoefficients(se) <- mc[, 1, drop = FALSE],
        regexp = paste("Expected", p, "but found 1")
    )
    ## Replace model coefficients with wrong dimension names = error
    mc2 <- mc
    rownames(mc2) <- rev(rownames(se))
    expect_error(
        scpModelCoefficients(se) <- mc2,
        regexp = "scpModelCoefficients.*do not match"
    )
    mc2 <- mc
    colnames(mc2) <- rev(colnames(md))
    expect_error(
        scpModelCoefficients(se) <- mc2,
        regexp = "scpModelCoefficients.*do not match"
    )
    ## Replace model coefficients for missing model = error
    expect_error(
        scpModelCoefficients(se, "missingModel") <- mc2,
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace model coefficients for empty model = error
    ## It makes no sense to add coefficients if no design is specified
    se2 <- se
    metadata(se2)[["empty"]] <- ScpModel()
    expect_error(
        scpModelCoefficients(se2, "empty") <- mc,
        regexp = "No available 'scpModelDesign'"
    )
    ## Replace model coefficients for existing model = replace
    scpModelCoefficients(se) <- mc + 10
    expect_identical(
        metadata(se)[["test"]]@scpModelCoefficients,
        mc + 10
    )
})

test_that("scpModelAssays<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    ## Value is not of class Assays = error
    expect_error(
        scpModelAssays(se) <- matrix(),
        regexp = "Assays.*not TRUE"
    )
    ## Replace model assays for missing model = error
    expect_error(
        scpModelAssays(se, "missingModel") <- Assays(),
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace model assays with wrong dimensions = error
    n <- 5
    m <- 10
    se <- SummarizedExperiment(assays = List(a1 = matrix(1, m, n)))
    dimnames(se) <- list(letters[1:m], LETTERS[1:n])
    metadata(se)[["test"]] <- model
    m2 <- matrix(1, m, n)
    ## Replace model assays without dimension names = error
    expect_error(
        scpModelAssays(se, "test") <- Assays(List(b1 = m2)),
        regexp = "'scpModelAssays' is missing feature names."
    )
    ## Replace model assays with wrong dimensions = error
    m2 <- matrix(1, m, n, dimnames = dimnames(se))
    m2 <- rbind(m2, z = 1)
    expect_error(
        scpModelAssays(se, "test") <- Assays(List(b1 = m2)),
        regexp = paste("Expected", m, "but found", m + 1)
    )
    m2 <- matrix(1, m, n, dimnames = dimnames(se))
    m2 <- cbind(m2, z = 1)
    expect_error(
        scpModelAssays(se, "test") <- Assays(List(b1 = m2)),
        regexp = paste("Expected", n, "but found", n + 1)
    )
    ## Replace model assays with wrong dimension names = error
    m2 <- matrix(1, m, n, dimnames = list(rev(rownames(se)), colnames(se)))
    expect_error(
        scpModelAssays(se, "test") <- Assays(List(b1 = m2)),
        regexp = "Names in 'scpModelAssays' do not match with the rownames in 'object'"
    )
    m2 <- matrix(1, m, n, dimnames = list(rownames(se), rev(colnames(se))))
    expect_error(
        scpModelAssays(se, "test") <- Assays(List(b1 = m2)),
        regexp = "Names in 'scpModelAssays' do not match with the colnames in 'object'"
    )
    ## Replace model assays for empty model = add
    a <- matrix(1, m, n)
    dimnames(a) <- dimnames(se)
    assays <- Assays(List(a1 = a, a2 = a))
    metadata(se)[["test"]] <- ScpModel()
    scpModelAssays(se, "test") <- assays
    expect_identical(
        metadata(se)[["test"]]@scpModelAssays,
        assays
    )
    ## Replace model assays for existing model = replace,
    ## no dimnames = message
    a <- matrix(1, m, n)
    assays <- Assays(List(a1 = a))
    expect_error(
        scpModelAssays(se) <- assays,
        regexp = "'scpModelAssays' is missing feature names."
    )
    ## Replace model assays for existing model = replace
    dimnames(a) <- dimnames(se)
    assays <- Assays(List(a1 = a))
    scpModelAssays(se) <- assays
    expect_identical(
        metadata(se)[["test"]]@scpModelAssays,
        assays
    )
})

test_that("scpModelResiduals<-", {
    ## TODO make sure the residuals must be named with the samlpe names
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    se <- SummarizedExperiment(assays = List(a1 = matrix(1, m, n)))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    ## Value is not of class matrix = error
    expect_error(
        scpModelResiduals(se) <- c(),
        regexp = "matrix.*not TRUE"
    )
    ## Replace model residuals for missing model = error
    expect_error(
        scpModelResiduals(se, "missingModel") <- matrix(),
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace model residuals for model without prior requirement
    ## No initialization = error
    expect_error(
        scpModelResiduals(se, "test") <- matrix(2, m + 1, n),
        regexp = "scpModelInputIndex.*scpModelPrepare"
    )
    metadata(se)[["test"]]@scpModelInputIndex <- 1
    ## No model estimation = error
    expect_error(
        scpModelResiduals(se, "test") <- matrix(2, m + 1, n),
        regexp = "scpModelCoefficients.*scpModelRun"
    )
    metadata(se)[["test"]]@scpModelCoefficients <- matrix(1, n, 10)
    ## no dimnames = error
    mr <- matrix(2, m, n)
    expect_error(
        scpModelResiduals(se) <- mr,
        regexp = "'scpModelAssays' is missing feature names."
    )
    rownames(mr) <- rownames(se)
    expect_error(
        scpModelResiduals(se) <- mr,
        regexp = "'scpModelAssays' is missing sample names."
    )
    ## Replace model residuals with wrong dimensions = error
    mr <- matrix(2, m + 1, n)
    dimnames(mr) <- list(letters[1:(m + 1)], letters[1:n])
    expect_error(
        scpModelResiduals(se, "test") <- mr,
        regexp = paste("Expected", m, "but found", m + 1)
    )
    mr <- matrix(2, m, n + 1)
    dimnames(mr) <- list(letters[1:m], letters[1:(n + 1)])
    expect_error(
        scpModelResiduals(se, "test") <- mr,
        regexp = paste("Expected", n, "but found", n + 1)
    )
    ## Replace model assays with wrong dimension names = error
    mr <- matrix(1, m, n)
    dimnames(mr) <- dimnames(se)
    colnames(mr) <- rev(colnames(mr))
    expect_error(
        scpModelResiduals(se, "test") <- mr,
        regexp = "'scpModelAssays' do not match"
    )
    dimnames(mr) <- dimnames(se)
    rownames(mr) <- rev(rownames(mr))
    expect_error(
        scpModelResiduals(se, "test") <- mr,
        regexp = "'scpModelAssays' do not match"
    )
    ## Replace model residuals for existing model = replace
    dimnames(mr) <- dimnames(se)
    scpModelResiduals(se) <- mr * 10
    expect_identical(
        getListElement(metadata(se)[["test"]]@scpModelAssays, "residuals"),
        mr * 10
    )
})

test_that("scpModelEffects<-", {
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    a <- matrix(1, m, n)
    se <- SummarizedExperiment(assays = List(a1 = a))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    ## Value is not of class Assays = error
    expect_error(
        scpModelEffects(se) <- matrix(),
        regexp = "Assays.*not TRUE"
    )
    ## Replace model effects for missing model = error
    expect_error(
        scpModelEffects(se, "missingModel") <- Assays(),
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace model effects for model without prior requirement
    ## No initialization = error
    expect_error(
        scpModelEffects(se, "test") <- Assays(),
        regexp = "scpModelInputIndex.*scpModelPrepare"
    )
    ## No model estimation = error
    metadata(se)[["test"]]@scpModelInputIndex <- 1
    expect_error(
        scpModelEffects(se, "test") <- Assays(),
        regexp = "scpModelCoefficients.*scpModelRun"
    )
    ## No formula = error
    metadata(se)[["test"]]@scpModelCoefficients <- matrix(1, n, 10)
    expect_error(
        scpModelEffects(se, "test") <- Assays(),
        regexp = "scpModelFormula.*scpModelPrepare"
    )
    ## Some effects are missing = error
    metadata(se)[["test"]]@scpModelFormula <- ~ 1 + var1 + var2
    var1 <- matrix(2, m, n)
    dimnames(var1) <- dimnames(se)
    expect_error(
        scpModelEffects(se) <- Assays(List(var1 = var1)),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## Some effects are not modelled = error
    expect_error(
        scpModelEffects(se) <- Assays(List(
            var1 = var1,
            var2 = var1,
            var3 = var1
        )),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## Replace model effects without dimnames = error
    dimnames(var1) <- NULL
    expect_error(
        scpModelEffects(se) <- Assays(List(var1 = var1, var2 = var1)),
        regexp = "'scpModelAssays' is missing feature names."
    )
    rownames(var1) <- rownames(se)
    expect_error(
        scpModelEffects(se) <- Assays(List(var1 = var1, var2 = var1)),
        regexp = "'scpModelAssays' is missing sample names."
    )
    ## Replace model effects with wrong dimensions = error
    var1 <- matrix(2, m, n, dimnames = dimnames(se))
    var1 <- rbind(var1, z = 1)
    expect_error(
        scpModelEffects(se) <- Assays(List(var1 = var1, var2 = var1)),
        regexp = paste0(
            "'scpModelAssays' does not contain the expected number ",
            "of features. Expected ", m, " but found ", m + 1, "."
        )
    )
    var1 <- matrix(2, m, n, dimnames = dimnames(se))
    var1 <- cbind(var1, z = 1)
    expect_error(
        scpModelEffects(se) <- Assays(List(var1 = var1, var2 = var1)),
        regexp = paste0(
            "'scpModelAssays' does not contain the expected number ",
            "of samples. Expected ", n, " but found ", n + 1, "."
        )
    )
    ## Replace model effects for empty model = add
    var1 <- matrix(2, m, n, dimnames = dimnames(se))
    scpModelEffects(se) <- Assays(List(var1 = var1, var2 = var1 * 10))
    expect_identical(
        getListElement(metadata(se)[["test"]]@scpModelAssays, "var1"),
        var1
    )
    expect_identical(
        getListElement(metadata(se)[["test"]]@scpModelAssays, "var2"),
        var1 * 10
    )
    ## Replace model effects for existing model = replace
    scpModelEffects(se) <- Assays(List(
        var1 = var1 * 20,
        var2 = var1 * 30
    ))
    expect_identical(
        getListElement(metadata(se)[["test"]]@scpModelAssays, "var1"),
        var1 * 20
    )
    expect_identical(
        getListElement(metadata(se)[["test"]]@scpModelAssays, "var2"),
        var1 * 30
    )
})

test_that("scpModelCovariance <-", {
    skip("todo")
})

test_that("scpModelSumsOfSquares<-", {
    skip("deprecated")
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    a <- matrix(1, m, n)
    se <- SummarizedExperiment(assays = List(a1 = a))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    model@scpModelFormula <- ~ 1 + var1 + var2
    effects <- c("var1", "var2")
    p <- length(effects)
    metadata(se)[["test"]] <- model
    ## Value is not a numeric vector = error
    expect_error(
        scpModelSumsOfSquares(se) <- c(1),
        regexp = "matrix.*not TRUE"
    )
    expect_error(
        scpModelSumsOfSquares(se) <- matrix("foo"),
        regexp = "numeric.*not TRUE"
    )
    expect_error(
        scpModelSumsOfSquares(se) <- matrix(TRUE),
        regexp = "numeric.*not TRUE"
    )
    ## Replace model effects for missing model = error
    expect_error(
        scpModelSumsOfSquares(se, "missingModel") <- matrix(1),
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## No coefficient estimation = error
    expect_error(
        scpModelSumsOfSquares(se, "test") <- matrix(1),
        regexp = "scpModelCoefficients.*scpModelRun"
    )
    metadata(se)[["test"]]@scpModelCoefficients <- matrix(1, m, p)
    ## No residuals = error
    expect_error(
        scpModelSumsOfSquares(se, "test") <- matrix(1),
        regexp = "scpModelResiduals.*scpModelRun"
    )
    metadata(se)[["test"]]@scpModelAssays <-
        setListElement(
            metadata(se)[["test"]]@scpModelAssays,
            "residuals", a
        )
    ## SS has wrong dimensions = error
    ssnames <- c("SSresiduals", "SStotal", effects)
    expect_error(
        scpModelSumsOfSquares(se) <- matrix(
            1, nrow = m-1, ncol = p + 2,
            dimnames = list(rownames(se)[-1], ssnames)
        ),
        regexp = paste0(
            "does not contain the expected number of features. ",
            "Expected ", m, " but found ", m - 1
        )
    )
    expect_error(
        scpModelSumsOfSquares(se) <- matrix(
            1, nrow = m, ncol = p + 1,
            dimnames = list(rownames(se), ssnames[-3])
        ),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## SS has wrong dimension names = error
    ## no names
    ss <- ssnonames <- matrix(1, nrow = m, ncol = p + 2)
    dimnames(ss) <- list(rownames(se), ssnames)
    expect_error(
        scpModelSumsOfSquares(se) <- ssnonames,
        regexp = "'scpModelSumsOfSquares' is missing names"
    )
    ## wrong rownames
    expect_error(
        scpModelSumsOfSquares(se) <- ss[rev(seq_len(nrow(ss))), ],
        regexp = "Names in 'scpModelSumsOfSquares' do not match"
    )
    ## Missing SS for one effect
    expect_error(
        scpModelSumsOfSquares(se) <- ss[, -3],
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## Missing SS for "TotalSS" or "ResidualSS"
    ss2 <- ss
    colnames(ss2)[1] <- "SSres"
    expect_error(
        scpModelSumsOfSquares(se) <- ss2,
        regexp = "SSresiduals"
    )
    ss2 <- ss
    colnames(ss2)[2] <- "SStot"
    expect_error(
        scpModelSumsOfSquares(se) <- ss2,
        regexp = "SStotal"
    )
    ## Replace model effects for empty model = add
    scpModelSumsOfSquares(se) <- ss
    expect_identical(
        metadata(se)[["test"]]@scpModelSumsOfSquares,
        ss
    )
    ## Replace model effects for existing model = replace
    scpModelSumsOfSquares(se) <- ss * 10
    expect_identical(
        metadata(se)[["test"]]@scpModelSumsOfSquares,
        ss * 10
    )
})

test_that("scpModelComponent<-", {
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    k <- 3
    se <- SummarizedExperiment(assays = List(a1 = matrix(1, m, n)))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    scores <- matrix(1, n, k, dimnames = list(letters[1:n], LETTERS[1:k]))
    eigenvectors <- matrix(1, m, k, dimnames = list(letters[1:m], LETTERS[1:k]))
    ca <- ScpModelComponent(scores, eigenvectors, 1:k, 1:k)
    ## Value is not of class List = error
    expect_error(
        scpModelComponent(se) <- matrix(),
        regexp = "inherits.value,.*ScpModelComponent.*is not TRUE"
    )
    ## Attempt to replace more than one component = error
    scores <- matrix(1, n, 5, dimnames = list(colnames(se), 1:5))
    eigenvectors <- matrix(1, m, 5, dimnames = list(rownames(se), 1:5))
    smc <- ScpModelComponent(scores, eigenvectors, 1:5, 1:5)
    expect_error(
        scpModelComponent(se, effect = 1:2) <- smc,
        regexp = "length.effect.*1.*not TRUE"
    )
    ## Replace component for missing model = error
    expect_error(
        scpModelComponent(se, "missingModel", effect = "residuals") <- smc,
        regexp = "missingModel.*not found.*scpModelPrepare"
    )
    ## Replace component for unknown effect = error
    metadata(se)[["test"]]@scpModelFormula <- ~ 1 + var1 + var2
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var3") <- smc,
        regexp = "unknown component.*ASCA_var3.*runComponentAnalysis"
    )
    ## Replace component for unknown method = error
    expect_error(
        scpModelComponent(se, method = "AXCA", effect = "var1") <- smc,
        regexp = "unknown component.*AXCA_var1.*runComponentAnalysis"
    )
    ## Replace component with wrong dimensions = error
    smcwrong <- ScpModelComponent(scores, scores, 1:5, 1:5)
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var1") <- smcwrong,
        regexp = paste0(
            "eigenvectors.*not contain.*features.*", m,
            " but found ", n
        )
    )
    smcwrong <- ScpModelComponent(eigenvectors, eigenvectors, 1:5, 1:5)
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var1") <- smcwrong,
        regexp = paste0(
            "scores.*not contain.*samples.*", n,
            " but found ", m
        )
    )
    ## Replace component with wrong dimension names = error
    smcwrong <- ScpModelComponent(
        scores[rev(seq_len(nrow(scores))), ],
        eigenvectors, 1:5, 1:5
    )
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var1") <- smcwrong,
        regexp = "Names in 'scores' do not match with the colnames in 'object'."
    )
    smcwrong <- ScpModelComponent(
        scores,
        eigenvectors[rev(seq_len(nrow(eigenvectors))), ], 1:5, 1:5
    )
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var1") <- smcwrong,
        regexp = "Names in 'eigenvectors' do not match with the rownames in 'object'."
    )
    ## Replace component for known effect and model but absent = add
    scpModelComponent(se, method = "ASCA", effect = "var1") <- smc
    expect_identical(
        metadata(se)[["test"]]@scpModelComponents,
        List(ASCA_var1 = smc)
    )
    ## Replace a second component for known effect and model = append
    scpModelComponent(se, method = "ASCA.E", effect = "var1") <- smc
    expect_identical(
        metadata(se)[["test"]]@scpModelComponents,
        List(ASCA_var1 = smc, ASCA.E_var1 = smc)
    )
    ## Same but for residuals (whatever the method) = add
    scpModelComponent(se, method = "ASCA", effect = "residuals") <- smc
    expect_identical(
        metadata(se)[["test"]]@scpModelComponents,
        List(ASCA_var1 = smc, ASCA.E_var1 = smc, residuals = smc)
    )
    scpModelComponent(se, effect = "residuals") <- smc ## missing method
    expect_identical(
        metadata(se)[["test"]]@scpModelComponents,
        List(ASCA_var1 = smc, ASCA.E_var1 = smc, residuals = smc)
    )
    ## Replace component for existing component = replace
    smc2 <- ScpModelComponent(scores*2, eigenvectors * 2, 1:5, 1:5)
    scpModelComponent(se, method = "ASCA.E", effect = "var1") <- smc2
    expect_identical(
        metadata(se)[["test"]]@scpModelComponents,
        List(ASCA_var1 = smc, ASCA.E_var1 = smc2, residuals = smc)
    )
})

## ---- Test internal utility functions ----

test_that(".defaultModelName", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## Error when no model
    expect_error(
        .defaultModelName(se),
        regexp = "No 'ScpModel' found in object. Use 'scpModelPrepare"
    )
    ## When 1 model, return first name
    metadata(se)[["test"]] <- ScpModel()
    expect_identical(.defaultModelName(se), "test")
    ## When 2 (or more) models, return first name
    metadata(se)[["test2"]] <- ScpModel()
    expect_identical(.defaultModelName(se), "test")
})

test_that(".checkn", {
    n <- 10
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(matrix(1, 10, n)))
    colnames(se) <- letters[1:n]
    ## nnames is NULL = error
    expect_error(
        .checkn(se, NULL, "test"),
        "'test' is missing sample names."
    )
    ## nnames has wrong length = error
    expect_error(
        .checkn(se, letters[1:(n + 1)], "test"),
        paste0("'test' does not contain the expected number of ",
               "samples. Expected ", n, " but found ", n + 1, ".")
    )
    ## nnames is not correct = error
    expect_error(
        .checkn(se, rev(letters[1:n]), "test"),
        "Names in 'test' do not match with the colnames in 'object'"
    )
    ## No error
    expect_identical(
        .checkn(se, letters[1:n], "test"),
        NULL
    )
})

test_that(".checkm", {
    m <- 10
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(matrix(1, m, 20)))
    rownames(se) <- letters[1:m]
    ## mnames is NULL = error
    expect_error(
        .checkm(se, "model", NULL, "test"),
        "'test' is missing feature names."
    )
    ## mnames has wrong length = error
    metadata(se)[["model"]] <- ScpModel()
    expect_error(
        .checkm(se, "model", letters[1:(m + 1)], "test"),
        paste0("'test' does not contain the expected number of ",
               "features. Expected ", m, " but found ", m + 1, ".")
    )
    ## mnames is not correct (without filter) = error
    expect_error(
        .checkm(se, "model", rev(letters[1:m]), "test"),
        "Names in 'test' do not match with the rownames in 'object'"
    )
    ## No error (without filter)
    expect_identical(
        .checkm(se, "model", letters[1:n], "test"),
        NULL
    )
    ## mnames is not correct (with filter) = error
    rowData(se)[["scpModelFilter_model"]] <-
        rep(c(TRUE, FALSE), each = m / 2)
    expect_error(
        .checkm(se, "model", rev(letters[1:(m / 2)]), "test"),
        "Names in 'test' do not match with the rownames in 'object'"
    )
    ## No error (with filer)
    expect_identical(
        .checkm(se, "model", letters[1:(m / 2)], "test"),
        NULL
    )
})

test_that(".checkp", {
    p <- 3
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(matrix(1, 10, 20)))
    colnames(se) <- letters[1:20]
    ## pnames is NULL = error
    expect_error(
        .checkp(se, "model", NULL, "test"),
        "'test' is missing parameter names."
    )
    ## pnames has wrong length = error
    model <- ScpModel()
    model@scpModelDesign <- matrix(1, 20, p)
    dimnames(model@scpModelDesign) <- list(colnames(se), letters[1:p])
    metadata(se)[["model"]] <- model
    expect_error(
        .checkp(se, "model", letters[1:(p + 1)], "test"),
        paste0("'test' does not contain the expected number of ",
               "parameters. Expected ", p, " but found ", p + 1, ".")
    )
    ## pnames is not correct = error
    expect_error(
        .checkp(se, "model", rev(letters[1:p]), "test"),
        "Names in 'test' do not match with the parameter names in 'object'"
    )
    ## No error
    expect_identical(
        .checkp(se, "model", letters[1:p], "test"),
        NULL
    )
})

test_that(".checkModelElement", {
    model <- ScpModel()
    ## NULL
    expect_error(
        .checkModelElement(NULL, "model1", "element1", "More info."),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelFormula = "formula",
    expect_error(
        .checkModelElement(
            model@scpModelFormula,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelInputIndex = "numeric",
    expect_error(
        .checkModelElement(
            model@scpModelInputIndex,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelDesign = "matrix",
    expect_error(
        .checkModelElement(
            model@scpModelDesign,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelCoefficients = "matrix",
    expect_error(
        .checkModelElement(
            model@scpModelCoefficients,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelAssays = "Assays",
    expect_error(
        .checkModelElement(
            model@scpModelAssays,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelSumsOfSquares = "numeric",
    expect_error(
        .checkModelElement(
            model@scpModelSumsOfSquares,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelComponentAnalysis = "SimpleList"
    expect_error(
        .checkModelElement(
            model@scpModelComponents,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
})

test_that(".checkScpModelFormula", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = List(a = matrix(1, 5, 5)))
    ## The formula has no variables = warning
    expect_warning(
        checkres <- .checkScpModelFormula(~ NULL, se),
        regexp = "You provided a formula with no variable to model. "
    )
    expect_identical(checkres, ~ NULL)
    expect_warning(
        checkres <- .checkScpModelFormula(~ 1, se),
        regexp = "You provided a formula with no variable to model. "
    )
    expect_identical(checkres, ~ 1)
    ## The colData is empty = error
    expect_error(
        .checkScpModelFormula(~ 1 + var1, se),
        regexp = "colData\\(object\\) is empty."
    )
    ## The colData is missing 1 variable = error
    se$var1 <- 1
    se$var2 <- 2
    expect_error(
        .checkScpModelFormula(~ 1 + var1 + var2 + var3, se),
        regexp = "missing one or more variables.*var3"
    )
    ## The colData is missing 2 variables = error
    expect_error(
        .checkScpModelFormula(~ 1 + var1 + var2 + var3 + var4, se),
        regexp = "missing one or more variables.*var3, var4"
    )
    ## The formula has no intercept = error
    expect_error(
        .checkScpModelFormula(~ 0 + var1, se),
        regexp = "The formula must contain an intercept"
    )
    ## The formula has a response variable = removed
    expect_warning(
        checkres <- .checkScpModelFormula(y ~ 1 + var1, se),
        regexp = "contains a response"
    )
    expect_warning(
        checkres <- .checkScpModelFormula(y + x ~ 1 + var1 + var2, se),
        regexp = "contains a response"
    )
    expect_identical(checkres, ~ 1 + var1 + var2)
    ## same but without intercept
    expect_warning(
        checkres <- .checkScpModelFormula(y ~ var1, se),
        regexp = "contains a response"
    )
    expect_identical(checkres, ~ var1)
    ## 1 variable is present in colData
    expect_identical(.checkScpModelFormula(~ 1 + var1, se), ~ 1 + var1)
    ## 2 variables are present in colData
    expect_identical(
        .checkScpModelFormula(~ 1 + var1 + var2, se),
        ~ 1 + var1 + var2
    )
    ## dot assignment works
    expect_identical(
        .checkScpModelFormula(~ 1 + ., se),
        ~ 1 + .
    )
})

test_that(".checkExplanatoryVariables", {
    ## No model variables = warning
    expect_warning(
        checkres <- .checkExplanatoryVariables(character(), "var1"),
        regexp = "You provided a formula with no variable to model. "
    )
    expect_identical(checkres, NULL)
    ## The colData is empty = error
    expect_error(
        .checkExplanatoryVariables("var1", character()),
        regexp = "colData\\(object\\) is empty."
    )
    ## The colData is missing 1 variable = error
    expect_error(
        .checkExplanatoryVariables(
            c("var1", "var2", "var3"),
            c("var1", "var2")
        ),
        regexp = "missing one or more variables.*var3"
    )
    ## The colData is missing 2 variables = error
    expect_error(
        .checkExplanatoryVariables(
            c("var1", "var2", "var3", "var4"),
            c("var1", "var2")
        ),
        regexp = "missing one or more variables.*var3, var4.$"
    )
    ## 1 variable is present in colData
    expect_identical(
        .checkExplanatoryVariables("var1", c("var1", "var2")),
        NULL
    )
    ## 2 variables are present in colData
    expect_identical(
        .checkExplanatoryVariables(
            c("var1", "var2"),
            c("var1", "var2")
        ),
        NULL
    )
    ## dot assignment works
    expect_identical(
        .checkExplanatoryVariables(".", c("var1", "var2")),
        NULL
    )
})

test_that(".checkEffects", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(matrix(1, 10, 20)))
    model <- ScpModel()
    model@scpModelFormula <- ~ 1 + var1 + var2
    metadata(se)[["model1"]] <- model
    ## Tested name is empty = error
    expect_error(
        .checkEffects(se, "model1", NULL),
        regexp = "must be named"
    )
    ## All tested names are unmodelled = error
    expect_error(
        .checkEffects(se, "model1", c("var3", "var4")),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## Some tested names are unmodelled = error
    expect_error(
        .checkEffects(se, "model1", c("var2", "var3")),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## Some tested names are missing = error
    expect_error(
        .checkEffects(se, "model1", "var2"),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## All tested names are modelled and present = ok
    expect_identical(
        .checkEffects(se, "model1", c("var1", "var2")),
        NULL
    )
    ## If residuals are not allowed = error
    expect_error(
        .checkEffects(se, "model1", c("var1", "var2", "residuals")),
        regexp = "Missing or unmodelled effects.*: var1, var2[.]"
    )
    ## If residuals are allowed = ok
    expect_identical(
        .checkEffects(
            se, "model1", c("var1", "var2", "residuals"), TRUE
        ),
        NULL
    )
    ## Test that when residuals are allowed, but some effects are
    ## missing, "residuals" appears in the error message
    expect_error(
        .checkEffects(
            se, "model1", c("var1", "residuals"), TRUE
        ),
        regexp = "Missing or unmodelled effects.*: var1, var2, residuals[.]"
    )
})

test_that(".checkComponentName", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## component name is "residuals" = always ok
    expect_identical(
        .checkComponentName(se, tested = "residuals"),
        NULL
    )
    ## Unknown effect = error
    model <- ScpModel()
    model@scpModelFormula <- ~ 1 + var1 + var2
    metadata(se)[["model1"]] <- model
    expect_error(
        .checkComponentName(se, tested = "ASCA_var3"),
        regexp = "unknown component analysis results: ASCA_var3"
    )
    ## Unknown method = error
    expect_error(
        .checkComponentName(se, tested = "AXCA_var1"),
        regexp = "unknown component analysis results: AXCA_var1"
    )
    ## Known method and effect = ok
    for (m in scpModelComponentMethods) {
        for (j in scpModelEffectNames(se)) {
            expect_identical(
                .checkComponentName(se, tested = paste0(m, "_", j)),
                NULL
            )
        }
    }
})

test_that(".removeResponseVariables", {
    ## There is a response and no explanatory variable
    expect_warning(
        test <- .removeResponseVariables(
            formula = y ~ 1,
            data = data.frame()),
        regexp = "The formula contains a response variable"
    )
    expect_identical(test, ~ 1)
    ## There is a response and explanatory variables
    expect_warning(
        test <- .removeResponseVariables(
            formula = y ~ 1 + var1 + var2,
            data = data.frame()
        ),
        regexp = "The formula contains a response variable"
    )
    expect_identical(test, ~ 1 + var1 + var2)
    ## There are multiple responses and explanatory variables
    expect_warning(
        test <- .removeResponseVariables(
            formula = y + x ~ 1 + var1 + var2,
            data = data.frame()
        ),
        regexp = "The formula contains a response variable"
    )
    expect_identical(test, ~ 1 + var1 + var2)
    ## There are multiple responses and explanatory is dot
    expect_warning(
        test <- .removeResponseVariables(
            formula = y + x ~ 1 + .,
            data = data.frame(var1 = NA, var2 = NA)
        ),
        regexp = "The formula contains a response variable"
    )
    expect_identical(test, ~ 1 + .)
    ## No response = no effect
    expect_identical(
        .removeResponseVariables(
            formula = ~ 1 + var1 + var2,
            data = data.frame()
        ),
        ~ 1 + var1 + var2
    )
    ## No response and no explanatory variable = no effect
    expect_identical(
        .removeResponseVariables(
            formula = ~ NULL,
            data = data.frame()
        ),
        ~ NULL
    )
})

test_that(".replaceDotVariable", {
    availableVariables <- paste0("var", 1:3)
    scpModelVariableNames <- "var1"
    ## No dot variable = no effect
    expect_identical(
        .replaceDotVariable(scpModelVariableNames, availableVariables),
        scpModelVariableNames
    )
    ## scpModelVariableNames is empty = no effect
    expect_identical(
        .replaceDotVariable(NULL, availableVariables),
        NULL
    )
    expect_identical(
        .replaceDotVariable(character(), availableVariables),
        character()
    )
    ## availabelVariables is empty = remove dot
    scpModelVariableNames <- c(".", "var1")
    expect_identical(
        .replaceDotVariable(scpModelVariableNames, character()),
        "var1"
    )
    ## modelVariable is only dot = becomes availableVariables
    expect_identical(
        .replaceDotVariable(".", availableVariables),
        availableVariables
    )
    ## modelVariable contains dot + variables (with overlap) = becomes
    ## availableVariables
    expect_identical(
        .replaceDotVariable(c(scpModelVariableNames, "."), availableVariables),
        availableVariables
    )
    ## modelVariable contains dot and partially overlaps = union
    scpModelVariableNames <- paste0("var", 2:5)
    expect_identical(
        .replaceDotVariable(c(scpModelVariableNames, "."), availableVariables),
        union(scpModelVariableNames, availableVariables)
    )
    ## modelVariable contains dot and no overlap = union
    scpModelVariableNames <- paste0("var", 4:6)
    expect_identical(
        .replaceDotVariable(c(scpModelVariableNames, "."), availableVariables),
        union(scpModelVariableNames, availableVariables)
    )
})

test_that(".indexToNumeric", {
    require(SummarizedExperiment)
    ## Value has wrong type = error
    expect_error(
        .indexToNumeric(NA, c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .indexToNumeric(factor(1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .indexToNumeric(data.frame(assay = 1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .indexToNumeric(matrix(1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    ## Value points to multiple assays = error
    expect_error(
        .indexToNumeric(c(TRUE, TRUE), c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    expect_error(
        .indexToNumeric(1:2, c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    expect_error(
        .indexToNumeric(c("test1", "test2"), c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    ## Value points to out of bound index = error
    expect_error(
        .indexToNumeric(c(FALSE, FALSE, TRUE), c("test1", "test2"), "mock"),
        regexp = "out of bounds"
    )
    expect_error(
        .indexToNumeric(3, c("test1", "test2"), "mock"),
        regexp = "out of bounds"
    )
    expect_error(
        .indexToNumeric("test3", c("test1", "test2"), "mock"),
        regexp = "'test3' not found."
    )
    ## Valid numeric = no conversion
    expect_identical(
        .indexToNumeric(1, c("test1", "test2"), "mock"),
        1
    )
    ## Convert from character
    expect_identical(
        .indexToNumeric("test2", c("test1", "test2"), "mock"),
        2L
    )
    ## Convert from logical
    expect_identical(
        .indexToNumeric(c(FALSE, TRUE), c("test1", "test2"), "mock"),
        2L
    )
})

test_that(".filterFromRowData", {
    skip("todo")
})

test_that(".generateFilterName", {
    skip("todo")
})

test_that(".generateComponentNames", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    model@scpModelFormula <- ~ 1 + var1 + var2
    ## No components = error
    metadata(se)[["model1"]] <- model
    expect_error(
        .generateComponentNames(se),
        regexp = "unavailable components.*none.*runComponentAnalysis()"
    )
    ## Select unknown effects = error
    compNames <- c(
        "ASCA_var1", "ASCA_var2", "APCA_var1",
        "APCA_var2", "residuals"
    )
    l <- List(1, 1, 1, 1, 1)
    names(l) <- compNames
    model@scpModelComponents <- l
    metadata(se)[["model1"]] <- model
    expect_error(
        .generateComponentNames(se, effect = "var3"),
        regexp = paste0(
            "unavailable components.*ASCA for var1.*ASCA ",
            "for var2.*APCA for var1.*APCA for var2.*residuals.*runComponentAnalysis()"
        )
    )
    ## Select unknown methods = error
    expect_error(
        .generateComponentNames(se, method = "AXCA"),
        regexp = paste0(
            "unavailable components.*ASCA for var1.*ASCA ",
            "for var2.*APCA for var1.*APCA for var2.*residuals.*runComponentAnalysis()"
        )
    )
    ## Select unknown combination of method and effect = error
    expect_error(
        .generateComponentNames(se, method = "AXCA", effect = "var1"),
        regexp = "method %in% scpModelComponentMethods. is not TRUE"
    )
    expect_error(
        .generateComponentNames(se, method = "ASCA", effect = "var3"),
        regexp = "effect %in% scpModelEffectNames.* is not TRUE"
    )
    ## Not specifying method or effect = select everything
    ## with residuals
    expect_identical(
        .generateComponentNames(se, residuals = TRUE),
        compNames
    )
    ## without residuals
    expect_identical(
        .generateComponentNames(se, residuals = FALSE),
        compNames[-5]
    )
    ## Specifying only method = select all effects
    ## with residuals
    expect_identical(
        .generateComponentNames(se, method = "ASCA", residuals = TRUE),
        c("ASCA_var1", "ASCA_var2", "residuals")
    )
    ## without residuals
    expect_identical(
        .generateComponentNames(se, method = "ASCA", residuals = FALSE),
        c("ASCA_var1", "ASCA_var2")
    )
    ## Specifying only effect = select all methods
    ## with residuals
    expect_identical(
        .generateComponentNames(se, effect = "var1", residuals = TRUE),
        c("ASCA_var1", "APCA_var1", "residuals")
    )
    ## without residuals
    expect_identical(
        .generateComponentNames(se, effect = "var1", residuals = FALSE),
        c("ASCA_var1", "APCA_var1")
    )
    ## Specifying both effect and method
    ## with residuals
    expect_identical(
        .generateComponentNames(se,
                                effect = "var1", method = "ASCA", residuals = TRUE
        ),
        c("ASCA_var1", "residuals")
    )
    ## without residuals
    expect_identical(
        .generateComponentNames(se,
                                effect = "var1", method = "ASCA", residuals = FALSE
        ),
        "ASCA_var1"
    )
})













####---- DUMP ----#####

test_that("scpModelEffectNames", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model formula = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelEffectNames(se),
        regexp = "scpModelFormula.*test1.*scpModelPrepare"
    )
    ## Formula with no variable
    model@scpModelFormula <- ~ 1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), character())
    ## Formula with 1 variable
    model@scpModelFormula <- ~ 0 + var1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), "var1")
    ## Formula with interecept + 1 variable
    model@scpModelFormula <- ~ 1 + var1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), "var1")
    ## Formula with interecept + 2 variables
    model@scpModelFormula <- ~ 1 + var1 + var2
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), c("var1", "var2"))
    ## Formula with interecept + interaction
    model@scpModelFormula <- ~ 1 + var1 : var2
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelEffectNames(se), c("var1:var2"))
    ## Formula with interecept + interaction + main effects
    model@scpModelFormula <- ~ 1 + var1 * var2
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelEffectNames(se),
        c("var1", "var2", "var1:var2")
    )
    ## Formula with .
    se <- SummarizedExperiment(assays = List(a = matrix(1, 5, 5)))
    metadata(se)[["test1"]] <- model
    se$var1 <- 1
    se$var2 <- 2
    se$var3 <- 3
    model@scpModelFormula <- ~ 1 + . + var1:var2
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelEffectNames(se),
        c("var1", "var2", "var3", "var1:var2")
    )
})

test_that("scpModelDesign", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model design = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelDesign(se),
        regexp = "scpModelDesign.*test1.*scpModelPrepare"
    )
    ## Retrieve model deisgn
    model@scpModelDesign <- matrix(1, 5, 5)
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelDesign(se),
        matrix(1, 5, 5)
    )
})

test_that("scpModelCoefficients", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model coefficients = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelCoefficients(se),
        regexp = "scpModelCoefficients.*test1.*scpModelRun"
    )
    ## Retrieve model coefficients
    model@scpModelCoefficients <- matrix(1, 5, 5)
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelCoefficients(se),
        matrix(1, 5, 5)
    )
})


test_that("scpModelSumsOfSquares", {
    skip("deprecated")
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no SS = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelSumsOfSquares(se),
        regexp = "scpModelSumsOfSquares.*test1.*computeSumsOfSquares"
    )
    ## Retrieve sums of squares
    ss <- matrix(1:10)
    model@scpModelSumsOfSquares <- ss
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelSumsOfSquares(se),
        ss
    )
})

test_that("scpModelNames", {
    require(SummarizedExperiment)
    ## SE metadata is empty = error
    se <- SummarizedExperiment()
    expect_error(
        scpModelNames(se),
        regexp = "No 'ScpModel' found in object."
    )
    ## SE metadata does not contain an ScpModel = error
    metadata(se)$foo <- "bar"
    expect_error(
        scpModelNames(se),
        regexp = "No 'ScpModel' found in object."
    )
    ## SE metadata contains empty model
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["model1"]] <- model
    expect_identical(scpModelNames(se), "model1")
    ## Additional element in SE metadata does not change output
    metadata(se)[["foo"]] <- "bar"
    expect_identical(scpModelNames(se), "model1")
    ## Works with multiple models
    metadata(se)[["model2"]] <- model
    expect_identical(scpModelNames(se), c("model1", "model2"))
})

test_that("scpModelComponent", {
    require(SummarizedExperiment)
    n <- 5
    m <- 10
    k <- 3
    se <- SummarizedExperiment(assays = List(a1 = matrix(1, m, n)))
    dimnames(se) <- list(letters[1:m], letters[1:n])
    rd <- DataFrame(rowd1 = letters[1:m], rowd2 = LETTERS[1:m])
    cd <- DataFrame(
        cold1 = letters[1:n], cold2 = LETTERS[1:n],
        row.names = colnames(se)
    )
    colData(se) <- cd
    rowData(se) <- rd
    model <- ScpModel()
    pcnames <- paste0("PC", 1:k)
    scores <- matrix(1, n, k, dimnames = list(colnames(se), pcnames))
    eigenvectors <- matrix(1, m, k, dimnames = list(rownames(se), pcnames))
    percentVar <- 1:k
    smc <- ScpModelComponent(scores, eigenvectors, 1:k, percentVar)
    smcs <- List(
        residuals = smc,
        ASCA_var1 = smc,
        ASCA_var2 = smc,
        APCA_var1 = smc,
        APCA_var2 = smc
    )
    model@scpModelComponents <- smcs
    model@scpModelFormula <- ~ var1 + var2
    metadata(se)[["test"]] <- model
    ## The scpModelComponent slot is empty = error
    metadata(se)[["empty"]] <- ScpModel()
    expect_error(
        scpModelComponent(se, "empty"),
        regexp = "unavailable components.*none.*runComponentAnalysis()"
    )
    ## Select unknown effects = error
    expect_error(
        scpModelComponent(se, effect = "var3"),
        regexp = paste0(
            "unavailable components.*residuals.*ASCA for var1.*ASCA ",
            "for var2.*APCA for var1.*APCA for var2.*runComponentAnalysis()"
        )
    )
    ## Select unknown methods = error
    expect_error(
        scpModelComponent(se, method = "AXCA"),
        regexp = paste0(
            "unavailable components.*residuals.*ASCA for var1.*ASCA ",
            "for var2.*APCA for var1.*APCA for var2.*runComponentAnalysis()"
        )
    )
    ## Select unknown combination of method and effect = error
    expect_error(
        scpModelComponent(se, method = "AXCA", effect = "var1"),
        regexp = "method %in% scpModelComponentMethods. is not TRUE"
    )
    expect_error(
        scpModelComponent(se, method = "ASCA", effect = "var3"),
        regexp = "effect %in% scpModelEffectNames.* is not TRUE"
    )
    ## Not specifying method or effect selects everything
    ## In sample-space
    target <- DataFrame(scores)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target,
            APCA_var1 = target,
            APCA_var2 = target,
            residuals = target
        )
    )
    ## In feature-space
    target <- DataFrame(eigenvectors)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, bySample = FALSE),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target,
            APCA_var1 = target,
            APCA_var2 = target,
            residuals = target
        )
    )
    ## Without residuals
    expect_identical(
        scpModelComponent(se, bySample = FALSE, residuals = FALSE),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target,
            APCA_var1 = target,
            APCA_var2 = target
        )
    )
    ## Specifying only method selects all effects
    ## In sample-space
    target <- DataFrame(scores)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, method = "ASCA"),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target,
            residuals = target
        )
    )
    ## In feature-space
    target <- DataFrame(eigenvectors)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, method = "ASCA", bySample = FALSE),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target,
            residuals = target
        )
    )
    ## Without residuals
    expect_identical(
        scpModelComponent(se,
                          method = "ASCA", bySample = FALSE,
                          residuals = FALSE
        ),
        List(
            ASCA_var1 = target,
            ASCA_var2 = target
        )
    )
    ## Specifying only effect selects all methods
    ## In sample-space
    target <- DataFrame(scores)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, effect = "var1"),
        List(
            ASCA_var1 = target,
            APCA_var1 = target,
            residuals = target
        )
    )
    ## In feature-space
    target <- DataFrame(eigenvectors)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, effect = "var1", bySample = FALSE),
        List(
            ASCA_var1 = target,
            APCA_var1 = target,
            residuals = target
        )
    )
    ## Without residuals
    expect_identical(
        scpModelComponent(se,
                          effect = "var1", bySample = FALSE,
                          residuals = FALSE
        ),
        List(
            ASCA_var1 = target,
            APCA_var1 = target
        )
    )
    ## Specifying both effect and method
    ## In sample-space
    target <- DataFrame(scores)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, effect = "var1", method = "ASCA"),
        List(
            ASCA_var1 = target,
            residuals = target
        )
    )
    ## In feature-space
    target <- DataFrame(eigenvectors)
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se,
                          effect = "var1", method = "ASCA",
                          bySample = FALSE
        ),
        List(
            ASCA_var1 = target,
            residuals = target
        )
    )
    ## Without residuals
    expect_identical(
        scpModelComponent(se,
                          effect = "var1", method = "ASCA",
                          bySample = FALSE, residuals = FALSE
        ),
        List(ASCA_var1 = target)
    )
    ## Testing adding annotations
    ## In sample-space
    target <- cbind(DataFrame(scores), colData(se))
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se, effect = "var1", withAnnotations = TRUE),
        List(
            ASCA_var1 = target,
            APCA_var1 = target,
            residuals = target
        )
    )
    ## In feature-space
    target <- cbind(DataFrame(eigenvectors), rowData(se))
    metadata(target)$proportionVariance <- percentVar
    expect_identical(
        scpModelComponent(se,
                          effect = "var1", bySample = FALSE,
                          withAnnotations = TRUE
        ),
        List(
            ASCA_var1 = target,
            APCA_var1 = target,
            residuals = target
        )
    )
})

test_that("parnames", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## When no model design = error
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    expect_error(
        parnames(se),
        regexp = "scpModelDesign.*test1.*initialiseScpModel"
    )
    ## Unnamed columns design
    m <- matrix(1, 10, 2)
    model@scpModelDesign <- m
    metadata(se)[["test1"]] <- model
    expect_identical(parnames(se), NULL)
    ## 2 variables design
    colnames(m) <- c("var1", "var2")
    model@scpModelDesign <- m
    metadata(se)[["test1"]] <- model
    expect_identical(parnames(se), c("var1", "var2"))
})
