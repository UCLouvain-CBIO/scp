## ---- Constructor ----

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

## ---- Test exported getters ----

## Internal function that creates a minimal SE object as expected by
## scplainer for unit testing ScpModel class methods
## @param nr Number of rows
.createMinimalData <- function(nr = 10) {
    require("SummarizedExperiment")
    nc <- nr * nr  # Set the number of columns to be sufficient for the largest row
    a <- matrix(NA, nr, nc)
    for (i in 1:nr) {
        a[i, 1:(i * i)] <- 1:(i * i)
    }
    rownames(a) <- paste0("row", 1:nr)
    colnames(a) <- paste0("col", 1:nc)
    se <- SummarizedExperiment(assays = List(assay = a))
    list(se = se, a = a)
}
## Internal function that creates a mock List of ScpModelFit objects
## for unit testing ScpModel class methods
## @param model An ScpModel object
## @param features A character() with the names of the features for
##     which to create mock ScpModelFit objects
## @param dfP A logical indicating whether to create coefficients for
## the purpose of having a p for test
## @param coefR A logical indicating whether to create coefficients
.addScpModelFitList <- function(model, features, dfP = FALSE, coefR = FALSE) {
    fitList <- as(lapply(1L:length(features), function(i) {
        smf <- ScpModelFit()
        if (dfP) smf@df <- i*(i-1)
        if (coefR) {
            smf@coefficients <- c(1, 0)
            names(smf@coefficients) <- c("(Intercept)", "conditionB")
        }
        smf
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
        regexp = "No available scpModelInputIndex for model 'test1'.*scpModelWorkflow"
    )
    ## Retrieve model input
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model@scpModelInputIndex <- 1
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelInput(se, filtered = FALSE), a)
    ## Test the 'filtered' argument
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    ## dfP = TRUE to be able to retrieve p to compute NP ratio
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
        regexp = "scpModelInputIndex.*test1.*scpModelWorkflow"
    )
    ## Retrieve NP ratio
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    exp <- as.numeric(1:nrow(a))
    names(exp) <- rownames(a)
    ## No filtering (threshold = 0), filtered = FALSE
    model@scpModelFilterThreshold <- 0
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_identical(scpModelFilterNPRatio(se, filtered = FALSE), exp)
    ## No filtering (threshold = 0), filtered = TRUE
    expect_identical(scpModelFilterNPRatio(se, filtered = TRUE),  exp)
    ## With filtering, filtered = TRUE
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
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
    ## No residuals = error
    l <- .createMinimalData(); se <- l$se; a <- l$a
    colData(se)$condition <- factor(rep(c("A", "B"), each = nrow(se)/2))
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE, coefR = TRUE)
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_error(
        scpModelResiduals(se, filtered = FALSE),
        regexp = "scpModelFormula.*test1.*scpModelWorkflow"
    )
    ## Retrieve residuals
    resids <- assay(se) - 1
    resids <- lapply(1:nrow(resids), function(i) {
        setNames(resids[i, ], colnames(resids))
    })
    for (i in seq_along(resids)) {
        if (sum(is.finite(resids[[i]])) <= 1) {
            resids[[i]][is.finite(resids[[i]])] <- NA
        }
    }
    resids <- as(resids, "List")
    ## No filtering, no joining
    metadata(se)[["test1"]] <- model
    scpModelFormula(se) <- ~ 1 + condition
    scpModelInputIndex(se) <- 1
    expect_identical(
        scpModelResiduals(se, join = FALSE, filtered = FALSE),
        resids
    )
    # No filtering, with joining
    joined_resids <- BiocGenerics::do.call(rbind, resids)
    rownames(joined_resids) <- rownames(se)
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = FALSE),
        joined_resids
    )
    # With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    scpModelFormula(se) <- ~ 1 + condition
    expect_identical(
        scpModelResiduals(se, join = FALSE, filtered = TRUE),
        resids[5:nrow(se)]
    )
    ## With filtering, with joining
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = TRUE),
        joined_resids[5:nrow(se),]
    )
    ## Test drop = FALSE
    model@scpModelFilterThreshold <- 10
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    scpModelFormula(se) <- ~ 1 + condition
    exp <- t(resids[[10]])
    rownames(exp) <- rownames(se)[10]
    expect_identical(
        scpModelResiduals(se, join = TRUE, filtered = TRUE),
        exp
    )
})

test_that("scpModelEffects", {
    require(SummarizedExperiment)
    model <- ScpModel()
    l <- .createMinimalData(); se <- l$se; a <- l$a
    colData(se)$condition <- factor(rep(c("A", "B"), each = nrow(se)/2))
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE, coefR = TRUE)
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_error(
        scpModelEffects(se, filtered = FALSE),
        regexp = "Formula.*test1.*scpModelWorkflow"
    )
    scpModelFormula(se) <- ~ 1 + condition
    ## Retrieve effects
    effects <- lapply(seq_len(nrow(se)), function(i) {
        out <- list( "condition" = structure(c(rep(0, i*i), rep(NA, 100 - i*i)),
                                             .Names = colnames(se)))
        out$condition <- out$condition[!is.na(out$condition)]
        as(out, "List")
    })
    names(effects) <- rownames(se)
    effects <- as(effects, "List")
    effects[[1]] <- numeric(0)
    effects[[2]] <- numeric(0)
    # No filtering, no joining
    # expect_identical(
    #     scpModelEffects(se, join = FALSE, filtered = FALSE),
    #     effects
    # )
    ## No filtering, with joining
    # model@scpModelFormula <- ~ 1 + Var1 + Var2
    # metadata(se)[["test1"]] <- model
    # eff_mat <- matrix(
    #     0, ncol = ncol(se), nrow = nrow(se), dimnames = dimnames(se)
    # )
    # expect_identical(
    #     scpModelEffects(se, join = TRUE, filtered = FALSE),
    #     List(Var1 = eff_mat, Var2 = eff_mat)
    # )
    # ## With filtering, no joining
    # model@scpModelFilterThreshold <- 5
    # metadata(se)[["test1"]] <- model
    # scpModelInputIndex(se) <- 1
    # expect_identical(
    #     scpModelEffects(se, join = FALSE, filtered = TRUE),
    #     effects[5:nrow(se)]
    # )
    # ## With filtering, with joining
    # expect_identical(
    #     scpModelEffects(se, join = TRUE, filtered = TRUE),
    #     List(Var1 = eff_mat[5:nrow(se), ], Var2 = eff_mat[5:nrow(se), ])
    # )
    # ## Test drop = FALSE
    # model@scpModelFilterThreshold <- 10
    # metadata(se)[["test1"]] <- model
    # scpModelInputIndex(se) <- 1
    # expect_identical(
    #     scpModelEffects(se, join = TRUE, filtered = TRUE),
    #     List(Var1 = eff_mat[10, , drop = FALSE],
    #          Var2 = eff_mat[10, , drop = FALSE])
    # )
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
        scpModelCoefficients(se),
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    fl <- model@scpModelFitList
    metadata(se)[["test1"]] <- model
    expect_identical(
        scpModelFitList(se, filtered = FALSE),
        fl
    )
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelFitElement(se, what = "foo"),
        regexp = "foo.*not a slot of an ScpModelFit"
    )
})

test_that("scpModelN", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model

    expect_error(
        scpModelN(se),
        regexp = "n.*test1.*scpModelWorkflow"
    )
    ## Retrieve N
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    n <- seq_len(nrow(se))^2
    names(n) <- rownames(se)
    ## No filtering
    expect_identical(scpModelN(se, filtered = FALSE), n)
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_identical(scpModelN(se, filtered = TRUE), n[5:nrow(se)])
})

test_that("scpModelP", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    model <- ScpModel()
    metadata(se)[["test1"]] <- model
    ## no model = error
    expect_error(
        scpModelP(se, filtered = FALSE),
        regexp = "scpModelInputIndex.*test1.*scpModelWorkflow"
    )
    ## Retrieve N
    l <- .createMinimalData(); se <- l$se; a <- l$a
    model <- .addScpModelFitList(model, rownames(se))
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_error(
        scpModelP(se, filtered = FALSE),
        regexp = "Df.*test1.*scpModelWorkflow"
    )
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    p <- as.numeric(seq_len(nrow(se)))
    names(p) <- rownames(se)
    ## No filtering
    expect_identical(scpModelP(se, filtered = FALSE), p)
    ## With filtering, no joining
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelCoefficients(se),
        regexp = "Coefficients.*test1.*scpModelWorkflow"
    )
    ## Retrieve coefficients
    coefs <- lapply(seq_len(nrow(se)), function(x) {
        structure(rep(0, x), .Names = paste0("param", 1:x))
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
    scpModelInputIndex(se) <- 1
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
        model@scpModelFitList[[i]]@df <- i * (i-1)
    }
    df <- structure(sapply(
        seq_len(nrow(se)), function(i) i * (i - 1)
        ), .Names = rownames(se))
    ## No filtering
    metadata(se)[["test1"]] <- model
    expect_identical(scpModelDf(se, filtered = FALSE), df)
    ## With filtering
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
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
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
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
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
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
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    metadata(se)[["test1"]] <- model
    expect_error(
        scpModelIntercept(se),
        regexp = "Coefficients.*test1.*scpModelWorkflow"
    )
    ## Retrieve coefficients
    coefs <- lapply(seq_len(nrow(se)), function(x) {
        names <- if(x > 1) {
            c("(Intercept)", paste0("param", seq_len(x - 1)))
            } else c("(Intercept)")
        structure(
            c(1, rep(0, x - 1)),
            .Names = names)
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
    scpModelInputIndex(se) <- 1
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
    model <- .addScpModelFitList(model, rownames(se), dfP = TRUE)
    model@scpModelFilterThreshold <- 0
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
    expect_identical(
        scpModelFeatureNames(se),
        rownames(se)
    )
    model@scpModelFilterThreshold <- 5
    metadata(se)[["test1"]] <- model
    scpModelInputIndex(se) <- 1
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
        scpModelFormula(se, "test") <- ~ 1 + var1,
        regexp = "empty"
    )
    se$var1 <- 1
    se$var2 <- 2
    ## Replace model formula for missing model = error
    expect_error(
        scpModelFormula(se, "missingModel") <- ~ 1 + var1,
        regexp = "missingModel.*not found.*scpModelWorkflow"
    )
    ## Variable in model is not in colData = error
    expect_error(
        scpModelFormula(se, "test") <- ~ 1 + var3,
        regexp = "missing.*var3"
    )
    ## The formula has no intercept = warning
    expect_warning(
        scpModelFormula(se, "test") <- ~ 0 + var1,
        regexp = "No intercept in the formula. It is added automatically."
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
        regexp = "The formula contains a response variable and is ignored"
    )
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1
    )
    ## dot works
    scpModelFormula(se) <-  ~ 1 + .
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1 + var2
    )
    scpModelFormula(se) <-  ~ 1 + var1 + .
    expect_identical(
        metadata(se)[["test"]]@scpModelFormula,
        ~ 1 + var1 + var2
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
    ## Object has no dimension names = error
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
        regexp = "missingModel.*not found.*scpModelWorkflow"
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
    assays(se) <- unname(assays(se))
    scpModelInputIndex(se) <- 1
    expect_identical(
        metadata(se)[["test"]]@scpModelInputIndex,
        1
    )
})

test_that("scpModelFitList<-", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = list(
        assay1 = matrix(1, 2, 2),
        assay2 = matrix(1, 2, 2)
    ))
    dimnames(se) <- list(LETTERS[1:2], letters[1:2])
    model <- ScpModel()
    metadata(se)[["test"]] <- model
    smFit <- ScpModelFit()
    ## Value has wrong type = is not a List = error
    expect_error(
        scpModelFitList(se) <- smFit,
        "inherits.value, .List.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- list(smFit), ## base list does not work
        "inherits.value, .List.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- "smFit",
        "inherits.value, .List.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- 1L,
        "inherits.value, .List.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- matrix(),
        "inherits.value, .List.. is not TRUE"
    )
    ## Value is a List with elements that are not ScpModelFit
    ## objects = error
    expect_error(
        scpModelFitList(se) <- List(A = "foo", B = 1L),
        "all.*value.*inherits.*ScpModelFit.* not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- List(A = smFit, B = 1L),
        "all.*value.*inherits.*ScpModelFit.* not TRUE"
    )
    ## Value has wrong length = error
    expect_error(
        scpModelFitList(se) <- List(),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- List(A = smFit),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- List(A = smFit, B = smFit, C = smFit),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    ## Value has wrong names = error
    expect_error(
        scpModelFitList(se) <- List(A = smFit, C = smFit),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    expect_error(
        scpModelFitList(se) <- List(smFit, smFit),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    ## Value has correct names but wrong order = error
    expect_error(
        scpModelFitList(se) <- List(B = smFit, A = smFit),
        "identical.rownames.object., names.value.. is not TRUE"
    )
    ## Replace scpModelFitList for empty model = add
    scpModelFitList(se) <- List(A = smFit, B = smFit)
    expect_identical(
        metadata(se)[["test"]]@scpModelFitList,
        List(A = smFit, B = smFit)
    )
    ## Replace input assay for non-empty model = replace
    smFit2 <- ScpModelFit()
    scpModelFitList(se) <- List(A = smFit, B = smFit2)
    expect_identical(
        metadata(se)[["test"]]@scpModelFitList,
        List(A = smFit, B = smFit2)
    )
})

## ---- Test internal utility functions ----

test_that(".defaultModelName", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## Error when no model
    expect_error(
        .defaultModelName(se),
        regexp = "No 'ScpModel' found in object. Use 'scpModelWorkflow"
    )
    ## When 1 model, return first name
    metadata(se)[["test"]] <- ScpModel()
    expect_identical(.defaultModelName(se), "test")
    ## When 2 (or more) models, return first name
    metadata(se)[["test2"]] <- ScpModel()
    expect_identical(.defaultModelName(se), "test")
})

test_that(".checkModelName", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## No model in object = error
    expect_error(
        .checkModelName(se, "foo"),
        "No 'ScpModel' found in object. Use 'scpModelWorkflow"
    )
    metadata(se)$foo <- "bar"
    metadata(se)$model1 <- ScpModel()
    metadata(se)$model2 <- ScpModel()
    metadata(se)$foo2 <- "bar"
    ## Name points to multiple models = error
    expect_error(
        .checkModelName(se, c("model1", "model2")),
        "length.name.*1 is not TRUE"
    )
    ## Name not in metadata = error
    expect_error(
        .checkModelName(se, "bar"),
        "Model name 'bar' not found in object.*scpModelWorkflow"
    )
    ## Name points to element that is not an ScpModel = error
    expect_error(
        .checkModelName(se, "foo"),
        "Model name 'foo' not found in object.*scpModelWorkflow"
    )
    ## Name points to a valid element
    expect_identical(
        .checkModelName(se, "model1"),
        "model1"
    )
    expect_identical(
        .checkModelName(se, "model2"),
        "model2"
    )
    ## No name = default model
    expect_identical(
        .checkModelName(se, "model1"),
        "model1"
    )
})

test_that(".checkModelElement", {
    model <- ScpModel()
    ## NULL
    expect_error(
        .checkModelElement(NULL, "model1", "element1", "More info."),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelFormula = error
    expect_error(
        .checkModelElement(
            model@scpModelFormula,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelInputIndex = error
    expect_error(
        .checkModelElement(
            model@scpModelInputIndex,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelFilterThreshold = error
    expect_error(
        .checkModelElement(
            model@scpModelFilterThreshold,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## Empty scpModelFitList = error
    expect_error(
        .checkModelElement(
            model@scpModelFitList,
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## x is an empty array = error
    expect_error(
        .checkModelElement(
            c(),
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## x is an empty List = error
    expect_error(
        .checkModelElement(
            List(),
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## x is a List or list with empty elements
    expect_error(
        .checkModelElement(
            List(a = c(), b = c()),
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    expect_error(
        .checkModelElement(
            list(a = c(), b = c()),
            "model1", "element1", "More info."
        ),
        regexp = "element1.*model1.*More info."
    )
    ## x is an non empty vector = NULL
    expect_identical(
        .checkModelElement(
            1:10,
            "model1", "element1", "More info."
        ),
        NULL
    )
    ## x is an non empty List or list = NULL
    expect_identical(
        .checkModelElement(
            List(a = 1:10, b = 1:10),
            "model1", "element1", "More info."
        ),
        NULL
    )
    expect_identical(
        .checkModelElement(
            list(a = 1:10, b = 1:10),
            "model1", "element1", "More info."
        ),
        NULL
    )
})

test_that(".checkScpModelFormula", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment(assays = List(a = matrix(1, 5, 5)))
    ## The formula has no variables = error
    expect_error(
        .checkScpModelFormula(~ NULL, se),
        regexp = "You provided a formula with no variable to model."
    )
    expect_error(
        .checkScpModelFormula(~ 1, se),
        regexp = "You provided a formula with no variable to model."
    )
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
    ## The formula has no intercept = warning
    expect_warning(
        .checkScpModelFormula(~ 0 + var1, se),
        regexp = "No intercept in the formula. It is added automatically."
    )
    ## The formula has a response variable = removed
    expect_warning(
        test <- .checkScpModelFormula(y ~ 1 + var1, se),
        regexp = "contains a response"
    )
    expect_warning(
        test <- .checkScpModelFormula(y + x ~ 1 + var1 + var2, se),
        regexp = "contains a response"
    )
    expect_identical(test, ~ 1 + var1 + var2)
    ## same but without intercept
    expect_warning(
        test <- .checkScpModelFormula(y ~ var1, se),
        regexp = "contains a response"
    )
    expect_identical(test, ~ 1 + var1)
    ## 1 variable is present in colData
    expect_identical(
        .checkScpModelFormula(~ 1 + var1, se),
        ~ 1 + var1
    )
    ## 2 variables are present in colData
    expect_identical(
        .checkScpModelFormula(~ 1 + var1 + var2, se),
        ~ 1 + var1 + var2
    )
    ## dot assignment works
    expect_identical(
        .checkScpModelFormula(~ 1 + ., se),
        ~ 1 + var1 + var2
    )
})

test_that(".removeResponseVariables", {
    ## response variable = warning
    ## response and no explanatory variable
    expect_warning(
        test <- .removeResponseVariables(y ~ 1, terms(y ~ 1)),
        "The formula contains a response variable and is ignored."
    )
    expect_identical(test, ~ 1)
    ## explanatory variable is retained
    expect_warning(
        test <- .removeResponseVariables(y ~ 1 + var1, terms(y ~ 1 + var1)),
        "The formula contains a response variable and is ignored."
    )
    expect_identical(test, ~ 1 + var1)
    ## dot "." variable is retained
    expect_warning(
        test <- .removeResponseVariables(y ~ 1 + ., terms(y ~ 1 + ., data = list(foo = 1))),
        "The formula contains a response variable and is ignored."
    )
    expect_identical(test, ~ 1 + .)
    ## There are multiple responses and explanatory variables
    expect_warning(
        test <- .removeResponseVariables(
            y + x ~ 1 + var1 + var2, terms(y + x ~ 1 + var1 + var2)
        ),
        regexp = "The formula contains a response variable"
    )
    expect_identical(test, ~ 1 + var1 + var2)
    ## When no response variable, input = output
    expect_identical(
        .removeResponseVariables(
            ~ 1 + var1 + var2, terms(~ 1 + var1 + var2)
        ),
        ~ 1 + var1 + var2
    )
})

test_that(".checkExplanatoryVariables", {
    ## No model variables = error
    f <- ~ NULL
    expect_error(
        .checkExplanatoryVariables(f, terms(f), "var1"),
        regexp = "You provided a formula with no variable to model."
    )
    ## The colData is empty = error
    f <- ~ 1 + var1
    expect_error(
        .checkExplanatoryVariables(f, terms(f), character()),
        regexp = "colData\\(object\\) is empty."
    )
    ## The colData is missing 1 variable = error
    f <- ~ 1 + var1 + var2 + var3
    expect_error(
        .checkExplanatoryVariables(f, terms(f), c("var1", "var2")),
        regexp = "missing one or more variables.*var3"
    )
    ## The colData is missing 2 variables = error
    f <- ~ 1 + var1 + var2 + var3 + var4
    expect_error(
        .checkExplanatoryVariables(f, terms(f), c("var1", "var2")),
        regexp = "missing one or more variables.*var3, var4.$"
    )
    ## 1 variable is present in colData
    f <- ~ 1 + var1
    expect_identical(
        .checkExplanatoryVariables(f, terms(f), c("var1", "var2")),
        ~ 1 + var1
    )
    ## 2 variables are present in colData
    f <- ~ 1 + var1 + var2
    expect_identical(
        .checkExplanatoryVariables(f, terms(f), c("var1", "var2")),
        ~ 1 + var1 + var2
    )
    ## dot assignment works
    f <- ~ 1 + .
    expect_identical(
        .checkExplanatoryVariables(
            f, terms(f, data = list(var1 = 1, var2 = 1)),
            c("var1", "var2")
        ),
        ~ 1 + var1 + var2
    )
})

test_that(".replaceDotVariable", {
    ## No dot = no effect
    expect_identical(
        .replaceDotVariable(letters[1:2], c()),
        letters[1:2]
    )
    expect_identical(
        ## .replaceDotVariable should not receive numeric but I here
        ## tested that modelVars is "untouched"
        .replaceDotVariable(1:2, c()),
        1:2
    )
    ## only dot = returns available vars
    expect_identical(
        .replaceDotVariable(".", letters[1:3]),
        letters[1:3]
    )
    ## dot + other non overlapping vars
    expect_identical(
        .replaceDotVariable(c(".", letters[1]), letters[2:3]),
        letters[1:3]
    )
    ## dot + other overlapping vars
    expect_identical(
        .replaceDotVariable(c(".", letters[1]), letters[1:3]),
        letters[1:3]
    )
})

test_that(".checkInputIndex", {
    ## Value has wrong type = error
    expect_error(
        .checkInputIndex(NA, c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .checkInputIndex(factor(1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .checkInputIndex(data.frame(assay = 1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    expect_error(
        .checkInputIndex(matrix(1), c("test1", "test2"), "mock"),
        regexp = "'mock' must be a character, numeric or logical"
    )
    ## Value points to multiple assays = error
    expect_error(
        .checkInputIndex(c(TRUE, TRUE), c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    expect_error(
        .checkInputIndex(1:2, c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    expect_error(
        .checkInputIndex(c("test1", "test2"), c("test1", "test2"), "mock"),
        regexp = "'mock' points to multiple input assays."
    )
    ## Value points to out of bound index = error
    expect_error(
        .checkInputIndex(c(FALSE, FALSE, TRUE), c("test1", "test2"), "mock"),
        regexp = "out of bounds"
    )
    expect_error(
        .checkInputIndex(3, c("test1", "test2"), "mock"),
        regexp = "out of bounds"
    )
    expect_error(
        .checkInputIndex("test3", c("test1", "test2"), "mock"),
        regexp = "'test3' not found."
    )
    ## Valid numeric = no conversion
    expect_identical(
        .checkInputIndex(1, c("test1", "test2"), "mock"),
        1
    )
    ## Convert from character
    expect_identical(
        .checkInputIndex("test2", c("test1", "test2"), "mock"),
        2L
    )
    ## Convert from logical
    expect_identical(
        .checkInputIndex(c(FALSE, TRUE), c("test1", "test2"), "mock"),
        2L
    )
})

test_that(".joinScpModelOutput", {
    require(SummarizedExperiment)
    se <- SummarizedExperiment()
    ## x is empty = error
    expect_error(
        .joinScpModelOutput(c(), se),
        "length.names.x.*not TRUE"
    )
    ## Complete case
    se <- SummarizedExperiment(assays = matrix(1, 2, 2))
    colnames(se) <- letters[1:2]
    x <- list(A = c(a = 1L, b = 3L), B = c(a = 2L, b = 4L))
    expect_identical(
        .joinScpModelOutput(x, se),
        matrix(1:4, 2, 2, dimnames = list(LETTERS[1:2], letters[1:2]))
    )
    ## With NA case
    x <- list(A = c(a = 1, b = NA), B = c(a = NA, b = 4))
    expect_identical(
        .joinScpModelOutput(x, se),
        matrix(c(1, NA, NA, 4), 2, 2, dimnames = list(LETTERS[1:2], letters[1:2]))
    )
})
