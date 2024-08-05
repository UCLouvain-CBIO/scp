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
            ScpModelFit()
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
    ## Test an intercept only model
    se <- .createMinimalData(nr = 1, nc = 6)
    se$condition <- as.factor(rep(1:2, each = 3))
    ## add missing values that make condition constant and hence will
    ## be dropped, leaving only the intercept
    assay(se)[, 1:3] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    expect_equal(
        scpModel(se)@scpModelFitList[[1]],
        new(
            "ScpModelFit",
            coefficients = structure(1, .Names = "(Intercept)"),
            residuals = structure(rep(0, 3), .Names = colnames(se)[4:6]),
            effects = List(),
            df = 2,
            var = 0,
            uvcov = matrix(0.3331, dimnames = list("(Intercept)", "(Intercept)")),
            levels = List()
        ),
        tolerance = 1E-3
    )
})

test_that(".runScpModel", {
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
        .runScpModel(se, "model", verbose = TRUE)
    )
    ## Test verbose
    expect_silent(.runScpModel(se, "model", verbose = FALSE))
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
    ## characters become factors
    se$condition <- LETTERS[1:5]
    test <- .checkAnnotations(se)
    expect_identical(
        levels(test$condition),
        LETTERS[1:5]
    )
    ## check that missing levels are dropped
    se$condition <- factor(LETTERS[1:5])
    levels(se$condition) <- LETTERS[1:7]
    test <- .checkAnnotations(se)
    expect_identical(
        levels(test$condition),
        LETTERS[1:5]
    )
})

test_that(".checkExperimentalDesignRank", {
    ## Singular design due to constant variable = error
    coldata <- DataFrame(condition = rep(2, 10))
    expect_error(
        .checkExperimentalDesignRank(~ 1 + condition, coldata),
        "The design matrix is .near. singular"
    )
    ## Singular design due to correlated variables = error
    coldata <- DataFrame(
        condition1 = factor(rep(1:2, each = 4)),
        condition2 = factor(rep(1:4, each = 2))
    )
    expect_error(
        .checkExperimentalDesignRank(~ 1 + condition1 + condition2, coldata),
        "The design matrix is .near. singular"
    )
    ## more params than observations = error
    coldata <- DataFrame(
        matrix(1:15, 3, 5, dimnames = list(NULL, paste0("foo", 1:5)))
    )
    expect_error(
        .checkExperimentalDesignRank(~ 1 + foo1 + foo2 + foo3 + foo4 + foo5, coldata),
        "The design matrix is underdetermined"
    )
    ## Test the check passes
    coldata <- DataFrame(condition = rep(1:2, 10))
    expect_identical(
        .checkExperimentalDesignRank(~ 1 + condition, coldata),
        NULL
    )
})

## ---- Adaptive modelling ----

test_that(".adaptModel", {
    ## Two observation or less = empty design
    exp <- matrix(nrow = 2, ncol = 0)
    attr(exp, "levels") <- List()
    expect_identical(
        .adaptModel(y = rep(1, 2), DataFrame(foo = 1:2), ~ 1 + foo),
        exp
    )
    ## no missing values: nrow coldata = nrow design
    expect_identical(
        nrow(.adaptModel(y = rep(1, 4), DataFrame(foo = 1:4), ~ 1 + foo)),
        4L
    )
    ## when missing values: nrow design = number of observations
    expect_identical(
        nrow(.adaptModel(y = c(NA, rep(1, 3)), DataFrame(foo = 1:4), ~ 1 + foo)),
        3L
    )
    ## Any variable not in the formula is dropped
    expect_identical(
        colnames(.adaptModel(y = c(NA, rep(1, 3)),
                             DataFrame(foo = 1:4, foo2 = 1), ~ 1 + foo)),
        c("(Intercept)", "foo")
    )
    ## Any constant variable is dropped
    ## factor
    expect_identical(
        colnames(.adaptModel(
            y = c(NA, rep(1, 3)),
            DataFrame(foo = factor(rep(1, 4))), ~ 1 + foo
        )),
        "(Intercept)"
    )
    ## numerical, cf issue: https://github.com/UCLouvain-CBIO/scp/issues/54
    expect_identical(
        colnames(.adaptModel(y = c(NA, rep(1, 3)),
                             DataFrame(foo = rep(1, 4)), ~ 1 + foo)),
        "(Intercept)"
    )
    ## Numerical variable are centred
    ## integer
    expect_identical(
        .adaptModel(y = c(NA, rep(1, 3)),
                    DataFrame(foo = 0:3), ~ 1 + foo)[, "foo"],
        structure(as.vector(scale(1:3)), .Names = 1:3)
    )
    ## numeric
    expect_identical(
        .adaptModel(y = c(NA, rep(1, 3)),
                    DataFrame(foo = c(0, 1.1, 2.2, 3.3)), ~ 1 + foo)[, "foo"],
        structure(as.vector(scale(c(1.1, 2.2, 3.3))), .Names = 1:3)
    )
    ## Factors levels are stored in the model matrix
    test <- .adaptModel(y = c(NA, rep(1, 3)),
                        DataFrame(foo = factor(1:4)), ~ 1 + foo)
    expect_identical(
        attr(test, "levels"),
        List(foo = as.character(2:4))
    )
    ## Factors are encoded as sum contrasts
    expect_identical(
        attr(test, "contrasts"),
        list(foo = "contr.sum")
    )
    ## Check model design content
    test <- .adaptModel(
        y = c(NA, rep(1, 3), NA, NA),
        coldata = DataFrame(numeric = 1:6,
                            dropped = 1,
                            twoLevels = factor(rep(1:2, 3)),
                            moreLevels = factor(rep(1:3, 2)),
                            moreLevelsButDropped = factor(c(1, 2, 2, 2, 1, 1)),
                            moreLevelsButAdapted = factor(c(3, 3, 2, 2, 1, 1))),
        formula =  ~ 1 + numeric + twoLevels + moreLevels + moreLevelsButDropped +moreLevelsButAdapted
    )
    exp <- as.matrix(data.frame(
        "(Intercept)" = 1,
        numeric = scale(2:4),
        twoLevels1 = c(-1, 1, -1),
        moreLevels1 = c(0, -1, 1),
        moreLevels2 = c(1, -1, 0),
        moreLevelsButAdapted1 = c(-1, 1, 1),
        check.names = FALSE,
        row.names = 1:3
    ))
    attr(exp, "assign") <- c(0L, 1L, 2L, 3L, 3L, 4L)
    attr(exp, "contrasts") <- list(
        twoLevels = "contr.sum",
        moreLevels = "contr.sum",
        moreLevelsButAdapted = "contr.sum"
    )
    attr(exp, "levels") <- List(twoLevels = as.character(1:2),
                                moreLevels = as.character(1:3),
                                moreLevelsButAdapted = as.character(2:3))
    expect_identical(
        exp,
        test
    )
})

test_that(".formatCategoricalVariablesAsFactor", {
    ## No categorical: input = output
    x <- DataFrame(foo = 1:10)
    expect_identical(
        .formatCategoricalVariablesAsFactor(x),
        x
    )
    ## logical
    x <- DataFrame(foo = rep(c(TRUE, FALSE), 3))
    factorVar <- factor(rep(c(TRUE, FALSE), 3))
    expect_identical(
        .formatCategoricalVariablesAsFactor(x),
        DataFrame(foo = factorVar)
    )
    ## character
    x <- DataFrame(foo = rep(c("bar1", "bar2"), 3))
    factorVar <- factor(rep(c("bar1", "bar2"), 3))
    expect_identical(
        .formatCategoricalVariablesAsFactor(x),
        DataFrame(foo = factorVar)
    )
    ## factor
    x <- DataFrame(foo = factor(rep(c("bar1", "bar2"), 3)))
    expect_identical(
        .formatCategoricalVariablesAsFactor(x),
        x
    )
    ## factor but drop levels
    x <- DataFrame(foo = factor(rep(c("bar1", "bar2"), 3),
                                levels = c("bar1", "bar2", "bar3")))
    expect_identical(
        .formatCategoricalVariablesAsFactor(x),
        DataFrame(foo = factor(rep(c("bar1", "bar2"), 3),
                               levels = c("bar1", "bar2")))
    )
})

test_that(".scaleNumericalVariables", {
    ## Not enough rows: input = output
    x <- DataFrame(foo = 1)
    expect_identical(
        .scaleNumericalVariables(x),
        x
    )
    ## Only categorical: input = output
    x <- DataFrame(foo = rep("bar1", 5))
    expect_identical(
        .scaleNumericalVariables(x),
        x
    )
    ## integer
    expect_identical(
        .scaleNumericalVariables(DataFrame(foo = 1L:5L)),
        DataFrame(foo = as.vector(scale(1L:5L)))
    )
    ## numeric
    expect_identical(
        .scaleNumericalVariables(DataFrame(foo = c(0, 1.1, 2.2, 3.3))),
        DataFrame(foo = as.vector(scale(c(0, 1.1, 2.2, 3.3))))
    )
    ## multiple numerical
    expect_identical(
        .scaleNumericalVariables(
            DataFrame(foo1 = c(0, 1.1, 2.2, 3.3),
                      foo2 = c(0, 1.1, 2.2, 3.3))
        ),
        DataFrame(foo1 = as.vector(scale(c(0, 1.1, 2.2, 3.3))),
                  foo2 = as.vector(scale(c(0, 1.1, 2.2, 3.3))))
    )
    ## numerical + categorical
    expect_identical(
        .scaleNumericalVariables(
            DataFrame(foo1 = c(0, 1.1, 2.2, 3.3),
                      foo2 = "bar")
        ),
        DataFrame(foo1 = as.vector(scale(c(0, 1.1, 2.2, 3.3))),
                  foo2 = "bar")
    )
})

test_that(".dropConstantVariablesFromFormula", {
    ## No constant variable: input = output
    x <- DataFrame(foo1 = 1:10, foo2 = rep(c("bar1", "bar2"), 5))
    formula <- ~ 1 + foo1 + foo2
    expect_identical(
        .dropConstantVariablesFromFormula(x, formula),
        formula
    )
    ## Constant numerical
    x <- DataFrame(foo1 = 1, foo2 = rep(c("bar1", "bar2"), 5))
    expect_identical(
        .dropConstantVariablesFromFormula(x, ~ 1 + foo1 + foo2),
        ~ 1 + foo2
    )
    ## Constant categorical
    x <- DataFrame(foo1 = 1:10, foo2 = factor("bar1"))
    expect_identical(
        .dropConstantVariablesFromFormula(x, ~ 1 + foo1 + foo2),
        ~ 1 + foo1
    )
    ## Constant numerical + categorical
    x <- DataFrame(foo1 = 1:10, foo2 = factor("bar1"), foo3 = 1)
    expect_identical(
        .dropConstantVariablesFromFormula(x, ~ 1 + foo1 + foo2 + foo3),
        ~ 1 + foo1
    )
    ## All variables are constant = intercept only model
    x <- DataFrame(foo1 = rep(1, 10), foo2 = factor(rep("bar1", 10)))
    expect_identical(
        .dropConstantVariablesFromFormula(x, ~ 1 + foo1 + foo2),
        ~ 1
    )
})

test_that(".modelContrasts", {
    ## No factor = return empty list
    expect_identical(
        .modelContrasts(DataFrame(foo = 1:10)),
        NULL
    )
    ## 1 factor
    expect_identical(
        .modelContrasts(DataFrame(foo = factor("bar"))),
        list(foo = "contr.sum")
    )
    ## 2 factors
    expect_identical(
        .modelContrasts(DataFrame(foo = factor("bar"), foo2 = factor("bar"))),
        list(foo = "contr.sum", foo2 = "contr.sum")
    )
    ## 1 factor + 1 other
    expect_identical(
        .modelContrasts(DataFrame(foo = factor("bar"), foo2 = 1:10)),
        list(foo = "contr.sum")
    )
})

test_that(".modelLevels", {
    ## No factors = empty list
    expect_identical(
        .modelLevels(DataFrame(foo = 1:10)),
        List()
    )
    ## 1 factor
    expect_identical(
        .modelLevels(DataFrame(foo = factor(rep(c("bar1", "bar2", "bar3"), 4)))),
        List(foo = c("bar1", "bar2", "bar3"))
    )
    expect_identical(
        .modelLevels(DataFrame(foo = factor(rep(c("bar1", "bar2"), 4)))),
        List(foo = c("bar1", "bar2"))
    )
    ## 2 factors
    expect_identical(
        .modelLevels(DataFrame(
            foo = factor(rep(c("bar1", "bar2", "bar3"), 4)),
            foo2 = factor(rep(c("bar4", "bar5", "bar6"), 4))
        )),
        List(foo = c("bar1", "bar2", "bar3"),
             foo2 = c("bar4", "bar5", "bar6"))
    )
    ## 1 factor + 1 other
    expect_identical(
        .modelLevels(DataFrame(foo = factor(c("bar1", "bar2")), foo2 = 1:10)),
        List(foo = c("bar1", "bar2"))
    )
})

## ---- Model estimation ----

test_that(".fitModel", {
    ## Empty design matrix = empty ScpModelFit object
    expect_identical(
        .fitModel(1:10, matrix(nrow = 10, ncol = 0), "foo"),
        ScpModelFit()
    )
    ## Test estimation
    y <- structure(rep(1, 10), .Names = letters[1:10])
    x <- matrix(c(rep(1, 10), rep(1, 5), rep(-1, 5)), ncol = 2,
                dimnames = list(letters[1:10], c("(Intercept)", "foo")))
    attr(x, "levels") <- List(foo = c("A", "B"))
    test <- .fitModel(y, x, "foo")

    expect_equal(
        test@coefficients,
        structure(c(1, 0), .Names = c("(Intercept)", "foo")),
        tolerance = 1E-3
    )
    expect_equal(
        test@residuals,
        structure(rep(0, 10), .Names = letters[1:10]),
        tolerance = 1E-3
    )
    expect_equal(
        test@effects,
        List(foo = structure(rep(0, 10), .Names = letters[1:10])),
        tolerance = 1E-3
    )
    expect_identical(
        test@df,
        8
    )
    expect_equal(
        test@var,
        0,
        tolerance = 1E-3
    )
    expect_equal(
        test@uvcov,
        matrix(c(0.1, 0, 0, 0.1), ncol = 2,
               dimnames = list(c("(Intercept)", "foo"), c("(Intercept)", "foo"))),
        tolerance = 1E-3
    )
    expect_identical(
        test@levels,
        List(foo = c("A", "B"))
    )
})

test_that(".fitRidge", {
    y <- structure(rep(1, 10), .Names = letters[1:10])
    x <- matrix(c(rep(1, 10), rep(1, 5), rep(-1, 5)), ncol = 2,
                dimnames = list(letters[1:10], c("(Intercept)", "foo")))
    test <- .fitRidge(y, x)
    expect_equal(
        test$coefficients,
        structure(c(1, 0), .Names = c("(Intercept)", "foo")),
        tolerance = 1E-3
    )
    expect_equal(
        test$residuals,
        structure(rep(0, 10), .Names = letters[1:10]),
        tolerance = 1E-3
    )
    expect_identical(
        test$df,
        8
    )
    expect_equal(
        test$var,
        0,
        tolerance = 1E-3
    )
    expect_equal(
        test$uvcov,
        matrix(c(0.1, 0, 0, 0.1), ncol = 2,
               dimnames = list(c("(Intercept)", "foo"), c("(Intercept)", "foo"))),
        tolerance = 1E-3
    )
})

test_that(".computeModelEffect", {
    df <- data.frame(
        fooClass = rep(LETTERS[1:3], 5),
        fooNum = 1:15,
        row.names = LETTERS[1:15]
    )
    ## Single numerical variable
    x <- model.matrix(~ 1 + fooNum, df)
    beta <- c("(Intercept)" = 2, fooNum = 1)
    expect_identical(
        .computeModelEffect(beta, x, effect = "fooNum"),
        (beta["fooNum", drop = FALSE] %*% t(x[, "fooNum", drop = FALSE]))[1, ]
    )
    ## Single categorical variable (3 levels)
    x <- model.matrix(~ 1 + fooClass, df, contrasts.arg = list(fooClass = "contr.sum"))
    beta <- c("(Intercept)" = 2, fooClass1 = 1, fooClass2 = 2)
    sel <- c("fooClass1", "fooClass2")
    expect_identical(
        .computeModelEffect(beta, x, effect = "fooClass"),
        (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ]
    )
    ## Categorical and numerical variables
    x <- model.matrix(
        ~ 1 + fooClass + fooNum, df,
        contrasts.arg = list(fooClass = "contr.sum")
    )
    beta <- c("(Intercept)" = 2, fooClass1 = 1, fooClass2 = 2, fooNum = 1)
    sel <- c("fooClass1", "fooClass2")
    expect_identical(
        .computeModelEffect(beta, x, effect = "fooClass"),
        (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ]
    )
    sel <- "fooNum"
    expect_identical(
        .computeModelEffect(beta, x, effect = "fooNum"),
        (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ]
    )
    ## Effect is absent = empty effects
    expect_identical(
        .computeModelEffect(beta, x, effect = "fooBar"),
        NULL
    )
    ## Works also for intercept. Although it should never be used, it
    ## is a test for variables with "()"
    sel <- "(Intercept)"
    expect_identical(
        .computeModelEffect(beta = beta, x, effect = "(Intercept)"),
        (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ]
    )
})

test_that(".computeModelEffects", {
    require("S4Vectors")
    df <- data.frame(
        fooClass = rep(LETTERS[1:3], 5),
        fooNum = 1:15,
        row.names = LETTERS[1:15]
    )
    ## Single variable
    x <- model.matrix(~ 1 + fooNum, df)
    beta <- c("(Intercept)" = 2, fooNum = 1)
    expect_identical(
        .computeModelEffects(beta = beta, x, effects = "fooNum"),
        List(fooNum = (beta["fooNum", drop = FALSE] %*% t(x[, "fooNum", drop = FALSE]))[1, ])
    )
    ## 2 variables
    x <- model.matrix(
        ~ 1 + fooClass + fooNum, df,
        contrasts.arg = list(fooClass = "contr.sum")
    )
    beta <- c("(Intercept)" = 2, fooClass1 = 1, fooClass2 = 2, fooNum = 1)
    sel <- c("fooClass1", "fooClass2")
    expect_identical(
        .computeModelEffects(beta = beta, x, effects = c("fooClass", "fooNum")),
        List(fooClass = (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ],
             fooNum = (beta["fooNum", drop = FALSE] %*% t(x[, "fooNum", drop = FALSE]))[1, ])
    )
    ## Effect is absent = empty List
    exp <- List()
    names(exp) <- character(0)
    expect_identical(
        .computeModelEffects(beta = beta, x, effects = "fooBar"),
        exp
    )
    ## One out of two effects are absent = missing effect is ignored
    expect_identical(
        .computeModelEffects(beta = beta, x, effects = c("fooClass", "fooBar")),
        List(fooClass = (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ])
    )
    ## Works also for intercept. Although it should never be used, it
    ## is a test for variables with "()"
    expect_identical(
        .computeModelEffects(beta = beta, x, effects = c("fooClass", "fooNum", "(Intercept)")),
        List(fooClass = (beta[sel, drop = FALSE] %*% t(x[, sel, drop = FALSE]))[1, ],
             fooNum = (beta["fooNum", drop = FALSE] %*% t(x[, "fooNum", drop = FALSE]))[1, ],
             "(Intercept)" = (beta["(Intercept)", drop = FALSE] %*% t(x[, "(Intercept)", drop = FALSE]))[1, ])
    )
    ## No effects = empty list
    expect_identical(
        .computeModelEffects(beta = beta, x = x, effects = c()),
        List()
    )
})

## ---- ScpModel filter plot ----

test_that("scpModelFilterPlot", {
    require(SummarizedExperiment)
    se <- .createMinimalData(10, 10)
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 2
    set.seed(1234)
    for (i in 1:nrow(se)) {
        assay(se)[i, sample(1:ncol(se), i)] <- NA
    }
    ## Standard case
    se <- scpModelWorkflow(object = se, formula = ~ 1 + condition)
    expect_message(
        pl <- scpModelFilterPlot(object = se),
        "To change the threshold, use:\nscpModelFilterThreshold"
    )
    vdiffr::expect_doppelganger(
        "scpModelFilterPlot for model1",
        pl
    )
    ## After adding threshold
    scpModelFilterThreshold(object = se) <- 2
    expect_message(
        pl <- scpModelFilterPlot(object = se),
        "To change the threshold, use:\nscpModelFilterThreshold"
    )
    vdiffr::expect_doppelganger(
        "scpModelFilterPlot after adding threshold",
        pl
    )
    ## Works with mulitple models
    se <- scpModelWorkflow(object = se, formula = ~ 1 + condition, name = "model2")
    expect_message(
        pl <- scpModelFilterPlot(object = se, name = "model2"),
        "To change the threshold, use:\nscpModelFilterThreshold"
    )
    vdiffr::expect_doppelganger(
        "scpModelFilterPlot for model2",
        pl
    )
})

test_that(".filterSubtitle", {
    ## Default filter threshold
    expect_identical(
        .filterSubtitle(npRatio = 0:10, threshold = 0),
        "Total features: 11 \nEstimated features (n/p >= 1): 10"
    )
    ## Filter removes about half the features threshold
    expect_identical(
        .filterSubtitle(npRatio = 0:10, threshold = 5),
        "Total features: 11 \nEstimated features (n/p >= 1): 10\nSelected features (n/p >= 5): 6"
    )
})

test_that(".filterPlot", {
    ## Default filter threshold
    vdiffr::expect_doppelganger(
        ".filterPlot with default filter threshold",
        .filterPlot(npRatio = 0:10, threshold = 0)
    )
    ## Filter removes about half the features threshold
    vdiffr::expect_doppelganger(
        ".filterPlot with filtering out half of the features",
        .filterPlot(npRatio = 0:10, threshold = 5)
    )
})
