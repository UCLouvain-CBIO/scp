## ---- Constructor ----

test_that("ScpModelFit", {
    ## n slot is required
    expect_error(
        ScpModelFit(p = 10L),
        "argument \"n\" is missing, with no default"
    )
    ## n must be integer
    expect_error(
        ScpModelFit(n = 10, p = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(n = "foo", p = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(n = list(1L), p = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(n = matrix(1), p = 10L),
        "should be or extend class \"integer\""
    )
    ## n must be positive
    expect_error(
        ScpModelFit(n = -1L, p = 10L),
        "n >= 0 is not TRUE"
    )
    ## n must have length 1
    expect_error(
        ScpModelFit(n = integer(), p = 10L),
        "length\\(n\\) == 1 is not TRUE"
    )
    expect_error(
        ScpModelFit(n = 1:2, p = 10L),
        "length\\(n\\) == 1 is not TRUE"
    )
    ## p slot is required
    expect_error(
        ScpModelFit(n = 10L),
        "argument \"p\" is missing, with no default"
    )
    ## p must be integer
    expect_error(
        ScpModelFit(p = 10, n = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(p = "foo", n = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(p = list(1L), n = 10L),
        "should be or extend class \"integer\""
    )
    expect_error(
        ScpModelFit(p = matrix(1), n = 10L),
        "should be or extend class \"integer\""
    )
    ## n must be positive
    expect_error(
        ScpModelFit(p = -1L, n = 10L),
        "p >= 0 is not TRUE"
    )
    ## n must have length 1
    expect_error(
        ScpModelFit(p = integer(), n = 10L),
        "length\\(p\\) == 1 is not TRUE"
    )
    expect_error(
        ScpModelFit(p = 1:2, n = 10L),
        "length\\(p\\) == 1 is not TRUE"
    )
    ## Check slots are correctly inialized
    x <- ScpModelFit(10L, 10L)
    expect_identical(x@n, 10L)
    expect_identical(x@p, 10L)
    expect_identical(x@coefficients, numeric())
    expect_identical(x@residuals, numeric())
    expect_identical(x@effects, List())
    expect_identical(x@df, numeric())
    expect_identical(x@var, numeric())
    expect_identical(x@uvcov, matrix(numeric(), 0, 0))
    expect_identical(x@levels, List())
})

## ---- Getters ----

test_that("scpModelFitN", {
    x <- ScpModelFit(10L, 5L)
    expect_identical(
        scpModelFitN(x),
        10L
    )
})

test_that("scpModelFitP", {
    x <- ScpModelFit(10L, 5L)
    expect_identical(
        scpModelFitP(x),
        5L
    )
})

test_that("scpModelFitCoefficients", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitCoefficients(x),
        numeric()
    )
    ## after assignment
    x@coefficients <- 1:10
    expect_identical(
        scpModelFitCoefficients(x),
        1:10
    )
})

test_that("scpModelFitResiduals", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitResiduals(x),
        numeric()
    )
    ## after assignment
    x@residuals <- 1:10
    expect_identical(
        scpModelFitResiduals(x),
        1:10
    )
})

test_that("scpModelFitEffects", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitEffects(x),
        List()
    )
    ## after assignment
    x@effects <- List(effect = 1:10)
    expect_identical(
        scpModelFitEffects(x),
        List(effect = 1:10)
    )
})

test_that("scpModelFitDf", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitDf(x),
        numeric()
    )
    ## after assignment
    x@df <- 10
    expect_identical(
        scpModelFitDf(x),
        10
    )
})

test_that("scpModelFitVar", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitVar(x),
        numeric()
    )
    ## after assignment
    x@var <- 10
    expect_identical(
        scpModelFitVar(x),
        10
    )
})

test_that("scpModelFitUvcov", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitUvcov(x),
        matrix(numeric(), 0, 0)
    )
    ## after assignment
    x@uvcov <- matrix(1:10, 2)
    expect_identical(
        scpModelFitUvcov(x),
        matrix(1:10, 2)
    )
})

test_that("scpModelFitVcov", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitVcov(x),
        matrix(numeric(), 0, 0)
    )
    ## after assignment
    x@uvcov <- matrix(1:10, 2)
    x@var <- 2
    expect_identical(
        scpModelFitVcov(x),
        matrix(1:10, 2) * 2
    )
})

test_that("scpModelFitLevels", {
    x <- ScpModelFit(10L, 5L)
    ## default
    expect_identical(
        scpModelFitLevels(x),
        List()
    )
    ## after assignment
    x@levels <- List(a = 1:10)
    expect_identical(
        scpModelFitLevels(x),
        List(a = 1:10)
    )
})

## ---- Setters ----

test_that("scpModelFitCoefficients<-", {
    x <- ScpModelFit(10L, 5L)
    ## value has wrong length = length is not p = error
    expect_error(
        scpModelFitCoefficients(x) <- 1:4,
        "length\\(value\\) == scpModelFitP\\(object\\) is not TRUE"
    )
    ## value is wrong type = error
    expect_error(
        scpModelFitCoefficients(x) <- rep("foo", 5),
        "an object of class .character. is not valid for @.coefficients"
    )
    expect_error(
        scpModelFitCoefficients(x) <- matrix(1:5),
        "an object of class .matrix. is not valid for @.coefficients"
    )
    expect_error(
        scpModelFitCoefficients(x) <- list(1, 2, 3, 4, 5),
        "an object of class .list. is not valid for @.coefficients"
    )
    ## Correct usage
    scpModelFitCoefficients(x) <- 1:5
    expect_identical(
        x@coefficients,
        1:5
    )
})

test_that("scpModelFitResiduals<-", {
    x <- ScpModelFit(10L, 5L)
    ## value has wrong length = length is not p = error
    expect_error(
        scpModelFitResiduals(x) <- 1:4,
        "length\\(value\\) == scpModelFitN\\(object\\) is not TRUE"
    )
    ## value is wrong type = error
    expect_error(
        scpModelFitResiduals(x) <- rep("foo", 10),
        "an object of class .character. is not valid for @.residuals"
    )
    expect_error(
        scpModelFitResiduals(x) <- matrix(1:10),
        "an object of class .matrix. is not valid for @.residuals"
    )
    expect_error(
        scpModelFitResiduals(x) <- list(1, 2, 3, 4, 5, 6, 7, 8, 9 , 10),
        "an object of class .list. is not valid for @.residuals"
    )
    ## Correct usage
    scpModelFitResiduals(x) <- 1:10
    expect_identical(
        x@residuals,
        1:10
    )
    ## Correct usage with empty array
    x <- ScpModelFit(0L, 5L)
    scpModelFitResiduals(x) <- numeric()
    expect_identical(
        x@residuals,
        numeric()
    )
})

test_that("scpModelFitEffects<-", {
    x <- ScpModelFit(10L, 5L)
    ## value has elements that are not a numeric vector = error
    el1 <- c(a = 1, b = 2, c = 3)
    el2 <- matrix(1:3)
    expect_error(
        scpModelFitEffects(x) <- List(el1, el2),
        "all.*value.*inherits.*numeric.*is not TRUE"
    )
    el2 <- list(1:3)
    expect_error(
        scpModelFitEffects(x) <- List(el1, el2),
        "all.*value.*inherits.*numeric.*is not TRUE"
    )
    ## value has elements with inconsistent names = error
    el2 <- c(a = 4, b = 5, d = 6)
    expect_error(
        scpModelFitEffects(x) <- List(el1, el2),
        "Effect vectors do not share identical names."
    )
    ## value has elements that are inconsistent with residual names = error
    x@residuals <- c(a = 7, b = 8, d = 9)
    el2 <- c(a = 4, b = 5, c = 6)
    expect_error(
        scpModelFitEffects(x) <- List(el1, el2),
        "Effects and residuals do not share identical names."
    )
    ## value has wrong type = error
    x@residuals <- c(a = 7, b = 8, c = 9)
    expect_error(
        scpModelFitEffects(x) <- list(el1, el2),
        "list. is not valid for @.effects"
    )
    ## Correct usage
    expect_identical(
        scpModelFitEffects(x) <- List(el1, el2),
        List(el1, el2)
    )
    ## Correct usage with empty list
    x@residuals <- numeric()
    expect_identical(
        scpModelFitEffects(x) <- List(),
        List()
    )
})

test_that("scpModelFitDf<-", {
    x <- ScpModelFit(10L, 5L)
    ## value has not length 1 = error
    expect_error(
        scpModelFitDf(x) <- 1:4,
        "length\\(value\\) == 1 is not TRUE"
    )
    ## value is wrong type = error
    expect_error(
        scpModelFitDf(x) <- "foo",
        "an object of class .character. is not valid for @.df"
    )
    expect_error(
        scpModelFitDf(x) <- matrix(1),
        "an object of class .matrix. is not valid for @.df"
    )
    expect_error(
        scpModelFitDf(x) <- list(1),
        "an object of class .list. is not valid for @.df"
    )
    ## value must be positive
    expect_error(
        scpModelFitDf(x) <- -1,
        "value >= 0 is not TRUE"
    )
    ## Correct usage
    scpModelFitDf(x) <- 5.2
    expect_identical(
        x@df,
        5.2
    )
})

test_that("scpModelFitVar<-", {
    x <- ScpModelFit(10L, 5L)
    ## value has not length 1 = error
    expect_error(
        scpModelFitVar(x) <- 1:4,
        "length\\(value\\) == 1 is not TRUE"
    )
    ## value is wrong type = error
    expect_error(
        scpModelFitVar(x) <- "foo",
        "an object of class .character. is not valid for @.var"
    )
    expect_error(
        scpModelFitVar(x) <- matrix(1),
        "an object of class .matrix. is not valid for @.var"
    )
    expect_error(
        scpModelFitVar(x) <- list(1),
        "an object of class .list. is not valid for @.var"
    )
    ## value must be positive
    expect_error(
        scpModelFitVar(x) <- -1,
        "value >= 0 is not TRUE"
    )
    ## Correct usage
    scpModelFitVar(x) <- 3.1
    expect_identical(
        x@var,
        3.1
    )
})

test_that("scpModelFitUvcov<-", {
    x <- ScpModelFit(10L, 5L)
    x@coefficients <- c(a = 1, b = 2, c = 3)
    ## value has no rownames  = error
    uvcov <- matrix(9, 3, 3, dimnames = list(NULL, letters[1:3]))
    expect_error(
        scpModelFitUvcov(x) <- uvcov,
        "identical.*names.*coef.*rownames.value.*is not TRUE"
    )
    ## value has no colnames  = error
    uvcov <- matrix(9, 3, 3, dimnames = list(letters[1:3], NULL))
    expect_error(
        scpModelFitUvcov(x) <- uvcov,
        "identical.*names.*coef.*colnames.value.*is not TRUE"
    )
    ## rownames in value do not match names in coefficients = error
    uvcov <- matrix(9, 3, 3, dimnames = list(letters[4:6], letters[1:3]))
    expect_error(
        scpModelFitUvcov(x) <- uvcov,
        "identical.*names.*coef.*rownames.value.*is not TRUE"
    )
    ## colnames in value do not match names in coefficients = error
    uvcov <- matrix(9, 3, 3, dimnames = list(letters[1:3], letters[4:6]))
    expect_error(
        scpModelFitUvcov(x) <- uvcov,
        "identical.*names.*coef.*colnames.value.*is not TRUE"
    )
    ## value has wrong type = error
    uvcov <- data.frame(matrix(9, 3, 3, dimnames = list(letters[1:3], letters[1:3])))
    expect_error(
        scpModelFitUvcov(x) <- uvcov,
        "data.frame. is not valid for @.uvcov"
    )
    ## Correct usage
    uvcov <- matrix(9, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
    scpModelFitUvcov(x) <- uvcov
    expect_identical(
        x@uvcov,
        uvcov
    )
})

test_that("scpModelFitLevels<-", {
    x <- ScpModelFit(10L, 5L)
    ## value elements have wrong type = error
    expect_error(
        scpModelFitLevels(x) <- List(rep("foo", 5), 1:5),
        "all.*value.*inherits.*character.*is not TRUE"
    )
    ## value is not named = error
    expect_error(
        scpModelFitLevels(x) <- List(rep("foo", 5), rep("foo", 5)),
        "List of levels must be named."
    )
    ## value has wrong type = error
    x@effects <- List(var1 = matrix())
    expect_error(
        scpModelFitLevels(x) <- list(var1 = rep("foo", 5)),
        "an object of class .list. is not valid for @.levels"
    )
    ## Element names are not present in effects
    expect_error(
        scpModelFitLevels(x) <- List(var1 = rep("foo", 5), var2 = rep("foo", 5)),
        "Some levels are not matched to effects."
    )
    expect_error(
        scpModelFitLevels(x) <- List(var2 = rep("foo", 5)),
        "Some levels are not matched to effects."
    )
    ## Correct usage
    scpModelFitLevels(x) <- List(var1 = rep("foo", 5))
    expect_identical(
        x@levels,
        List(var1 = rep("foo", 5))
    )
    x@effects <- List(var1 = matrix(), var2 = rep("foo", 5))
    scpModelFitLevels(x) <- List(var1 = rep("foo", 5))
    expect_identical(
        x@levels,
        List(var1 = rep("foo", 5))
    )
    scpModelFitLevels(x) <- List(var1 = rep("foo", 5), var2 = rep("foo", 5))
    expect_identical(
        x@levels,
        List(var1 = rep("foo", 5), var2 = rep("foo", 5))
    )
    ## Also works when no levels = empty List
    scpModelFitLevels(x) <- List()
    expect_identical(
        x@levels,
        List()
    )
})
