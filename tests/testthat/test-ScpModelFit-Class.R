## ---- Constructor ----

test_that("ScpModelFit", {

    ## Check slots are correctly inialized
    x <- ScpModelFit()
    expect_identical(x@coefficients, numeric())
    expect_identical(x@df, numeric())
    expect_identical(x@var, numeric())
    expect_identical(x@uvcov, matrix(numeric(), 0, 0))
    expect_identical(x@levels, List())
})

## ---- Getters ----

test_that("scpModelFitCoefficients", {
    x <- ScpModelFit()
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

test_that("scpModelFitDf", {
    x <- ScpModelFit()
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
    x <- ScpModelFit()
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
    x <- ScpModelFit()
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
    x <- ScpModelFit()
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
    x <- ScpModelFit()
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
    x <- ScpModelFit()
    ## value is wrong type = error
    expect_error(
        scpModelFitCoefficients(x) <- rep("foo", 5),
        "an object of class .character. is not valid for @.coefficients"
    )
    expect_error(
        scpModelFitCoefficients(x) <- matrix(1:5),
        "an object of class.*matrix.*is not valid for @.coefficients"
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

test_that("scpModelFitDf<-", {
    x <- ScpModelFit()
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
        "an object of class.*matrix.*is not valid for @.df"
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
    x <- ScpModelFit()
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
        "an object of class.*matrix.*is not valid for @.var"
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
    x <- ScpModelFit()
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
    x <- ScpModelFit()
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
