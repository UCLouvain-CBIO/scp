library("vdiffr")
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

test_that("scpDifferentialAnalysis", {
    se <- .createMinimalData(nc = 6)
    ## Both contrast and coefficient is missing = error
    expect_error(
        scpDifferentialAnalysis(se),
        "'contrasts' and 'coefficients' cannot be both NULL."
    )
    ## Missing name is identical to default name
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    estimable <- c("a", "c", "d", "e", "h", "i", "j")
    expect_identical(
        scpDifferentialAnalysis(
            se, contrast = list(c("condition", "1", "2"))
        ),
        scpDifferentialAnalysis(
            se, contrast = list(c("condition", "1", "2")), name = "model"
        )
    )

    ## Test contrasts with 2 levels
    exp <- List(condition_1_vs_2 = DataFrame(
        feature = estimable,
        Estimate = structure((0:9), .Names = rownames(se))[estimable],
        SE = structure(rep(0, 7), .Names = estimable),
        Df = rowSums(!is.na(assay(se)[estimable, ])) - 2,
        tstatistic = structure(c(0.7072246, 3452.7200380, 1206.4382131, 1778.6625013, 1990.7885585, 3534.6922866, 2058.9077196), .Names = estimable),
        pvalue = structure(c(0.5527, rep(0, 6)), .Names = estimable),
        padj = structure(c(0.5527, rep(0, 6)), .Names = estimable)
    ))
    metadata(exp[[1]])$contrast <- c("condition", "1", "2")
    expect_equal(
        scpDifferentialAnalysis(
            se, contrast = list(c("condition", "1", "2"))
        ),
        exp,
        tolerance = 1E-2
    )
    ## Test 2 contrasts with 3 leveled-group
    se <- .createMinimalData(nc = 9)
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    estimable <- c("c", "d", "f", "g", "h", "j")
    c1vs2 <- DataFrame(
        feature = estimable,
        Estimate = structure((0:9), .Names = rownames(se))[estimable],
        SE = structure(rep(0, 6), .Names = estimable),
        Df = rowSums(!is.na(assay(se)[estimable, ])) - 3,
        tstatistic = structure(c(824.8086, 1777.5883, 1975.3086, 1088.6531, 4759.9688, 2614.4243), .Names = estimable),
        pvalue = structure(rep(0, 6), .Names = estimable),
        padj = structure(rep(0, 6), .Names = estimable)
    )
    metadata(c1vs2)$contrast <- c("condition", "1", "2")
    estimable <- c("a", "b", "c", "d", "e", "f", "g", "h", "j")
    c1vs3 <- DataFrame(
        feature = estimable,
        Estimate = structure(rep(0, 9), .Names = estimable),
        SE = structure(rep(0, 9), .Names = estimable),
        Df = sapply(scpModelFitList(se), scpModelFitDf)[estimable],
        tstatistic = structure(c(-0.7072246, 0.3333611, -0.2378374, -0.3530348, 0.3333611, -0.5921999, -0.2115692, -0.7239608, -0.3951838), .Names = estimable),
        pvalue = structure(c(0.5527268, 0.7951513, 0.8513498, 0.7474114, 0.7951513, 0.6137490, 0.8672682, 0.5091646, 0.7308728), .Names = estimable),
        padj = structure(rep(0.8672, 9), .Names = estimable)
    )
    metadata(c1vs3)$contrast <- c("condition", "1", "3")
    expect_equal(
        scpDifferentialAnalysis(
            se,
            contrast = list(
                c("condition", "1", "2"), ## there is change = positive control
                c("condition", "1", "3") ## there is no change = negative control
            )
        ),
        List(condition_1_vs_2 = c1vs2,
             condition_1_vs_3 = c1vs3),
        tolerance = 1E-2
    )
    ## Test 1 coefficient
    se <- .createMinimalData(nc = 9)
    set.seed(1234)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se$numeric <- 1:ncol(se)
    se <- scpModelWorkflow(se, formula = ~ 1 + numeric, name = "model2")
    estimable <- c("a", "b", "c", "d", "e", "g", "h", "i", "j")
    expect_equal(
        scpDifferentialAnalysis(
            se,
            coefficients = "numeric",
            name = "model2"
        ),
        List(numeric = DataFrame(
            feature = estimable,
            Estimate = structure(rep(0, 9), .Names = estimable),
            SE = structure(rep(0, 9), .Names = estimable),
            Df = sapply(scpModelFitList(se)[estimable], scpModelFitDf),
            tstatistic = structure(rep(0, 9), .Names = estimable),
            pvalue = structure(rep(1, 9), .Names = estimable),
            padj = structure(rep(1, 9), .Names = estimable)
        )),
        tolerance = 1E-2
    )
    ## Test contrast and coefficient
    ## Use this case also to assess positive controls. Simulated data
    ## with different intercepts, different slopes for numeric and different FC for contrasts
    ## Use more sample to better test the accuracy
    se <-  SummarizedExperiment(assays = List(
        matrix(0, 10, 100, dimnames = list(paste0("row", 1:10),
                                           paste0("col", 1:100)))
    ))
    ## add missing values
    set.seed(1234)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/5)] <- NA
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    se$numeric <- scale(1:ncol(se))
    for (i in 0:9) {
        ## add group effect
        assay(se)[i+1, se$condition == 2] <- assay(se)[i + 1, se$condition == 2] + i
        ## add intercept
        assay(se)[i+1, ] <- assay(se)[i+1, ] + i
        ## add numeric effect
        assay(se)[i+1, ] <- assay(se)[i + 1, ] + i * se$numeric
        ## add little noise
        assay(se)[i+1, ] <- assay(se)[i+1, ] + rnorm(ncol(se), sd = 0.01)
    }
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric, name = "model2")
    test <- scpDifferentialAnalysis(
        se,
        contrasts = list(
            c("condition", "1", "2"), ## there is change = positive control
            c("condition", "1", "3") ## there is no change = negative control
        ),
        coefficients = "numeric",
        name = "model2"
    )
    expect_equal(
        test$condition_1_vs_2$Estimate,
        structure(0:9, .Names = rownames(se)),
        tolerance = 1E-2
    )
    expect_equal(
        test$condition_1_vs_3$Estimate,
        structure(rep(0, nrow(se)), .Names = rownames(se)),
        tolerance = 1E-2
    )
    expect_equal(
        test$numeric$Estimate,
        structure(0:9, .Names = rownames(se)),
        tolerance = 1E-1
    )
    expect_equal(
        test$condition_1_vs_2$pvalue,
        structure(c(0.8379897, rep(0, 9)), .Names = rownames(se)),
        tolerance = 1E-6
    )
    expect_equal(
        test$condition_1_vs_3$Estimate,
        structure(rep(0, nrow(se)), .Names = rownames(se)),
        tolerance = 1E-2
    )
    expect_equal(
        test$numeric$Estimate,
        structure(0:9, .Names = rownames(se)),
        tolerance = 1E-1
    )
})

test_that(".scpDifferentialAnalysisOnContrast", {
    ## test 1 contrast
    se <- .createMinimalData(nc = 6)
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    estimable <- c("a", "c", "d", "e", "h", "i", "j")
    exp <- DataFrame(
        feature = estimable,
        Estimate = structure((0:9), .Names = rownames(se))[estimable],
        SE = structure(rep(0, 7), .Names = estimable),
        Df = rowSums(!is.na(assay(se)[estimable, ])) - 2,
        tstatistic = structure(c(0.7072246, 3452.7200380, 1206.4382131, 1778.6625013, 1990.7885585, 3534.6922866, 2058.9077196), .Names = estimable),
        pvalue = structure(c(0.5527, rep(0, 6)), .Names = estimable),
        padj = structure(c(0.5527, rep(0, 6)), .Names = estimable)
    )
    metadata(exp)$contrast <- c("condition", "1", "2")
    expect_equal(
        .scpDifferentialAnalysisOnContrast(
            se, contrast = list(c("condition", "1", "2")), name = "model"
        ),
        list(condition_1_vs_2 = exp),
        tolerance = 1E-2
    )
    ## test 2 contrasts
    se <- .createMinimalData(nc = 9)
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    estimable <- c("c", "d", "f", "g", "h", "j")
    c1vs2 <- DataFrame(
        feature = estimable,
        Estimate = structure((0:9), .Names = rownames(se))[estimable],
        SE = structure(rep(0, 6), .Names = estimable),
        Df = rowSums(!is.na(assay(se)[estimable, ])) - 3,
        tstatistic = structure(c(824.8086, 1777.5883, 1975.3086, 1088.6531, 4759.9688, 2614.4243), .Names = estimable),
        pvalue = structure(rep(0, 6), .Names = estimable),
        padj = structure(rep(0, 6), .Names = estimable)
    )
    metadata(c1vs2)$contrast <- c("condition", "1", "2")
    estimable <- c("a", "b", "c", "d", "e", "f", "g", "h", "j")
    c1vs3 <- DataFrame(
        feature = estimable,
        Estimate = structure(rep(0, 9), .Names = estimable),
        SE = structure(rep(0, 9), .Names = estimable),
        Df = sapply(scpModelFitList(se), scpModelFitDf)[estimable],
        tstatistic = structure(c(-0.7072246, 0.3333611, -0.2378374, -0.3530348, 0.3333611, -0.5921999, -0.2115692, -0.7239608, -0.3951838), .Names = estimable),
        pvalue = structure(c(0.5527268, 0.7951513, 0.8513498, 0.7474114, 0.7951513, 0.6137490, 0.8672682, 0.5091646, 0.7308728), .Names = estimable),
        padj = structure(rep(0.8672, 9), .Names = estimable)
    )
    metadata(c1vs3)$contrast <- c("condition", "1", "3")
    expect_equal(
        .scpDifferentialAnalysisOnContrast(
            se,
            contrast = list(
                c("condition", "1", "2"), ## there is change = positive control
                c("condition", "1", "3") ## there is no change = negative control
            ),
            name = "model"
        ),
        list(condition_1_vs_2 = c1vs2,
             condition_1_vs_3 = c1vs3),
        tolerance = 1E-2
    )
})

test_that(".checkContrasts", {
    se <- .createMinimalData(nc = 9)
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    se$numeric <- 1:ncol(se)
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric)
    ## No contrasts = return NULL
    expect_identical(
        .checkContrasts(se, contrasts = NULL, name = "model"),
        NULL
    )
    ## Contrasts is not a list = error
    expect_error(
        .checkContrasts(se, c("condition", "1", "2"), "model"),
        "'contrasts' must be a list."
    )
    ## Some contrast don't have length 3 = error
    expect_error(
        .checkContrasts(se, list(c("condition", "1", "2"), c("condition", "1")), "model"),
        "'contrasts' must be a list with elements of length 3."
    )
    expect_error(
        .checkContrasts(se, list(c("condition", "1", "2"), c("condition", "1", "2", "3")), "model"),
        "'contrasts' must be a list with elements of length 3."
    )
    ## First element is not a modelled variable = error
    expect_error(
        .checkContrasts(se, list(c("foo", "1", "2")), "model"),
        "Effect.s. 'foo' not found."
    )
    ## Variable is not categorical = error
    expect_error(
        .checkContrasts(se, list(c("numeric", "1", "2")), "model"),
        "Provide 'contrasts' only for categorical variables. Problematic variable.s.: numeric."
    )
    ## Missing level in contrast  = error
    expect_error(
        .checkContrasts(se, list(c("condition", "foo", "2")), "model"),
        "Level.s. not found for effect 'condition': foo.",
    )
    expect_error(
        .checkContrasts(se, list(c("condition", "foo", "bar")), "model"),
        "Level.s. not found for effect 'condition': foo, bar.",
    )
    ## Function returns a named contrast list
    expect_identical(
        .checkContrasts(se, list(c("condition", "1", "2"),
                                 c("condition", "2", "1"),
                                 c("condition", "1", "3")),
                        "model"),
        list(condition_1_vs_2 = c("condition", "1", "2"),
             condition_2_vs_1 = c("condition", "2", "1"),
             condition_1_vs_3 = c("condition", "1", "3"))
    )
    ## Duplicate contrasts are removed
    expect_identical(
        .checkContrasts(se, list(c("condition", "1", "2"),
                                 c("condition", "1", "2"),
                                 c("condition", "2", "1")),
                        "model"),
        list(condition_1_vs_2 = c("condition", "1", "2"),
             condition_2_vs_1 = c("condition", "2", "1"))
    )
})

test_that(".contrastToEstimates", {
    ## test contrast with 2-level categorical variable
    se <- .createMinimalData(nc = 6)
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    estimable <- c("a", "c", "d", "e", "h", "i", "j")
    expect_equal(
        .contrastToEstimates(
            se, contrast = c("condition", "1", "2"), name = "model"
        ),
        list(logFc = c(a = 0.000333111231419991, c = 1.99941682634665, d = 2.99825093701587, e = 3.99925006254683, h = 6.99850024999997, i = 7.99800049987503, j = 8.99800037496873),
             se = c(a = 0.000471011958497654, c = 0.000579084549100373, d = 0.00248520886070925, e = 0.0022484591987983, h = 0.003515441265855, i = 0.00226271478572325, j = 0.00437027861387649)),
        tolerance = 1E-2
    )
    ## test contrast with 2-level categorical variable, inverse reference contrast group
    expect_equal(
        .contrastToEstimates(
            se, contrast = c("condition", "2", "1"), name = "model"
        ),
        list(
            logFc = structure(-(0:9), .Names = rownames(se))[estimable],
            se = structure(rep(0, 7), .Names = estimable)
        ),
        tolerance = 1E-2
    )
    ## test contrast with 3-level categorical variable
    se <- .createMinimalData(nc = 9)
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "1", "2"), "model"
        ),
        list(logFc = c(a = NA, b = NA, c = 1.998001998002, d = 2.99733549827914, e = NA, f = 4.99683566486252, g = 5.99483749669705, h = 6.99733442546727, i = NA, j = 8.99550224887556),
             se = c(a = NA, b = NA, c = 0.00242238261001802, d = 0.00168618089019642, e = NA, f = 0.00252964813529043, g = 0.00550665558411383, h = 0.00147003786198391, i = NA, j = 0.00344072014574582)),
        tolerance = 1E-10 ## high precision for first test to assess for reproducibility
    )
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "3", "2"), "model"
        ),
        list(
            logFc = c(a = NA, b = NA, c = 2, d = 3, e = NA, f = 5, g = 6, h = 7, i = 8, j = 9),
            se = c(a = NA, b = NA, c = 0, d = 0, e = NA, f = 0, g = 0, h = 0, i = 0, j = 0)
        ),
        tolerance = 1E-2
    )
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "1", "3"), "model"
        ),
        list(
            logFc = c(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = NA, j = 0),
            se = c(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = NA, j = 0)
        ),
        tolerance = 1E-2
    )
    ## test contrast with 3-level categorical variable, inverse reference contrast group
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "2", "1"), "model"
        ),
        list(
            logFc = -c(a = NA, b = NA, c = 2, d = 3, e = NA, f = 5, g = 6, h = 7, i = NA, j = 9),
            se = c(a = NA, b = NA, c = 0, d = 0, e = NA, f = 0, g = 0, h = 0, i = NA, j = 0)
        ),
        tolerance = 1E-2
    )
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "2", "3"), "model"
        ),
        list(
            logFc = -c(a = NA, b = NA, c = 2, d = 3, e = NA, f = 5, g = 6, h = 7, i = 8, j = 9),
            se = c(a = NA, b = NA, c = 0, d = 0, e = NA, f = 0, g = 0, h = 0, i = 0, j = 0)
        ),
        tolerance = 1E-2
    )
    expect_equal(
        .contrastToEstimates(
            se, c("condition", "3", "1"), "model"
        ),
        list(
            logFc = c(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = NA, j = 0),
            se = c(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = NA, j = 0)
        ),
        tolerance = 1E-2
    )
    ## Test when a level of interest is dropped for a 2-level variable
    se <- .createMinimalData(nr = 2, nc = 10)
    se$condition1 <- as.factor(rep(1:2, length.out = ncol(se)))
    se$condition2 <- as.factor(rep(1:2, each = 5))
    assay(se)[, se$condition1 == 2] <- 10 * assay(se)[, se$condition1 == 2]
    ## Drop one level of condition2 for first row, second row is
    ## positive control.
    assay(se)[1, se$condition2 == 2] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition1 + condition2)
    expect_equal(
        .contrastToEstimates(
            se, contrast = c("condition2", "1", "2"), name = "model"
        ),
        list(logFc = c(a = NA, b = 0.0001874609),
             se = c(a = NA, b  = 0.0005527098)),
        tolerance = 1E-6
    )
    ## Test when a level of interest is dropped for a 3-level variable
    se <- .createMinimalData(nr = 2, nc = 12)
    se$condition1 <- as.factor(rep(1:2, length.out = ncol(se)))
    se$condition2 <- as.factor(rep(1:3, each = 4))
    assay(se)[, se$condition1 == 2] <- 10 * assay(se)[, se$condition1 == 2]
    ## Drop one level of condition2 for first row, second row is
    ## positive control.
    assay(se)[1, se$condition2 == 2] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition1 + condition2)
    expect_equal(
        .contrastToEstimates(
            se, contrast = c("condition2", "1", "2"), name = "model"
        ),
        list(logFc = c(a = NA, b = 0),
             se = c(a = NA, b  = 0.0005126847)),
        tolerance = 1E-6
    )
    ## Make sure the contrast for remaining levels is still estimated
    expect_equal(
        .contrastToEstimates(
            se, contrast = c("condition2", "1", "3"), name = "model"
        ),
        list(logFc = c(a = 0, b = 0),
             se = c(a = 0.0007943138, b = 0.0005127487 )),
        tolerance = 1E-6
    )
})

test_that(".levelsToContrastMatrix", {
    ## 1 level or less = NULL
    expect_identical(
        .levelsToContrastMatrix(contrast = c("condition", "A", "B"),
                                levels = "A"),
        NULL
    )
    expect_identical(
        .levelsToContrastMatrix(contrast = c("condition", "A", "B"),
                                levels = character()),
        NULL
    )
    ## 1 level to test is not in available levels = NA
    expect_identical(
        .levelsToContrastMatrix(contrast = c("condition", "A", "foo"),
                                levels = c("A", "B")),
        NULL
    )
    ## 2-level variable
    expect_identical(
        .levelsToContrastMatrix(c("condition", "A", "B"),
                                levels = c("A", "B")),
        matrix(-2, nrow = 1, dimnames = list(NULL, "condition1"))
    )
    ## 2-level variable, inverse the reference contrast group
    expect_identical(
        .levelsToContrastMatrix(c("condition", "B", "A"),
                                levels = c("A", "B")),
        matrix(2, nrow = 1, dimnames = list(NULL, "condition1"))
    )
    ## multi-level variable
    expect_identical(
        .levelsToContrastMatrix(c("condition", "A", "B"),
                                levels = c("A", "B", "C", "D")),
        matrix(c(-1, 1, 0), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
    expect_identical(
        .levelsToContrastMatrix(c("condition", "B", "C"),
                                levels = c("A", "B", "C", "D")),
        matrix(c(0, -1, 1), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
    ## multi-level variable, include the reference class in contrast
    ## (the last level is encoded as a reference class when coding
    ## categorical variable as a sum contrast)
    expect_identical(
        .levelsToContrastMatrix(c("condition", "B", "D"),
                                levels = c("A", "B", "C", "D")),
        matrix(c(-1, -2, -1), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
    ## multi-level variable, inverse the reference contrast group
    expect_identical(
        .levelsToContrastMatrix(c("condition", "B", "A"),
                                levels = c("A", "B", "C", "D")),
        matrix(-c(-1, 1, 0), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
    expect_identical(
        .levelsToContrastMatrix(c("condition", "C", "B"),
                                levels = c("A", "B", "C", "D")),
        matrix(-c(0, -1, 1), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
    expect_identical(
        .levelsToContrastMatrix(c("condition", "D", "B"),
                                levels = c("A", "B", "C", "D")),
        matrix(-c(-1, -2, -1), nrow = 1, dimnames = list(NULL, paste0("condition", 1:3)))
    )
})

test_that(".scpDifferentialAnalysisOnCoefficient", {
    se <-  SummarizedExperiment(assays = List(
        matrix(0, 10, 100, dimnames = list(paste0("row", 1:10),
                                           paste0("col", 1:100)))
    ))
    set.seed(1234)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/5)] <- NA
    se$numeric <- scale(1:ncol(se)) ## positive control
    se$numeric2 <- runif(ncol(se)) ## negative control
    for (i in 0:9) {
        ## add intercept
        assay(se)[i+1, ] <- assay(se)[i+1, ] + i
        ## add numeric effect
        assay(se)[i+1, ] <- assay(se)[i + 1, ] + i * se$numeric
        ## add little noise
        assay(se)[i+1, ] <- assay(se)[i+1, ] + rnorm(ncol(se), sd = 0.01)
    }
    se <- scpModelWorkflow(se, formula = ~ 1 + numeric + numeric2, name = "model")
    ## Missing coefficient = error
    expect_error(
        .scpDifferentialAnalysisOnCoefficient(
            se, coefficients = "foo",  name = "model"
        ),
        "Some coefficients not found: foo.",
    )
    ## Test 1 coefficient
    exp1 <- DataFrame(
        feature = rownames(se),
        Estimate = c(row1 = 0.000468779211816878, row2 = 0.999256248458486, row3 = 1.98000310009773, row4 = 2.90865795861583, row5 = 4.11359911775159, row6 = 4.90657593661401, row7 = 6.20062598462335, row8 = 7.0746719228876, row9 = 8.10859321582726, row10 = 8.99816838925607),
        SE = c(row1 = 0.00116545745645793, row2 = 0.0011158363022493, row3 = 0.00121486498566538, row4 = 0.001267820701554, row5 = 0.000925222377619779, row6 = 0.00107855077845707, row7 = 0.00106710103503816, row8 = 0.00121263097901738, row9 = 0.00105141615528807, row10 = 0.00099029700050081),
        Df = c(row1 = 78, row2 = 80, row3 = 79, row4 = 75, row5 = 76, row6 = 80, row7 = 80, row8 = 67, row9 = 80, row10 = 75),
        tstatistic = c(row1 = 0.402227648224584, row2 = 895.522261145459, row3 = 1629.81329074464, row4 = 2294.21869752609, row5 = 4446.06531062749, row6 = 4549.2303511506, row7 = 5810.72061691106, row8 = 5834.15073942804, row9 = 7712.06831381209, row10 = 9086.33307452768),
        pvalue = c(row1 = 0.688616634304459, row2 = 8.03183556631787e-162, row3 = 1.40808176700861e-180, row4 = 1.69792709140424e-183, row5 = 1.53036452302809e-207, row6 = 2.73865195054442e-218, row7 = 8.59403606066981e-227, row8 = 6.92241457155388e-193, row9 = 1.25587061385632e-236, row10 = 2.5037770708955e-228),
        padj = c(row1 = 0.688616634304459, row2 = 8.92426174035319e-162, row3 = 1.76010220876076e-180, row4 = 2.42561013057748e-183, row5 = 3.06072904605617e-207, row6 = 6.84662987636104e-218, row7 = 2.86467868688994e-226, row8 = 1.15373576192565e-192, row9 = 1.25587061385632e-235, row10 = 1.25188853544775e-227)
    )
    expect_equal(
        .scpDifferentialAnalysisOnCoefficient(
            se, coefficients = "numeric", name = "model"
        ),
        list(numeric = exp1),
        tolerance = 1E-2
    )
    ## Test 2 coefficients
    expect_equal(
        .scpDifferentialAnalysisOnCoefficient(
            se, coefficients = c("numeric", "numeric2"), name = "model"
        ),
        list(numeric = exp1,
             numeric2 = DataFrame(
                 feature = c("row1", "row2", "row3", "row4", "row5", "row6",
                             "row7", "row8", "row9", "row10"),
                 Estimate = c(row1 = -0.00157904250034698, row2 = 0.000705397503822341, row3 = -0.000455517960428316, row4 = 0.000136301698053502, row5 = -0.00151376184245366, row6 = -0.00102897824142928, row7 = 0.000330019702798889, row8 = 0.0038180174935476, row9 = 0.000345141172191415, row10 = -0.00126466517281232),
                 SE = c(row1 = 0.00116545745645793, row2 = 0.00111583630224931, row3 = 0.00121486498566538, row4 = 0.001267820701554, row5 = 0.000925222377619779, row6 = 0.00107855077845707, row7 = 0.00106710103503816, row8 = 0.00121263097901738, row9 = 0.00105141615528807, row10 = 0.000990297000500811),
                 Df = c(row1 = 78, row2 = 80, row3 = 79, row4 = 75, row5 = 76, row6 = 80, row7 = 80, row8 = 67, row9 = 80, row10 = 75),
                 tstatistic = c(row1 = -1.35486927609183, row2 = 0.632169344553856, row3 = -0.37495356751831, row4 = 0.107508654722576, row5 = -1.63610595578974, row6 = -0.954037827408824, row7 = 0.309267531342136, row8 = 3.14854028934789, row9 = 0.328263143433299, row10 = -1.27705645091599),
                 pvalue = c(row1 = 0.179370436752218, row2 = 0.529078800285441, row3 = 0.708699733514142, row4 = 0.914672538466102, row5 = 0.105952378364963, row6 = 0.342938375045117, row7 = 0.757922079636902, row8 = 0.00245056736923267, row9 = 0.743570356188176, row10 = 0.205521868083594),
                 padj = c(row1 = 0.513804670208984, row2 = 0.842135644041002, row3 = 0.842135644041002, row4 = 0.914672538466102, row5 = 0.513804670208984, row6 = 0.685876750090234, row7 = 0.842135644041002, row8 = 0.0245056736923267, row9 = 0.842135644041002, row10 = 0.513804670208984)
             )),
        tolerance = 1E-10
    )
})


test_that(".coefficientToEstimates", {
    se <- .createMinimalData(nc = 6)
    se$numeric <- scale(1:ncol(se))
    for (i in 0:9) {
        assay(se)[i+1, ] <- assay(se)[i + 1, ] + i * se$numeric
    }
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + numeric)
    expect_equal(
        .coefficientToEstimates(
            se, coefficient = "numeric", name = "model"
        ),
        list(
            beta = c(a = 0, c = 2.21626476138346, d = 3.33642362977023, e = 3.26435414663759, h = 7.78498846946388, i = 9.23452613165679, j = 4.80829820532098),
            se = c(a = 0.000204005125241858, c = 0.000339941667740344, d = 0.00167301187722831, e = 0.0019871773566072, h = 0.00394651632004783, i = 0.00218822549009997, j = 0.00412383501963632)),
        tolerance = 1E-10
    )
})

test_that(".computeTTest", {
    ## When < 3000 features, use BH for multiple testing
    n <- 300
    beta <- structure(seq(0, 10, length.out = n), .Names = 1:n)
    se <- seq(0.1, 1E-2, length.out = n)
    df <- rep(20, n)
    weights <- rnorm(n)
    expTstat <- beta / se
    expPval <- 2 * pt(abs(expTstat), df, lower.tail = FALSE)
    expPadj <- p.adjust(expPval, method = "BH")
    expect_identical(
        .computeTTest(beta = beta, se = se, df = df, weights = weights),
        DataFrame(
            feature = as.character(1:n),
            Estimate = beta, SE = se, Df = df,
            tstatistic = expTstat, pvalue = expPval,
            padj = structure(expPadj, .Names = 1:n)
        )
    )
    ## When >= 3000 features, use IHW for multiple testing
    n <- 3000
    beta <- structure(seq(0, 10, length.out = n), .Names = 1:n)
    se <- seq(0.1, 1E-2, length.out = n)
    df <- rep(20, n)
    set.seed(123)
    weights <- rnorm(n)
    expTstat <- beta / se
    expPval <- 2 * pt(abs(expTstat), df, lower.tail = FALSE)
    expPadj <- adj_pvalues(ihw(pvalues = expPval, covariates = weights,
                               alpha = 0.05, adjustment_type = "BH"))
    expect_identical(
        .computeTTest(beta = beta, se = se, df = df, weights = weights),
        DataFrame(
            feature = as.character(1:n),
            Estimate = beta, SE = se, Df = df,
            tstatistic = expTstat, pvalue = expPval,
            padj = structure(expPadj, .Names = 1:n)
        )
    )
})

test_that(".adjustPvalue", {
    ## No weights = BH
    n <- 100
    expect_identical(
        .adjustPvalue(pval = seq(0, 1, length.out = n), weights = NULL),
        p.adjust(seq(0, 1, length.out = n), method = "BH")
    )
    ## Less then 3000 features = BH
    weights <- runif(n)
    expect_identical(
        .adjustPvalue(pval = seq(0, 1, length.out = n), weights = NULL),
        p.adjust(seq(0, 1, length.out = n), method = "BH")
    )
    ## Else = IHW
    ## example inspired from ihw's man page
    n <- 3000
    set.seed(123)
    weights   <- runif(n, min=0, max=2.5)
    h   <- rbinom(n, 1, 0.1) # hypothesis true or false
    z   <- rnorm(n, h * weights) # Z-score
    pval <- 1-pnorm(z)
    expect_identical(
        .adjustPvalue(pval = pval, weights = weights),
        adj_pvalues(ihw(
            pvalues = pval, covariates = weights, alpha = 0.05
        ))
    )
})

test_that("scpDifferentialAggregate", {
    se <- .createMinimalData(nc = 6)
    set.seed(1234)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/5)] <- NA
    se$condition <- as.factor(rep(1:3, length.out = ncol(se)))
    se$numeric <- scale(1:ncol(se))
    for (i in 0:9) {
        assay(se)[i+1, se$condition == 2] <- assay(se)[i + 1, se$condition == 2] + i
        assay(se)[i+1, ] <- assay(se)[i+1, ] + i
        assay(se)[i+1, ] <- assay(se)[i + 1, ] + i * se$numeric
        assay(se)[i+1, ] <- assay(se)[i+1, ] + rnorm(ncol(se), sd = 0.01)
    }
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric)
    daRes <- scpDifferentialAnalysis(
        se, contrasts = list(c("condition", "1", "2")),
        coefficients = "numeric"
    )
    ## fcol not in table = error
    expect_error(
        scpDifferentialAggregate(differentialList = daRes, fcol = "prot"),
        "fcol %in% colnames.*is not TRUE"
    )
    ## Test scpDifferentialAggregate
    daRes <- endoapply(daRes, function(x) {
        x$prot <- c(rep(1:2, each = 4), c(NA, NA)) ## NAs rows will be removed
        x$removed <- 1:10 ## annotation columns not redundant within fcol are removed
        x$keep <- c(rep(c("a", "b"), each = 4), c("c", "c")) ## annotation columns redundant within fcol are kept
        x
    })
    library("metapod")
    metapodRes1 <- combineGroupedPValues(
        daRes[[1]]$pvalue[-(9:10)], daRes[[1]][["prot"]][-(9:10)]
    )
    metapodRes2 <- combineGroupedPValues(
        daRes[[2]]$pvalue[-(9:10)], daRes[[2]][["prot"]][-(9:10)]
    )
    expect_identical(
        scpDifferentialAggregate(daRes, "prot"),
        DataFrameList(
            condition_1_vs_2 = DataFrame(
                feature = 1:2,
                Estimate = daRes[[1]]$Estimate[metapodRes1$representative],
                pvalue = metapodRes1$p.value,
                padj = p.adjust(metapodRes1$p.value, method = "BH"),
                prot = 1:2,
                keep = c("a", "b"),
                .n = c(4L, 4L)
            ),
            numeric = DataFrame(
                feature = 1:2,
                Estimate = daRes[[2]]$Estimate[metapodRes2$representative],
                pvalue = metapodRes2$p.value,
                padj = p.adjust(metapodRes2$p.value, method = "BH"),
                prot = 1:2,
                keep = c("a", "b"),
                .n = c(4L, 4L)
            )
        )
    )
    ## Change metapod defaults
    metapodRes1 <- combineGroupedPValues(
        daRes[[1]]$pvalue[-(9:10)], daRes[[1]][["prot"]][-(9:10)],
        method = "berger", weights = runif(10)
    )
    metapodRes2 <- combineGroupedPValues(
        daRes[[2]]$pvalue[-(9:10)], daRes[[2]][["prot"]][-(9:10)],
        method = "berger", weights = runif(10)
    )
    expect_identical(
        scpDifferentialAggregate(daRes, "prot", method = "berger", weights = runif(10)),
        DataFrameList(
            condition_1_vs_2 = DataFrame(
                feature = 1:2,
                Estimate = daRes[[1]]$Estimate[metapodRes1$representative],
                pvalue = metapodRes1$p.value,
                padj = p.adjust(metapodRes1$p.value, method = "BH"),
                prot = 1:2,
                keep = c("a", "b"),
                .n = c(4L, 4L)
            ),
            numeric = DataFrame(
                feature = 1:2,
                Estimate = daRes[[2]]$Estimate[metapodRes2$representative],
                pvalue = metapodRes2$p.value,
                padj = p.adjust(metapodRes2$p.value, method = "BH"),
                prot = 1:2,
                keep = c("a", "b"),
                .n = c(4L, 4L)
            )
        )
    )
})

test_that("scpVolcanoPlot", {
    se <-  SummarizedExperiment(assays = List(
        matrix(0, 100, 100, dimnames = list(paste0("row", 1:100),
                                           paste0("col", 1:100)))
    ))
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    assay(se)[, se$condition == 2] <- 1:nrow(se) - nrow(se) / 2 + assay(se)[, se$condition == 2]
    set.seed(124)
    assay(se) <- assay(se) + rnorm(length(assay(se)))
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition)
    daRes <- scpDifferentialAnalysis(
        se, contrast = list(c("condition", "1", "2"))
    )
    ## default plot
    set.seed(124) ## ggrepel is stochastic
    expect_doppelganger(
        "scpVolcanoPlot default",
        scpVolcanoPlot(daRes)

    )
    ## change FDR line
    expect_doppelganger(
        "scpVolcanoPlot change fdrLine",
        scpVolcanoPlot(daRes, fdrLine = 1E-5)

    )
    ## change number of labels
    expect_doppelganger(
        "scpVolcanoPlot change top",
        scpVolcanoPlot(daRes, top = 30, labelParams = list(max.overlaps = 100))

    )
    expect_doppelganger(
        "scpVolcanoPlot change top is zero",
        scpVolcanoPlot(daRes, top = 0)
    )
    ## change label filter
    ## label filter is absent = error
    expect_error(
        scpVolcanoPlot(daRes, by = "foo"),
        "'foo' not found in differentialList tables."
    )
    expect_doppelganger(
        "scpVolcanoPlot change by",
        scpVolcanoPlot(daRes, by = "Estimate")
    )
    ## change ordering direction
    expect_doppelganger(
        "scpVolcanoPlot decreasing",
        scpVolcanoPlot(daRes, by = "Estimate", decreasing = TRUE)
    )
    ## change labelling variable
    ## labelling variable is absent = error
    expect_error(
        scpVolcanoPlot(daRes, textBy = "foo"),
        "'foo' not found in results. Use scpAnnotateResults.. to add custom annotations."
    )
    expect_doppelganger(
        "scpVolcanoPlot textBy",
        scpVolcanoPlot(daRes, textBy = "Df")
    )
    ## change point params
    expect_doppelganger(
        "scpVolcanoPlot pointParams",
        scpVolcanoPlot(daRes, pointParams = list(aes(size = Df)))
    )
    ## change label params
    expect_doppelganger(
        "scpVolcanoPlot labelParams",
        scpVolcanoPlot(daRes, labelParams = list(aes(colour = Df)))
    )
    expect_doppelganger(
        "scpVolcanoPlot labelParams change label",
        scpVolcanoPlot(daRes, labelParams = list(aes(label = Df)))
    )
})

test_that(".filterDifferentialData", {
    x <- data.frame(criteria = 1:100)
    ## By is absent = error
    expect_error(
        .filterDifferentialData(df = x, by = "foo"),
        "'foo' not found in differentialList tables."
    )
    ## top is zero = return empty table
    expect_identical(
        .filterDifferentialData(df = x, top = 0, by = "criteria"),
        data.frame(criteria = integer())
    )
    ## top is greater than nrows table = return table without filtering
    expect_identical(
        .filterDifferentialData(df = x, top = 1000, by = "criteria"),
        x
    )
    ## test by
    x <- data.frame(criteria1 = 1:100, criteria2 = 100:1)
    expect_identical(
        .filterDifferentialData(df = x, top = 10, by = "criteria1",
                                decreasing = FALSE),
        x[1:10, ]
    )
    expect_identical(
        .filterDifferentialData(df = x, top = 10, by = "criteria2",
                                decreasing = FALSE),
        x[100:91, ]
    )
    ## Test decreasing
    expect_identical(
        .filterDifferentialData(df = x, top = 10, by = "criteria1",
                                decreasing = TRUE),
        x[100:91, ]
    )
    expect_identical(
        .filterDifferentialData(df = x, top = 10, by = "criteria2",
                                decreasing = TRUE),
        x[1:10, ]
    )
})

test_that(".plotVolcano", {
    x <- data.frame(
        Estimate = -10:10,
        padj = runif(21)/10,
        names = paste("feat", 1:21)
    )
    set.seed(124) ## ggrepel is stochastic
    ## Default
    expect_doppelganger(
        ".plotVolcano default",
        .plotVolcano(x, pointParams = list(), labelParams = list(),
                     textBy = "names")
    )
    ## Change FDR line
    expect_doppelganger(
        ".plotVolcano change fdr",
        .plotVolcano(x, pointParams = list(), labelParams = list(),
                     textBy = "names", fdrLine = 0.01)
    )
    ## Change contrast
    expect_doppelganger(
        ".plotVolcano change contrast",
        .plotVolcano(x, pointParams = list(), labelParams = list(),
                     textBy = "names", contrast = c("condition", "A", "B"))
    )
    ## Change point and label params
    expect_doppelganger(
        ".plotVolcano change aes params",
        .plotVolcano(x, pointParams = list(aes(col = padj), size = 5),
                     labelParams = list(aes(size = -padj), colour = "red"),
                     textBy = "names", contrast = c("condition", "A", "B"))
    )
})

test_that(".annotateDirection", {
    ## Annotate axis with 2 directions, that is logFoldChange contains
    ## positive and negative values
    expect_identical(
        .annotateDirection(-10:10, c("Group", "foo", "bar")),
        xlab("foo <-   log2(Fold change)   -> bar")
    )
    ## Reverse direction
    expect_identical(
        .annotateDirection(-10:10, c("Group", "bar", "foo")),
        xlab("bar <-   log2(Fold change)   -> foo")
    )
    ## Only positive direction = right annotation
    expect_identical(
        .annotateDirection(10:0, c("Group", "foo", "bar")),
        xlab("log2(Fold change)   -> bar")
    )
    ## Only negative direction = left annotation
    expect_identical(
        .annotateDirection(-10:0, c("Group", "foo", "bar")),
        xlab("foo <-   log2(Fold change)")
    )
    ## Only zero = no direction annotation
    expect_identical(
        .annotateDirection(0, c("Group", "foo", "bar")),
        xlab("log2(Fold change)")
    )
})
