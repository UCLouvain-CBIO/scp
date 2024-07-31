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

test_that("scpVarianceAnalysis", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    estimable <- c("a", "c", "e", "h", "i")
    ## Standard use case
    expect_equal(
        scpVarianceAnalysis(se, "model"),
        List(
            Residuals = DataFrame(
                feature = estimable,
                SS = structure(rep(0, 5), .Names = estimable),
                df = structure(rep(1, 5), .Names = estimable),
                percentExplainedVar = structure(rep(0, 5), .Names = estimable)
            ),
            condition = DataFrame(
                feature = estimable,
                SS = structure(c(1, 0.997, 0.748, 0.187, 0.748), .Names = estimable),
                df = structure(rep(1, 5), .Names = estimable),
                percentExplainedVar = structure(c(100, 100, 100, 50, 100), .Names = estimable)
            ),
            batch = DataFrame(
                feature = estimable,
                SS = structure(c(0, 0, 0, 0.1872, 0), .Names = estimable),
                df = structure(rep(1, 5), .Names = estimable),
                percentExplainedVar = structure(c(0, 0, 0, 50, 0), .Names = estimable)
            )
        ),
        tolerance = 1E-3
    )
    ## Also works without providing name
    expect_identical(
        scpVarianceAnalysis(se, "model"),
        scpVarianceAnalysis(se)
    )
})

test_that(".computeSumsOfSquares", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    estimable <- c("a", "c", "e", "h", "i")
    ## Standard use case
    expect_equal(
        .computeSumsOfSquares(se, "model"),
        matrix(
            c(1, 0.9986, 0.7493, 0.7494, 0.7491,
              rep(0, 5),
              1, 0.997, 0.748, 0.187, 0.748,
              0, 0, 0, 0.1872, 0),
            ncol = 4,
            dimnames = list(estimable, c("SStotal", "Residuals", "condition", "batch"))
        ),
        tolerance = 1E-3
    )
    expect_equal(
        .computeSumsOfSquares(se, "model"),
        .computeSumsOfSquares(se)
    )
})

test_that(".computeTotalSS", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    estimable <- c("a", "c", "e", "h", "i")
    ## Standard use case
    ## Note the intercept = grand mean = 1.5
    y <- assay(se)[estimable, ]
    exp <- rowSums((y - 1.5)^2, na.rm = TRUE)
    expect_equal(
        .computeTotalSS(se),
        exp,
        tolerance = 1E-3
    )
})

test_that(".computeResidualSS", {
    se <- .createMinimalData()
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + batch)
    estimable <- c("a", "c", "e", "h", "i")
    ## Standard use case
    expect_equal(
        .computeResidualSS(se),
        structure(rep(0, 5), .Names = estimable),
        tolerance = 1E-5
    )
})

test_that(".computeEffectSS", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + batch + condition)
    estimable <- c("a", "c", "e", "h", "i")
    ## Standard use case
    expect_equal(
        .computeEffectSS(se),
        matrix(
            c(0, 0, 0, 0.1872, 0,
              1, 0.9975, 0.7488, 0.1872, 0.7481),
            ncol = 2,
            dimnames = list(estimable, c("batch", "condition"))
        ),
        tolerance = 1E-3
    )
})

test_that(".explainedVarianceDenominator", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    ss <- .computeSumsOfSquares(se)
    ## Used the total SS
    expect_identical(
        .explainedVarianceDenominator(ss, useTotalSS = TRUE),
        ss[, "SStotal"]
    )
    ## Don't use the total SS but sum all SS
    expect_identical(
        .explainedVarianceDenominator(ss, useTotalSS = FALSE),
        rowSums(ss[, -1])
    )
})

test_that("scpVarianceAggregate", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    ## Standard use case
    varRes <- scpVarianceAnalysis(se, "model")
    varRes <- endoapply(varRes, function(x) {
        x$prot <- c("b", "b", "b", "a", "a")
        x
    })
    expect_equal(
        scpVarianceAggregate(varRes, fcol = "prot"),
        List(
            Residuals = DataFrame(
                feature = c("a", "b"),
                SS = structure(
                    c(sum(varRes$Residuals$SS[4:5]), sum(varRes$Residuals$SS[1:3])),
                    .Names = c("a", "b")
                ),
                df = c(mean(varRes$Residuals$df[4:5]), mean(varRes$Residuals$df[1:3])),
                percentExplainedVar = structure(
                    c(mean(varRes$Residuals$percentExplainedVar[4:5]),
                      mean(varRes$Residuals$percentExplainedVar[1:3])),
                    .Names = c("a", "b")
                ),
                prot = c("a", "b"),
                .n = c(2, 3),
                row.names = c("a", "b")
            ),
            condition = DataFrame(
                feature = c("a", "b"),
                SS = structure(
                    c(sum(varRes$condition$SS[4:5]), sum(varRes$condition$SS[1:3])),
                    .Names = c("a", "b")
                ),
                df = c(mean(varRes$condition$df[4:5]), mean(varRes$condition$df[1:3])),
                percentExplainedVar = structure(
                    c(mean(varRes$condition$percentExplainedVar[4:5]),
                      mean(varRes$condition$percentExplainedVar[1:3])),
                    .Names = c("a", "b")
                ),
                prot = c("a", "b"),
                .n = c(2, 3),
                row.names = c("a", "b")
            ),
            batch = DataFrame(
                feature = c("a", "b"),
                SS = structure(
                    c(sum(varRes$batch$SS[4:5]), sum(varRes$batch$SS[1:3])),
                    .Names = c("a", "b")
                ),
                df = c(mean(varRes$batch$df[4:5]), mean(varRes$batch$df[1:3])),
                percentExplainedVar = structure(
                    c(mean(varRes$batch$percentExplainedVar[4:5]),
                      mean(varRes$batch$percentExplainedVar[1:3])),
                    .Names = c("a", "b")
                ),
                prot = c("a", "b"),
                .n = c(2, 3),
                row.names = c("a", "b")
            )
        ),
        tolerance = 1E-7
    )
})

test_that("scpVariancePlot", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    varRes <- endoapply(varRes, function(x) {
        x$prot <- c("b", "b", "b", "a", "a")
        x
    })
    ## Combined default
    vdiffr::expect_doppelganger(
        "scpVariancePlot default",
        scpVariancePlot(varRes)
    )
    ## Combined, apply filtering
    vdiffr::expect_doppelganger(
        "scpVariancePlot combined change filtering",
        scpVariancePlot(varRes, effect = "batch", by = "SS", top = 2,
                        decreasing = FALSE)
    )
    ## Combined, fcol has no effect
    expect_equal(
        scpVariancePlot(varRes, effect = "batch", by = "SS", top = 2,
                        decreasing = FALSE),
        scpVariancePlot(varRes, effect = "batch", by = "SS", top = 2,
                        decreasing = FALSE, fcol = "prot"),
        tolerance = 1E-10
    )
    ## Individual default
    vdiffr::expect_doppelganger(
        "scpVariancePlot individual default",
        scpVariancePlot(varRes, combined = FALSE)
    )
    ## Individual, apply filtering
    vdiffr::expect_doppelganger(
        "scpVariancePlot individual change filtering",
        scpVariancePlot(varRes, combined = FALSE, effect = "condition",
                        by = "SS", top = 4, decreasing = FALSE)
    )
    ## Individual, apply filtering and fcol
    vdiffr::expect_doppelganger(
        "scpVariancePlot combined with fcol",
        scpVariancePlot(varRes, combined = FALSE, effect = "condition",
                        by = "SS", top = 4, decreasing = FALSE, fcol = "prot")
    )
    ## Change colour seed
    vdiffr::expect_doppelganger(
        "scpVariancePlot colour seed",
        scpVariancePlot(varRes, colourSeed = 1)
    )
})

test_that(".gatherVarianceData", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    expect_identical(
        .gatherVarianceData(varRes),
        rbind(
            cbind(varRes$Residuals, effectName = "Residuals"),
            cbind(varRes$condition, effectName = "condition"),
            cbind(varRes$batch, effectName = "batch")
        )
    )
})

test_that(".filterVarianceData", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    df <- .gatherVarianceData(varRes)
    ## Test a first set of filtering parameters
    topFeatures <- c("i", "c", "e")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(varianceTable = df, top = 3, by = "SS",
                            effect = "Residuals", decreasing = TRUE),
        exp
    )
    ## Test top, top < number of features
    topFeatures <- c("i", "c", "e", "h")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(df, 4, "SS", "Residuals", TRUE),
        exp
    )
    ## Test top, top > number of features
    topFeatures <- c("i", "c", "e", "h", "a")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(df, 10, "SS", "Residuals", TRUE),
        exp
    )
    ## Test by
    topFeatures <- c("h", "i", "e")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(df, 3, "percentExplainedVar", "Residuals", TRUE),
        exp
    )
    ## Test effect
    topFeatures <- c("h", "i", "c")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(df, 3, "SS", "batch", TRUE),
        exp
    )
    ## Test decreasing
    topFeatures <- c("a", "h", "e")
    exp <- df[df$feature %in% topFeatures, ]
    exp$feature <- factor(exp$feature, levels = topFeatures)
    expect_identical(
        .filterVarianceData(df, 3, "SS", "Residuals", FALSE),
        exp
    )
})

test_that(".plotExplainedVarianceByFeature", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    df <- .gatherVarianceData(varRes)
    df$prot <- rep(c("b", "b", "b", "a", "a"), 3)
    ## Without fcol
    vdiffr::expect_doppelganger(
        ".plotExplainedVarianceByFeature without fcol",
        .plotExplainedVarianceByFeature(varianceTable = df)
    )
    ## With fcol
    vdiffr::expect_doppelganger(
        ".plotExplainedVarianceByFeature with fcol",
        .plotExplainedVarianceByFeature(df, fcol = "prot")
    )
})

test_that(".plotExplainedVarianceCombined", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    df <- .gatherVarianceData(varRes)
    vdiffr::expect_doppelganger(
        ".plotExplainedVarianceCombined",
        .plotExplainedVarianceCombined(df)
    )
})

test_that(".prettifyVariancePlot", {
    se <- .createMinimalData()
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$batch <- as.factor(c(1, 1, 2, 2, 2))
    assay(se)[, se$condition == 2] <- 2
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + batch)
    varRes <- scpVarianceAnalysis(se)
    df <- .gatherVarianceData(varRes)
    ## Prettify .plotExplainedVarianceByFeature
    pl <- .plotExplainedVarianceByFeature(df)
    set.seed(123)
    vdiffr::expect_doppelganger(
        ".prettifyVariancePlot by feature",
        .prettifyVariancePlot(pl)
    )
    ## Prettify .plotExplainedVarianceCombined
    pl <- .plotExplainedVarianceCombined(df)
    set.seed(123)
    vdiffr::expect_doppelganger(
        ".prettifyVariancePlot combined",
        .prettifyVariancePlot(pl)
    )
    ## Test with > 8 effects (only 8 colors available)
    require("SummarizedExperiment")
    a <- matrix(rnorm(1000), 10, 100)
    rownames(a) <- letters[1:10]
    colnames(a) <- paste("col", 1:100)
    se <- SummarizedExperiment(assays = List(assay = a))
    for (i in 1:10) {## Create 10 mock variables = 10 effects
        colData(se)[[paste0("mock", i)]] <- c(i:10, seq_len(i-1))
    }
    se <- scpModelWorkflow(
        se, formula = ~mock1 + mock2 + mock3 + mock4 + mock5 + mock6 + mock7 + mock8 + mock9
        )
    varRes <- scpVarianceAnalysis(se)
    df <- .gatherVarianceData(varRes)
    pl <- .plotExplainedVarianceCombined(df)
    set.seed(123)
    vdiffr::expect_doppelganger(
        ".prettifyVariancePlot with 9 colors",
        .prettifyVariancePlot(pl)
    )
})