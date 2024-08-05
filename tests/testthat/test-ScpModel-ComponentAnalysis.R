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

## ---- SCP Component Analysis ----

test_that("scpComponentAnalysis", {
    require("SummarizedExperiment")
    m <- matrix(
        1, 10, 15, dimnames = list(paste0("row", 1:10), paste0("col", 1:15))
    )
    se <- SummarizedExperiment(assays = List(assay = m))
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$numeric <- as.vector(scale(1:ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    assay(se) <- sweep(assay(se), 2, se$numeric * 2, "+")
    set.seed(124)
    assays(se)[[2]] <- assays(se)[[1]]
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric,
                           name = "withNA", i = 1)
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric,
                           name = "noNA", i = 2)
    ## effects not modelled = error
    expect_error(
        scpComponentAnalysis(se, effect = "foo"),
        "'foo' is/are not modelled effects."
    )
    expect_error(
        scpComponentAnalysis(se, effect = c("foo1", "foo2")),
        "'foo1', 'foo2' is/are not modelled effects."
    )
    ## method not known = error
    expect_error(
        scpComponentAnalysis(se, method = "foo"),
        "Allowed values for 'method' are 'APCA', 'ASCA', 'ASCA.E'."
    )
    ## No method, no residuals, no unmodelled = message
    expect_message(expect_identical(
        scpComponentAnalysis(se, unmodelled = FALSE, residuals = FALSE),
        List()
    ), "No components were computed. Make sure to provide 'method' or set 'unmodelled = TRUE' or 'residuals = TRUE'.")
    ## svd (model with no missing values) on residuals and unmodelled data
    expect_identical(
        test <- scpComponentAnalysis(se, name = "noNA"), ## auto = svd when no NA
        scpComponentAnalysis(se, pcaFUN = "svd", name = "noNA")
    )
    res <- .componentsToTable(.svdWrapper(scpModelResiduals(se, name = "noNA")))
    unmod <- .componentsToTable(.svdWrapper(assay(se, 2)))
    expect_identical(
        test,
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample

            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature
            )
        )
    )
    ## nipals (model with missing values, default) on residuals and unmodelled data
    expect_identical(
        test <- scpComponentAnalysis(se), ## auto = svd when no NA
        scpComponentAnalysis(se, pcaFUN = "nipals")
    )
    res <- .componentsToTable(.nipalsWrapper(scpModelResiduals(se)))
    unmod <- .componentsToTable(.nipalsWrapper(assay(se, 1)))
    expect_identical(
        test,
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample

            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature
            )
        )
    )
    ## APCA using nipals for one effect
    apcaCondition <- .componentsToTable(.nipalsWrapper(
        scpModelEffects(se)[["condition"]] + scpModelResiduals(se)
    ))
    expect_identical(
        scpComponentAnalysis(se, method = "APCA", effect = "condition"),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample,
                APCA_condition = apcaCondition$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature,
                APCA_condition = apcaCondition$byFeature
            )
        )
    )
    ## APCA using nipals (model with no missing values) for two effects
    apcaNumeric <- .componentsToTable(.nipalsWrapper(
        scpModelEffects(se)[["numeric"]] + scpModelResiduals(se)
    ))
    expect_identical(
        scpComponentAnalysis(se, method = "APCA", effect = c("condition", "numeric")),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample,
                APCA_condition = apcaCondition$bySample,
                APCA_numeric = apcaNumeric$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature,
                APCA_condition = apcaCondition$byFeature,
                APCA_numeric = apcaNumeric$byFeature
            )
        )
    )
    ## Missing effect = all effects
    expect_identical(
        scpComponentAnalysis(se, method = "APCA", effect = c("condition", "numeric")),
        scpComponentAnalysis(se, method = "APCA")
    )
    ## ASCA
    ascaCondition <- .componentsToTable(.nipalsWrapper(
        scpModelEffects(se)[["condition"]]
    ))
    expect_identical(
        scpComponentAnalysis(se, method = "ASCA", effect = "condition"),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample,
                ASCA_condition = ascaCondition$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature,
                ASCA_condition = ascaCondition$byFeature
            )
        )
    )
    ## ASCA-E
    ascaeCondition <- .componentsToTable(.runASCA.E(
        scpModelEffects(se)[["condition"]], scpModelResiduals(se), .nipalsWrapper
    ))
    expect_identical(
        scpComponentAnalysis(se, method = "ASCA.E", effect = "condition"),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample,
                ASCA.E_condition = ascaeCondition$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature,
                ASCA.E_condition = ascaeCondition$byFeature
            )
        )
    )
    ## All methods
    expect_identical(
        scpComponentAnalysis(se, method = c("APCA", "ASCA", "ASCA.E"),
                             effect = "condition"),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                residuals = res$bySample,
                APCA_condition = apcaCondition$bySample,
                ASCA_condition = ascaCondition$bySample,
                ASCA.E_condition = ascaeCondition$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                residuals = res$byFeature,
                APCA_condition = apcaCondition$byFeature,
                ASCA_condition = ascaCondition$byFeature,
                ASCA.E_condition = ascaeCondition$byFeature
            )
        )
    )
    ## no residuals
    expect_identical(
        scpComponentAnalysis(se, method = "ASCA", effect = "condition",
                             residuals = FALSE),
        List(
            bySample = SimpleList(
                unmodelled = unmod$bySample,
                ASCA_condition = ascaCondition$bySample
            ),
            byFeature = SimpleList(
                unmodelled = unmod$byFeature,
                ASCA_condition = ascaCondition$byFeature
            )
        )
    )
    ## no unmodelled
    expect_identical(
        scpComponentAnalysis(se, method = "ASCA", effect = "condition",
                             unmodelled = FALSE),
        List(
            bySample = SimpleList(
                residuals = res$bySample,
                ASCA_condition = ascaCondition$bySample
            ),
            byFeature = SimpleList(
                residuals = res$byFeature,
                ASCA_condition = ascaCondition$byFeature
            )
        )
    )
})

test_that(".getPcaFunction", {
    ## "nipals" returns .nipalsWrapper
    expect_identical(
        .getPcaFunction(pcaFun = "nipals"),
        .nipalsWrapper
    )
    ## "svd" returns .svdWrapper
    expect_identical(
        .getPcaFunction(pcaFun = "svd"),
        .svdWrapper
    )
    ## "auto" without missing values returns .svdWrapper
    se <- .createMinimalData()
    se$foo <- runif(ncol(se))
    se <- scpModelWorkflow(se, ~ 1 + foo)
    expect_identical(
        .getPcaFunction(pcaFun = "auto", object = se, name = "model"),
        .svdWrapper
    )
    ## "auto" with missing values returns .nipalsWrapper
    se <- .createMinimalData()
    assay(se)[1] <- NA
    se$foo <- runif(ncol(se))
    se <- scpModelWorkflow(se, ~ 1 + foo)
    expect_identical(
        .getPcaFunction(pcaFun = "auto", se, "model"),
        .nipalsWrapper
    )
    ## other = error
    expect_error(
        .getPcaFunction("foo"),
        "Available PCA functions are: 'nipals' or 'svd'"
    )
})

test_that(".addPcaToList", {
    set.seed(1234)
    scores <- matrix(
        runif(20), ncol = 2,
        dimnames = list(paste0("col", 1:10), paste0("PC", 1:2))
    )
    eigenvectors <- matrix(
        runif(30), ncol = 2,
        dimnames = list(paste0("row", 1:15), paste0("PC", 1:2))
    )
    eigenvalues <- structure(runif(2), .Names = paste0("PC", 1:2))
    proportionVariance <- structure(runif(2), .Names = paste0("PC", 1:2))
    bySample <- DataFrame(
        scores, cell = paste0("col", 1:10)
    )
    byFeature <- DataFrame(
        eigenvectors, feature = paste0("row", 1:15)
    )
    metadata(bySample)$proportionVariance <-
        metadata(byFeature)$proportionVariance <-
        structure(proportionVariance, .Names = paste0("PC", 1:2))
    exp <- SimpleList(bySample = bySample, byFeature = byFeature)
    ## Append to empty list
    expect_identical(
        test <- .addPcaToList(
            pcaResults = list(
                scores = scores, eigenvectors = eigenvectors,
                eigenvalues = eigenvalues, proportionVariance = proportionVariance
            ),
            pcaList = List()
        ),
        exp
    )
    ## Test prefix and add to non empty list
    expect_identical(
        test <- .addPcaToList(pcaResults = list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues, proportionVariance = proportionVariance
        ),
        pcaList = test, prefix = "new_"),
        exp <- c(exp, List(new_bySample = bySample, new_byFeature = byFeature))
    )
})

test_that(".componentsToTable", {
    scores <- matrix(
        runif(20), ncol = 2,
        dimnames = list(paste0("col", 1:10), paste0("PC", 1:2))
    )
    eigenvectors <- matrix(
        runif(30), ncol = 2,
        dimnames = list(paste0("row", 1:15), paste0("PC", 1:2))
    )
    eigenvalues <- structure(runif(2), .Names = paste0("PC", 1:2))
    proportionVariance <- structure(runif(2), .Names = paste0("PC", 1:2))
    bySample <- DataFrame(
        scores, cell = paste0("col", 1:10)
    )
    byFeature <- DataFrame(
        eigenvectors, feature = paste0("row", 1:15)
    )
    metadata(bySample)$proportionVariance <-
        metadata(byFeature)$proportionVariance <-
        structure(proportionVariance, .Names = paste0("PC", 1:2))
    exp <- DataFrameList(bySample = bySample, byFeature = byFeature)
    ## Append to empty list
    expect_identical(
        .componentsToTable(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues, proportionVariance = proportionVariance
        )),
        exp
    )
})

test_that(".checkPcaResults", {
    ## x must be a list
    expect_error(
        .checkPcaResults(x = 1:10),
        "PCA results must be a list."
    )
    ## x must a named list with expected elements
    expect_error(
        .checkPcaResults(list()),
        "PCA results must at least contain the following elements: scores, eigenvectors, eigenvalues, proportionVariance."
    )
    expect_error(
        .checkPcaResults(list(foo = 1:10, bar = 1:10)),
        "PCA results must at least contain the following elements: scores, eigenvectors, eigenvalues, proportionVariance."
    )
    expect_error(
        .checkPcaResults(list(scores = 1:10, eigenvectors = 1:10)),
        "PCA results must at least contain the following elements: scores, eigenvectors, eigenvalues, proportionVariance."
    )
    ## Scores, eigenvectors, eigenvalues and proportion variance have
    ## no names = error
    scores <- matrix(runif(20), ncol = 2)
    eigenvectors <- matrix(runif(30), ncol = 2)
    eigenvalues <- runif(2)
    proportionVariance <- runif(2)
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues, proportionVariance = proportionVariance
        )),
        "Components must be named."
    )
    ## Names do not match = error
    ## No match because no name in score, eigenvectors or prop var
    names(eigenvalues) <- paste0("PC", 1:2)
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance
        )),
        "Components in scores do not match the eigenvalues."
    )
    colnames(scores) <- paste0("PC", 1:2)
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance
        )),
        "Components in eigenvectors do not match the eigenvalues."
    )
    colnames(eigenvectors) <- paste0("PC", 1:2)
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance
        )),
        "The proportions of variance explained do not match the eigenvalues."
    )
    names(proportionVariance) <- paste0("PC", 1:2)
    ## Wrong dimensions = error
    expect_error(
        .checkPcaResults(list(
            scores = scores[, c(1, 1, 2)], eigenvectors = eigenvectors,
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance
        )),
        "Components in scores do not match the eigenvalues."
    )
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors[, c(1, 1, 2)],
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance
        )),
        "Components in eigenvectors do not match the eigenvalues."
    )
    expect_error(
        .checkPcaResults(list(
            scores = scores, eigenvectors = eigenvectors,
            eigenvalues = eigenvalues,
            proportionVariance = proportionVariance[c(1, 1, 2)]
        )),
        "The proportions of variance explained do not match the eigenvalues."
    )
})

test_that(".runAPCA", {
    set.seed(3)
    ## with svd
    m1 <- matrix(runif(100), 10, 10)
    m2 <- matrix(runif(100), 10, 10)
    expect_identical(
        .runAPCA(m1, m2, .svdWrapper),
        .svdWrapper(m1 + m2)
    )
    ## with nipals
    i <- sample(1:length(m1), length(m1) / 2)
    m1[i] <- m2[i] <- NA
    exp <- expect_warning(
        .nipalsWrapper(m1 + m2),
        regexp = "Stopping after 500 iterations"
    )
    test <- expect_warning(
        .runAPCA(m1, m2, .nipalsWrapper),
        regexp = "Stopping after 500 iterations"
    )
    expect_identical(test, exp)
    ## Test ... arguments
    exp <- expect_warning(
        .nipalsWrapper(m1 + m2, maxiter = 3),
        regexp = "Stopping after 3 iterations"
    )
    test <- expect_warning(
        .runAPCA(m1, m2, .nipalsWrapper, maxiter = 3),
        regexp = "Stopping after 3 iterations"
    )
    expect_identical(test, exp)
})

test_that(".runASCA", {
    set.seed(1)
    ## with svd
    m1 <- matrix(runif(100) * c(1,10), 10, 10) ## add some signal to avoid warnings that are difficult to reproduce
    m2 <- matrix(runif(100), 10, 10)
    expect_identical(
        .runASCA(m1, m2, .svdWrapper),
        .svdWrapper(m1)
    )
    ## with nipals
    i <- sample(1:length(m1), length(m1) / 2)
    m1[i] <- m2[i] <- NA
    expect_warning(
        expect_identical(
            .runASCA(m1, m2, .nipalsWrapper),
            .nipalsWrapper(m1)
        ),
        regexp = "Stopping after 500 iterations"
    )
    ## Test ... arguments
    expect_warning(
        expect_identical(
            .runASCA(m1, m2, .nipalsWrapper, maxiter = 3),
            .nipalsWrapper(m1, maxiter = 3)
        ),
        regexp = "Stopping after 3 iterations"
    )
})

test_that(".runASCA.E", {
    set.seed(1234)
    ## with svd
    m1 <- matrix(runif(100), 10, 10)
    m2 <- matrix(runif(100), 10, 10)
    exp <- .svdWrapper(m1)
    m1res <- m1 + m2
    m1res[is.na(m1res)] <- 0
    exp$scores <- crossprod(m1res, exp$eigenvectors)
    exp$eigenvalues <- structure(rep(NA, length(exp$eigenvalues)), .Names = c("PC1", "PC2"))
    exp$proportionVariance <- structure(rep(NA, length(exp$proportionVariance)), .Names = c("PC1", "PC2"))
    expect_identical(
        .runASCA.E(m1, m2, .svdWrapper),
        exp
    )
    ## with nipals
    i <- sample(1:length(m1), length(m1) / 2)
    m1[i] <- m2[i] <- NA
    exp <- .nipalsWrapper(m1)
    m1res <- m1 + m2
    m1res[is.na(m1res)] <- 0
    exp$scores <- crossprod(m1res, exp$eigenvectors)
    exp$eigenvalues <- structure(rep(NA, length(exp$eigenvalues)), .Names = c("PC1", "PC2"))
    exp$proportionVariance <- structure(rep(NA, length(exp$proportionVariance)), .Names = c("PC1", "PC2"))
    expect_identical(
        .runASCA.E(m1, m2, .nipalsWrapper),
        exp
    )
    ## Test ... arguments
    expect_warning(
        exp <- .nipalsWrapper(m1, ncomp = 10, maxiter = 30),
        regexp = "Stopping after 30 iterations"
    )
    m1res <- m1 + m2
    m1res[is.na(m1res)] <- 0
    exp$scores <- crossprod(m1res, exp$eigenvectors)
    exp$eigenvalues <- structure(rep(NA, length(exp$eigenvalues)), .Names = paste0("PC", 1:10))
    exp$proportionVariance <- structure(rep(NA, length(exp$proportionVariance)), .Names = paste0("PC", 1:10))
    expect_warning(
        expect_identical(
            .runASCA.E(m1, m2, .nipalsWrapper, ncomp = 10, maxiter = 30),
            exp
        ),
        regexp = "Stopping after 30 iterations"
    )
})

test_that(".nipalsWrapper", {
    set.seed(1234)
    m <- matrix(runif(120), 12, 10,
                dimnames = list(paste0("row", 1:12), paste0("col", 1:10)))
    ## ncomp larger than max dimension = error
    expect_error(
        .nipalsWrapper(m, ncomp = 11),
        "'ncomp' cannot exceeded number of features or samples, whichever is smallest."
    )
    ## NIPALS without missing values, default parameters
    require(nipals)
    exp <- nipals(t(m), startcol = function(x) sum(!is.na(x)),
                  scale = FALSE, ncomp = 2)
    expect_identical(
        .nipalsWrapper(m),
        list(
            scores = matrix(exp$scores %*% diag(exp$eig), ncol = 2,
                            dimnames = list(colnames(m), c("PC1", "PC2"))),
            eigenvectors = exp$loadings,
            eigenvalues = structure(exp$eig^2 / (nrow(m) - 1),
                                    .Names = c("PC1", "PC2")),
            proportionVariance = structure(exp$R2, .Names = c("PC1", "PC2"))
        )
    )
    ## NIPALS with missing values, default parameters
    i <- sample(1:length(m), length(m) / 2)
    m[i] <- NA
    expect_warning(
        exp <- nipals(t(m), startcol = function(x) sum(!is.na(x)),
                      scale = FALSE, ncomp = 2),
        "Stopping after 500 iterations for PC 1."
    )
    expect_warning(expect_identical(
        .nipalsWrapper(m),
        list(
            scores = matrix(exp$scores %*% diag(exp$eig), ncol = 2,
                            dimnames = list(colnames(m), c("PC1", "PC2"))),
            eigenvectors = exp$loadings,
            eigenvalues = structure(exp$eig ^ 2 / (nrow(m) - 1),
                                    .Names = c("PC1", "PC2")),
            proportionVariance = structure(exp$R2, .Names = c("PC1", "PC2"))
        )
    ),"Stopping after 500 iterations for PC 1.")
    ## NIPALS with missing values, one row is all NA
    m[1, ] <- NA
    exp <- nipals(t(m[-1, ]), startcol = function(x) sum(!is.na(x)),
                  scale = FALSE, ncomp = 2)
    expect_identical(
        .nipalsWrapper(m),
        list(
            scores = matrix(exp$scores %*% diag(exp$eig), ncol = 2,
                            dimnames = list(colnames(m), c("PC1", "PC2"))),
            eigenvectors = exp$loadings,
            eigenvalues = structure(exp$eig ^ 2 / (nrow(m[-1, ]) - 1),
                                    .Names = c("PC1", "PC2")),
            proportionVariance = structure(exp$R2, .Names = c("PC1", "PC2"))
        )
    )
    ## NIPALS with missing values, change default parameters
    expect_warning(
        exp <- nipals(t(m[-1, ]), startcol = function(x) sum(!is.na(x)),
                      scale = TRUE, center = FALSE, maxiter = 3,
                      ncomp = 2, tol = 1E-2),
        "Stopping after 3 iterations for PC 1."
    )
    expect_warning(expect_identical(
        .nipalsWrapper(m, scale = TRUE, center = FALSE, maxiter = 3,
                       ncomp = 2, tol = 1E-2),
        list(
            scores = matrix(exp$scores %*% diag(exp$eig), ncol = 2,
                            dimnames = list(colnames(m), c("PC1", "PC2"))),
            eigenvectors = exp$loadings,
            eigenvalues = structure(exp$eig ^ 2 / (nrow(m[-1, ]) - 1),
                                    .Names = c("PC1", "PC2")),
            proportionVariance = structure(exp$R2, .Names = c("PC1", "PC2"))
        )
    ), "Stopping after 3 iterations for PC 1.")
})

test_that(".svdWrapper", {
    set.seed(1234)
    m <- matrix(runif(120), 12, 10,
                dimnames = list(paste0("row", 1:12), paste0("col", 1:10)))
    ## missing values = error
    mna <- m
    mna[1] <- NA
    expect_error(
        .svdWrapper(mna),
        "svd cannot deal with missing values. Use 'algorithm = \"nipals\"' instead."
    )
    ## ncomp larger than max dimension = error
    expect_error(
        .svdWrapper(m, ncomp = 11),
        "'ncomp' cannot exceeded number of features or samples, whichever is smallest."
    )
    ## SVD, default parameters
    exp <- svd(scale(t(m), center = TRUE, scale = FALSE))
    eig <- exp$d^2 / (nrow(m) - 1)
    expect_identical(
        .svdWrapper(m),
        list(
            scores = matrix((exp$u %*% diag(exp$d))[, 1:2], ncol = 2,
                            dimnames = list(colnames(m), c("PC1", "PC2"))),
            eigenvectors = matrix((exp$v)[, 1:2], ncol = 2,
                                  dimnames = list(rownames(m), c("PC1", "PC2"))),
            eigenvalues = structure(eig[1:2],
                                    .Names = c("PC1", "PC2")),
            proportionVariance = structure(eig[1:2] / sum(eig),
                                           .Names = c("PC1", "PC2"))
        )
    )
    ## SVD, change parameters
    exp <- svd(scale(t(m), center = FALSE, scale = TRUE))
    eig <- exp$d^2 / (nrow(m) - 1)
    expect_identical(
        .svdWrapper(m, ncomp = 3, center = FALSE, scale = TRUE),
        list(
            scores = matrix((exp$u %*% diag(exp$d))[, 1:3], ncol = 3,
                            dimnames = list(colnames(m), c("PC1", "PC2", "PC3"))),
            eigenvectors = matrix((exp$v)[, 1:3], ncol = 3,
                                  dimnames = list(rownames(m), c("PC1", "PC2", "PC3"))),
            eigenvalues = structure(eig[1:3],
                                    .Names = c("PC1", "PC2", "PC3")),
            proportionVariance = structure(eig[1:3] / sum(eig),
                                           .Names = c("PC1", "PC2", "PC3"))
        )
    )
})

test_that(".formatPcaList", {
    ## Empty list return an empty list
    expect_identical(
        .formatPcaList(list()),
        List()
    )
    ## List with wrong names (do not end with bySample or byFeature) = error
    expect_error(
        .formatPcaList(List(foo = 1)),
        "Unexpected names in PCA result list."
    )
    expect_error(
        .formatPcaList(List(feature = 1)),
        "Unexpected names in PCA result list."
    )
    expect_error(
        .formatPcaList(List(byFeature = 1, foo = 1)),
        "Unexpected names in PCA result list."
    )
    expect_error(
        .formatPcaList(List(byFeature_var1 = 1, bySample_var1 = 1)),
        "Unexpected names in PCA result list."
    )
    ## Reformat list with only 1 PCA result, ie 1 bySample and 1 byFeature
    expect_identical(
        .formatPcaList(List(Var1_byFeature = 1, Var1_bySample = 2)),
        List(bySample = List(Var1 = 2),
             byFeature = List(Var1 = 1))
    )
    ## Reformat list with 2 PCA results
    expect_identical(
        .formatPcaList(List(Var1_byFeature = 1, Var1_bySample = 2,
                            Var2_byFeature = 3, Var2_bySample = 4)),
        List(bySample = List(Var1 = 2, Var2 = 4),
             byFeature = List(Var1 = 1, Var2 = 3))
    )
})

## ---- Aggregation functions ----


test_that("scpComponentAggregate", {
    table1 = DataFrame(PC1 = 1:10, PC2 = 1:10, prot = rep(1:2, each = 5),
                       annotationDrop = paste("foo", 1:10),
                       annotationKeep = paste("foo", rep(1:2, each = 5)))
    metadata(table1)$foo <- "bar"
    table2 = DataFrame(PC1 = 1:10, PC2 = 1:10, prot = rep(1:2, each = 5),
                       annotationDrop = paste("foo", 1:10),
                       annotationKeep = paste("foo", rep(1:2, each = 5)))
    ex <- List(
        table1 = table1,
        table2 = table2
    )
    ## fcol not found = error
    expect_error(
        scpComponentAggregate(componentList = ex, fcol = "foo"),
        "'foo' not found in list element.s.: table1, table2."
    )
    ## Default usage
    exp1 <- exp2 <- DataFrame(
        PC1 = c(3, 8), PC2 = c(3, 8), prot = 1:2,
        annotationKeep = c("foo 1", "foo 2"), .n = c(5L, 5L),
        row.names = 1:2)
    metadata(exp1)$foo <- "bar"
    expect_message(expect_identical(
        scpComponentAggregate(componentList = ex, fcol = "prot"),
        List(
            table1 = exp1,
            table2 = exp2
        )
    ), "Components may no longer be orthogonal after aggregation.")
    ## Change aggregation function
    setToOne <- function(x, warning = FALSE) {
        if (warning) warning("This is a warning")
        return(rep(1, ncol(x)))
    }
    exp1 <- exp2 <- DataFrame(
        PC1 = c(1, 1), PC2 = c(1, 1), prot = 1:2,
        annotationKeep = c("foo 1", "foo 2"), .n = c(5L, 5L),
        row.names = 1:2)
    metadata(exp1)$foo <- "bar"
    expect_message(expect_identical(
        scpComponentAggregate(componentList = ex, fcol = "prot", fun = setToOne),
        List(
            table1 = exp1,
            table2 = exp2
        )
    ), "Components may no longer be orthogonal after aggregation.")
    ## Test ... argument, eg throw a warning through setToOne
    expect_warning(
        expect_message(
            expect_identical(
                scpComponentAggregate(componentList = ex, fcol = "prot", fun = setToOne,
                                      warning = TRUE),
                scpComponentAggregate(componentList = ex, fcol = "prot", fun = setToOne)
            ), "Components may no longer be orthogonal after aggregation."
        ), "This is a warning"
    )
})

## ---- Plotting functions ----

test_that("scpComponentPlot", {
    m <- matrix(
        1, 10, 15, dimnames = list(paste0("row", 1:10), paste0("col", 1:15))
    )
    se <- SummarizedExperiment(assays = List(assay = m))
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$numeric <- as.vector(scale(1:ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    assay(se) <- sweep(assay(se), 2, se$numeric * 2, "+")
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric)
    caRes <- scpComponentAnalysis(se, method = "APCA", ncomp = 5 )
    ## Select components out of bound = error
    expect_error(
        scpComponentPlot(caRes$bySample, comp = 5:6),
        "'comp' is out of bounds."
    )
    ## Default plot
    vdiffr::expect_doppelganger(
        "scpComponentPlot bySample unmodelled",
        scpComponentPlot(caRes$bySample)[[1]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentPlot bySample residuals",
        scpComponentPlot(caRes$bySample)[[2]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentPlot bySample APCA_condition",
        scpComponentPlot(caRes$bySample)[[3]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentPlot byFeature unmodelled",
        scpComponentPlot(caRes$byFeature)[[1]]
    )
    ## Change components
    vdiffr::expect_doppelganger(
        "scpComponentPlot change comp",
        scpComponentPlot(caRes$bySample, comp = 4:5)[[3]]
    )
    ## Point params
    se$cell <- colnames(se)
    caRes$bySample <- scpAnnotateResults(caRes$bySample, colData(se), by = "cell")
    vdiffr::expect_doppelganger(
        "scpComponentPlot change pointParams",
        scpComponentPlot(caRes$bySample,
                         pointParams = list(aes(colour = condition,
                                                size = numeric))
        )[[3]]
    )
    ## maxLevel
    se$conditionMore <- factor(1:ncol(se))
    caRes$bySample <- scpAnnotateResults(caRes$bySample, colData(se), by = "cell")
    vdiffr::expect_doppelganger(
        "scpComponentPlot change maxLevel",
        scpComponentPlot(caRes$bySample,
                         pointParams = list(aes(colour = conditionMore,
                                                size = numeric)),
                         maxLevels = 3
        )[[3]]
    )
})

test_that(".plotPCA", {
    ## comp is not length 2 = error
    expect_error(
        .plotPCA(comp = 1:3),
        "length.comp. == 2 is not TRUE"
    )
    expect_error(
        .plotPCA(comp = 1),
        "length.comp. == 2 is not TRUE"
    )
    pcs <- DataFrame(PC1 = rep(1:5, 3),
                     PC2 = rep(5:1, 3) + rep((1:5)/10, each = 3),
                     PC3 = rep(1:3, 5),
                     colourVar = factor(rep(1:3, 5)),
                     sizeVar = 1:15)
    ## Initial test
    vdiffr::expect_doppelganger(
        ".plotPCA initial config",
        .plotPCA(pcs, comp = 1:2, proportionVariance = c(0.5, 0.25),
                 pointParams = list())
    )
    ## Other components
    vdiffr::expect_doppelganger(
        ".plotPCA change components",
        .plotPCA(pcs, c(1, 3), c(0.15, 0.35), list())
    )
    ## Proportio variance is NA = no percentage displayed on axis
    vdiffr::expect_doppelganger(
        ".plotPCA pVar is NA",
        .plotPCA(pcs, comp = 1:2, proportionVariance = c(NA, NA),
                 pointParams = list())
    )
    ## Change pointParams
    vdiffr::expect_doppelganger(
        ".plotPCA change pointParams",
        .plotPCA(pcs, comp = 1:2, proportionVariance = c(0.5, 0.25),
                 pointParams = list(aes(colour = colourVar, size = sizeVar),
                                    alpha = 0.5))
    )
})

test_that(".pcaAxisLabels", {
    ## standard use
    expect_identical(
        .pcaAxisLabels(c(0.5, 0.25), 1:2),
        list(x = "PC1 (50%)",
             y = "PC2 (25%)")
    )
    ## proportion is NA
    expect_identical(
        .pcaAxisLabels(c(NA, NA), 1:2),
        list(x = "PC1",
             y = "PC2")
    )
    ## change comp
    expect_identical(
        .pcaAxisLabels(c(0.7, 0.1, 0.04, 0.01), 4:3),
        list(x = "PC4 (1%)",
             y = "PC3 (4%)")
    )
    expect_identical(
        .pcaAxisLabels(rep(NA, 4), 4:3),
        list(x = "PC4",
             y = "PC3")
    )
    ## small fraction
    expect_identical(
        .pcaAxisLabels(c(0.01, 0.001), 1:2),
        list(x = "PC1 (1%)",
             y = "PC2 (0.1%)")
    )
})

test_that(".trimPlotLevels", {
    df <- data.frame(
        x = 1:10,
        y = 1:10,
        numeric = 1:10,
        categorical = as.factor(1:10)
    )
    ## No colour = no effect
    pl <- ggplot(df) + aes(x = x, y = y) + geom_point()
    expect_identical(
        .trimPlotLevels(pl = pl, maxLevels = 1),
        pl
    )
    ## No categorcial variable = no effect
    pl <- ggplot(df) + aes(x = x, y = y, colour = numeric) + geom_point()
    expect_identical(
        .trimPlotLevels(pl, 1),
        pl
    )
    ## maxLevel smaller than number of levels = no effect
    pl <- ggplot(df) + aes(x = x, y = y, colour = categorical) + geom_point()
    expect_identical(
        .trimPlotLevels(pl, 100),
        pl
    )
    ## trim levels
    vdiffr::expect_doppelganger(
        ".trimPlotLevels",
        .trimPlotLevels(pl, 3)
    )
})

test_that("scpComponentBiplot", {
    m <- matrix(
        1, 10, 15, dimnames = list(paste0("row", 1:10), paste0("col", 1:15))
    )
    se <- SummarizedExperiment(assays = List(assay = m))
    se$condition <- as.factor(rep(1:2, length.out = ncol(se)))
    se$numeric <- as.vector(scale(1:ncol(se)))
    assay(se)[, se$condition == 2] <- 1:10 * assay(se)[, se$condition == 2]
    assay(se) <- sweep(assay(se), 2, se$numeric * 2, "+")
    set.seed(124)
    assay(se)[sample(1:length(assay(se)), length(assay(se))/2)] <- NA
    se <- scpModelWorkflow(se, formula = ~ 1 + condition + numeric)
    caRes <- scpComponentAnalysis(se, method = "APCA", ncomp = 5 )
    caRes$bySample <- scpAnnotateResults(
        caRes$bySample, data.frame(cell = colnames(se),
                                   colData(se)),
        by = "cell")
    caRes$byFeature <- scpAnnotateResults(
        caRes$byFeature, data.frame(feature = paste0("row", 1:10),
                                    annot = 1:10),
        by = "feature")
    ## Test default params
    vdiffr::expect_doppelganger(
        "scpComponentBiplot unmodelled default",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature
        )[["unmodelled"]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentBiplot residuals default",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature
        )[["residuals"]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentBiplot APCA_condition default",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature
        )[["APCA_condition"]]
    )
    vdiffr::expect_doppelganger(
        "scpComponentBiplot APCA_numeric default",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature
        )[["APCA_numeric"]]
    )
    ## Change comp
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change comp",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            comp = c(1, 3)
        )[["unmodelled"]]
    )
    ## Change pointParams
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change pointParams",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            pointParams = list(aes(colour = condition, size = numeric))
        )[["unmodelled"]]
    )
    ## Change arrowParams
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change arrowParams",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            arrowParams = list(mapping = aes(colour = annot, linewidth = annot))
        )[["unmodelled"]]
    )
    ## Change labelParams
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change labelParams",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            labelParams = list(mapping = aes(colour = annot, size = annot))
        )[["unmodelled"]]
    )
    ## Change textBy
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change textBy",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            textBy = "annot"
        )[["unmodelled"]]
    )
    ## Change top
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change top",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            top = 5
        )[["unmodelled"]]
    )
    ## Change maxLevels
    vdiffr::expect_doppelganger(
        "scpComponentBiplot change maxLevels",
        scpComponentBiplot(
            scoreList = caRes$bySample, eigenvectorList = caRes$byFeature,
            pointParams = list(aes(colour = condition)), maxLevels = 1
        )[["unmodelled"]]
    )
})

test_that(".scaleComponentsToUnity", {
    table <- List(
        table1 = DataFrame(PC1 = 1:10, PC2 = 1:10,
                           annotation = "foo"),
        table2 = DataFrame(PC1 = seq(0, 1, 0.1), PC2 = seq(0, 1, 0.1),
                           annotation = "foo")
    )
    expect_equal(
        .scaleComponentsToUnity(table, compNames = c("PC1", "PC2")),
        List(
            table1 = DataFrame(PC1 = (1:10)/sqrt(200), PC2 = (1:10)/sqrt(200),
                               annotation = "foo"),
            table2 = DataFrame(PC1 = seq(0, 1, 0.1) / sqrt(2), PC2 = seq(0, 1, 0.1) / sqrt(2),
                               annotation = "foo")
        ),
        tolerance = 1E-15
    )
})

test_that(".addEigenArrows", {
    df <- data.frame(
        x = seq(-1, 1, length.out = 10),
        y = seq(1, -1, length.out = 10),
        numeric = 1:10,
        categorical = as.factor(1:10)
    )
    eigenvectors <- DataFrame(
        PC1 = sin(seq(-1, 1, length.out = 15)),
        PC2 = cos(seq(-1, 1, length.out = 15)),
        PC3 = seq(-1, 1, length.out = 15)^2,
        name = paste0("foo", 1:15),
        name2 = paste0("bar", 1:15)
    )
    pl <- ggplot(df) + aes(x = x, y = y) + geom_point()
    ## Standard usage
    vdiffr::expect_doppelganger(
        ".addEigenArrows standard case",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = 1:2, textBy = "name", top = 10,
                        arrowParams = list(), labelParams = list())
    )
    ## Change comp
    vdiffr::expect_doppelganger(
        ".addEigenArrows change comp",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name", top = 10,
                        arrowParams = list(), labelParams = list())
    )
    ## Change textBy
    vdiffr::expect_doppelganger(
        ".addEigenArrows change textBy",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name2", top = 10,
                        arrowParams = list(), labelParams = list())
    )
    ## Change top
    vdiffr::expect_doppelganger(
        ".addEigenArrows change top",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name", top = 5,
                        arrowParams = list(), labelParams = list())
    )
    vdiffr::expect_doppelganger(
        ".addEigenArrows top is zero",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name", top = 0,
                        arrowParams = list(), labelParams = list())
    )
    vdiffr::expect_doppelganger(
        ".addEigenArrows top is max",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name", top = 100,
                        arrowParams = list(), labelParams = list())
    )
    ## Change arrowParams
    vdiffr::expect_doppelganger(
        ".addEigenArrows change arrowParams",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name2", top = 10,
                        arrowParams = list(aes(colour = name), linewidth = 2),
                        labelParams = list())
    )
    ## Change labelParams
    vdiffr::expect_doppelganger(
        ".addEigenArrows change labelParams",
        .addEigenArrows(pl = pl, eigenvectors = eigenvectors,
                        comp = c(1, 3), textBy = "name", top = 10,
                        arrowParams = list(),
                        labelParams = list(aes(colour = name), size = 2))
    )
})

test_that(".formatEigenVectors", {
    eigenvectors <- DataFrame(
        PC1 = sin(seq(-1, 1, length.out = 15)) / (1:15),
        PC2 = cos(seq(-1, 1, length.out = 15)) / (1:15),
        PC3 = (1:15)^2,
        name = paste0("foo", 1:15),
        name2 = paste0("bar", 1:15)
    )
    ## Text by not in table = error
    expect_error(
        .formatEigenVectors(
            eigenvectors = eigenvectors, textBy = "foo",
        ),
        "'foo' not found in component tables. Use 'scpAnnotateResults..' to add custom annotations."
    )
    ## Test an initial config
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = 1:2, textBy = "name",
            top = 10
        ),
        data.frame(eigenvectors[1:10,], label = eigenvectors$name[1:10])
    )
    ## change comp
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = c(1, 3), textBy = "name",
            top = 10
        ),
        data.frame(eigenvectors[15:6,], label = eigenvectors$name[15:6])
    )
    ## change textby
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = 1:2, textBy = "name2",
            top = 10
        ),
        data.frame(eigenvectors[1:10,], label = eigenvectors$name2[1:10])
    )
    ## change top
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = 1:2, textBy = "name",
            top = 5
        ),
        data.frame(eigenvectors[1:5,], label = eigenvectors$name[1:5])
    )
    ## top is zero = return empty table
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = 1:2, textBy = "name",
            top = 0
        ),
        data.frame(eigenvectors[0, , drop = FALSE])
    )
    ## top is higher or equal than nrow input table = no effect
    expect_identical(
        .formatEigenVectors(
            eigenvectors = eigenvectors, comp = 1:2, textBy = "name",
            top = 100
        ),
        data.frame(eigenvectors, label = eigenvectors$name)
    )
})

test_that(".formatParamsMapping", {
    ## aesParams is not a list = error
    expect_error(
        .formatParamsMapping(aesParams = aes(x = 1)),
        "'labelParams' and 'arrowParams' must be a list."
    )
    ## x xend y yend are user provided = error
    ## outside aes()
    expect_error(
        .formatParamsMapping(aesParams = list(x = 1)),
        "'x' cannot be user-provided."
    )
    expect_error(
        .formatParamsMapping(aesParams = list(x = 1, y = 1, xend = 1, yend = 1)),
        "'x', 'xend', 'y', 'yend' cannot be user-provided."
    )
    ## inside aes()
    expect_error(
        .formatParamsMapping(aesParams = list(aes(x = 1))),
        "'x' cannot be user-provided."
    )
    expect_error(
        .formatParamsMapping(aesParams = list(aes(x = 1, y = 1, xend = 1, yend = 1))),
        "'x', 'xend', 'y', 'yend' cannot be user-provided."
    )
    ## empty list = return list with imposed aes
    ## for arrowParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(), comp = 1:2, isArrow = TRUE),
        list(mapping = aes(x = 0, y = 0, xend = .data[["PC1"]],
                           yend = .data[["PC2"]]))
    )
    ## for labelParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(), 1:2, FALSE),
        list(mapping = aes(x = .data[["PC1"]], y = .data[["PC2"]],
                           label = .data$label))
    )
    ## non empty list only adapt the mapping
    ## for arrowParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(colour = "red"), comp = 1:2, isArrow = TRUE),
        list(colour = "red",
             mapping = aes(x = 0, y = 0, xend = .data[["PC1"]],
                           yend = .data[["PC2"]])
        )
    )
    ## for labelParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(colour = "red"), 1:2, FALSE),
        list(colour = "red",
             mapping = aes(x = .data[["PC1"]], y = .data[["PC2"]],
                           label = .data$label))
    )
    ## non empty list and other mappings are present
    ## for arrowParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(
            colour = "red",
            mapping = aes(size = annot)),
            comp = 1:2, isArrow = TRUE),
        list(colour = "red",
             mapping = aes(x = 0, y = 0, xend = .data[["PC1"]],
                           yend = .data[["PC2"]], size = annot)
        )
    )
    ## for labelParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(aesParams = list(
            colour = "red",
            mapping = aes(size = annot)),
            1:2, FALSE),
        list(colour = "red",
             mapping = aes(x = .data[["PC1"]], y = .data[["PC2"]],
                           label = .data$label, size = annot))
    )
    ## if no names, first is "mapping" and other are empty names
    ## for arrowParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(list(aes(colour = annot), "red"), 1:2, TRUE),
        list(mapping = aes(x = 0, y = 0, xend = .data[["PC1"]],
                           yend = .data[["PC2"]], colour = annot),
             "red"
        )
    )
    ## for labelParams
    expect_equal( ## equal because quosures have different environment, as expected
        .formatParamsMapping(list(aes(colour = annot), "red"), 1:2, FALSE),
        list(mapping = aes(x = .data[["PC1"]], y = .data[["PC2"]],
                           label = .data$label, colour = annot),
             "red"
        )
    )
})
