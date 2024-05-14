
## ---- ScpModel workflow ----


##' @name ScpModel-Workflow
##'
##' @title Modelling single-cell proteomics data
##'
##' @description
##'
##' Function to estimate a linear model for each feature (peptide or
##' protein) of a single-cell proteomics data set.
##'
##' @section Input data:
##'
##' The main input is `object` that inherits from the
##' `SingleCellExperiment` class. The quantitative data will be
##' retrieve using `assay(object)`. If `object` contains multiple
##' assays, you can specify which assay to take as input thanks to the
##' argument `i`, the function will then assume `assay(object, i)` as
##' quantification input .
##'
##' The objective of modelling single-cell proteomics data is to
##' estimate, for each feature (peptide or protein), the effect of
##' known cell annotations on the measured intensities. These annotations
##' may contain biological information such as the cell line,
##' FACS-derived cell type, treatment, etc. We also highly recommend
##' including technical information, such as the MS acquisition run
##' information or the chemical label (in case of multiplexed
##' experiments). These annotation must be available from
##' `colData(object)`. `formula` specifies which annotations to use
##' during modelling.
##'
##' @section Data modelling workflow:
##'
##' The modelling worflow starts with generating a model matrix for
##' each feature given the `colData(object)` and `formula`. The model
##' matrix for peptide \eqn{i}, denoted \eqn{X_i}, is adapted to the
##' pattern of missing values (see section below). Then, the functions
##' fits the model matrix against the quantitative data. In other
##' words, the function determines for each feature \eqn{i} (row in
##' the input data) the contribution of each variable in the model.
##' More formally, the general model definition is:
##'
##' \deqn{Y_i = \beta_i X^T_{(i)} + \epsilon_i}
##'
##' where \eqn{Y} is the feature by cell quantification matrix,
##' \eqn{\beta_i} contains the estimated coefficients for feature
##' \eqn{i} with as many coefficients as variables to estimate,
##' \eqn{X^T_{(i)}} is the model matrix generated for feature \eqn{i},
##' and \eqn{\epsilon} is the feature by cell matrix with
##' residuals.
##'
##' The coefficients are estimated using penalized least squares
##' regression. Next, the function computes the residual matrix and
##' the effect matrices. An effect matrix contains the data that is
##' captured by a given cell annotation. Formally, for each feature
##' \eqn{i}:
##'
##' \deqn{\hat{M^f_i} = \hat{\beta^f_i} X^{fT}_{(i)} }
##'
##' where \eqn{\hat{M^f}} is a cell by feature matrix containing the
##' variables associated to annotation \eqn{f}, \eqn{\hat{\beta^f_i}}
##' are the estimated coefficients associated to annotation \eqn{f}
##' and estimated for feature \eqn{i}, and \eqn{X^{fT}_{(i)}} is the
##' model matrix for peptide \eqn{i} containing only the variables to
##' annotation \eqn{f}.
##'
##' All the results are stored in an [ScpModel] object which is stored
##' in the `object`'s metadata. Note that multiple models can be
##' estimated for the same `object`. In that case, provide the `name`
##' argument to store the results in a separate `ScpModel`.
##'
##' @section Feature filtering:
##'
##' The proportion of missing values for each features is high in
##' single-cell proteomics data. Many features can typically contain
##' more coefficients to estimate than observed values. These features
##' cannot be estimated and will be ignored during further steps.
##' These features are identified by computing the ratio between the
##' number of observed values and the number of coefficients to
##' estimate. We call it the **n/p ratio**. Once the model is
##' estimated, use `scpModelFilterPlot(object)` to explore the
##' distribution of n/p ratios across the features. You can also
##' extract the n/p ratio for each feature using
##' `scpModelFilterNPRatio(object)`. By default, any feature that has
##' an n/p ratio lower than 1 is ignored. However, feature with an
##' n/p ratio close to 1 may lead to unreliable outcome because there
##' are not enough observed data. You could consider the n/p ratio as
##' the average number of replicate per coefficient to estimate.
##' Therefore, you may want to increase the n/p threshold. You can do
##' so using `scpModelFilter(object) <- npThreshold`.
##'
##' @section About missing values:
##'
##' The data modelling workflow is designed to take the presence of
##' missing values into account. We highly recommend to **not impute**
##' the data before modelling. Instead, the modelling approach will
##' ignore missing values and will generate a model matrix using only
##' the observed values for each feature. However, the model matrices
##' for some features may contain highly correlated variables, leading
##' to near singular designs. We include a small ridge penalty to
##' reduce numerical instability associated to correlated variables.
##'
##' @seealso
##'
##' - [ScpModel-class] for functions to extract information from the
##'   `ScpModel` object
##' - [ScpModel-VarianceAnalysis], [ScpModel-DifferentialAnalysis],
##'   [ScpModel-ComponentAnalysis] to explore the model results
##' - [scpKeepEffect] and [scpRemoveBatchEffect] to perform batch
##'   correction for downstream analyses.
##'
##' @author Christophe Vanderaa, Laurent Gatto
##'
##' @example inst/examples/examples_ScpModel-Workflow.R
##'
NULL

##' @name ScpModel-Workflow
##'
##' @param object An object that inherits from the
##'     `SingleCellExperiment` class.
##'
##' @param formula A `formula` object controlling which variables are
##'     to be modelled.
##'
##' @param i A `logical`, `numeric` or `character` indicating which
##'     assay of `object` to use as input for modelling. Only a single
##'     assay can be provided. Defaults to the first assays.
##'
##' @param name A `character(1)` providing the name to use to store or
##'     retrieve the modelling results. When retrieving a model and
##'     `name` is missing, the name of the first model found in
##'     `object` is used.
##'
##' @param verbose A `logical(1)` indicating whether to print progress
##'     to the console.
##'
##' @export
scpModelWorkflow <- function(object, formula,
                             i = 1,
                             name = "model",
                             verbose = TRUE) {
    if (name %in% names(metadata(object))) {
        warning("An element called '", name, "' is already present ",
                "in the metadata. The associated content will be ",
                "overwritten.")
    }
    metadata(object)[[name]] <- ScpModel()
    scpModelFormula(object, name) <- formula
    scpModelInputIndex(object, name) <- i
    scpModelFilterThreshold(object, name) <- 1
    scpModelFitList(object, name) <- .runScpModel(object, name, verbose)
    object
}

## Internal functions that runs an ScpModel contained in an SE
## object. Formula and input data are determined based on a
## pre-initialised ScpModel object contained in the metadata of the
## SE object. For each row of the input data, the function will
## generate an scpModelFit object containing the fitting results. The
## run process consist of 2 steps:
## 1. adapt the model design (= generate the design matrix) based on
## the pattern of missing values
## 2. fit the model
## The functions returns a list of scpModelFit objects
## @param object A SummarizedExperiment object containing a
##     pre-initialised ScpModel object in the metadata.
## @param name A `character(1)` providing the name to use to retrieve
##     the pre-initialised ScpModel object in the metadata. When
##     missing, the name of the first model found in `object` is used.
## @param verbose A logical(1) indicating whether to show a progress
##     bar on the console.
##' @importFrom utils txtProgressBar setTxtProgressBar
##' @importFrom S4Vectors List
.runScpModel <- function(object, name, verbose) {
    coldata <- .checkAnnotations(object, name)
    Y <- scpModelInput(object, name, filtered = FALSE)
    if (verbose) pb <- txtProgressBar(max = nrow(Y), style = 3)
    scpModelFits <- lapply(seq_len(nrow(Y)), function(i) {
        if (verbose) setTxtProgressBar(pb, i)
        design <- .adaptModel(
            Y[i, ], coldata, scpModelFormula(object, name)
        )
        out <- .fitModel(
            Y[i, ], design, scpModelEffectNames(object, name)
        )
        out
    })
    if (verbose) close(pb)
    names(scpModelFits) <- rownames(Y)
    as(scpModelFits, "List")
}

## Internal functions that checks whether the ScpModel formula and the
##
##' @importFrom SummarizedExperiment colData
.checkAnnotations <- function(object, name) {
    formula <- scpModelFormula(object, name)
    coldata <- colData(object)[, all.vars(formula), drop = FALSE]
    if (is.null(rownames(coldata)))
        stop("'colData(object)' must have row names.")
    if (any(sapply(coldata, function(x) any(is.na(x)))))
        stop("Sample annotations (colData) cannot contain missing values.")
    .checkExperimentalDesignRank(formula, coldata)
    .formatCategoricalVariablesAsFactor(colData(object))
}

## Internal function that makes sure that the design matrix is not
## singular, either because there are more coefficients than
## observations or because of singular designs. If the design matrix
## does not pass the checks, an informative error is thrown. Otherwise,
## the function returns NULL. This functions checks the design matrix
## for the full data set, irrespective of missing values. It  makes
## sure there are no obvious experimental design mistakes.
## @param formula A formula object used to generate a design matrix.
## @param coldata A DataFrame object used to generate a design matrix.
##' @importFrom stats model.matrix
.checkExperimentalDesignRank <- function(formula, coldata) {
    X <- model.matrix(formula, data = droplevels(coldata))
    if (ncol(X) > nrow(X)) {
        stop("The design matrix is underdetermined, i.e. there are ",
             "more coefficients to fit than observations. Solve this ",
             "issue by simplifying the model formula.")
    }
    if (any(svd(X)$d < 1E-10)) {
        stop("The design matrix is (near) singular. This indicates ",
             "the experimental design contains colinear variables. ",
             "You should randomize your SCP experiment (cf Gatto ",
             "et al. 2023, Nat Comm). Solve this issue by ",
             "simplifying the model formula.")
    }
    return(NULL)
}

## ---- Adaptive modelling ----

## Internal function that generates the design matrix and adapts it to
## the pattern of missing values.
##
## Adaptations consists of:
## 1. Remove rows associated to missing observations.
## 2. Remove unused levels for variable encoded as factors.
## 3. Remove any variable that became constant upon removal of missing
##    values.
## 4. Scale any numeric variable.
##
## If after missing value removal, the data consists of 2 observations
## or less, an empty design matrix is returned. Otherwise, the
## function returns the adapted design matrix with factor levels
## stored in the attributes for later retrieval during inference.
##
## @param y A `numeric()` providing the response variable used for
##     assessing which observation is missing values.
## @param coldata A `DataFrame` object containing the sample
##     annotations used to generate the model matrix.
## @param formula A `formula` object used to generate the model matrix.
##
.adaptModel <- function(y, coldata, formula) {
    coldata <- coldata[!is.na(y), all.vars(formula), drop = FALSE]
    if (nrow(coldata) <= 2) {
        out <- matrix(nrow = nrow(coldata), ncol = 0)
        attr(out, "levels") <- List()
        return(out)
    }
    coldata <- droplevels(coldata)
    formula <- .dropConstantVariablesFromFormula(coldata, formula)
    coldata <- coldata[, all.vars(formula), drop = FALSE]
    coldata <- .scaleNumericalVariables(coldata)
    out <- model.matrix(
        formula, data = coldata, contrasts.arg = .modelContrasts(coldata)
    )
    attr(out, "levels") <- .modelLevels(coldata)
    out
}

## Internal function that formats any categorical variables (encoded
## as a factor, character or logical) as factors and drops unused
## levels. The functions returns the input table with format
## variables.
## @param x A DataFrame.
.formatCategoricalVariablesAsFactor <- function(x) {
    categoricalClass <- c("factor", "character", "logical")
    for (i in which(sapply(x, class) %in% categoricalClass)) {
        x[, i] <- factor(x[, i])
    }
    x
}

## Internal function that scales any numerical (numeric or integer)
## variable in a table. The functions returns a table with scaled
## numerical variables. When less than 2 observations, values cannot
## be scaled and the function return the input table.
## @param x A DataFrame.
.scaleNumericalVariables <- function(x) {
    if (nrow(x) < 2 || ncol(x) == 0) return(x)
    for (i in which(sapply(x, is.numeric))) {
        x[, i] <- scale(x[, i])
    }
    x
}

## Internal function that removes variable from a formula that are
## constant in a table. The function returns an updated formula.
##
## - A numerical variable is considered constant if its variance is 0.
## - A categorical variable is considered constant if it has less than
##   2 levels.
##
## @param x A DataFrame
## @param formula A formula
##
##' @importFrom stats reformulate terms var
.dropConstantVariablesFromFormula <- function(x, formula) {
    isConstant <- sapply(x, function(xx) {
        (is.numeric(xx) && var(xx) == 0) ||
        (is.factor(xx) && nlevels(xx) < 2)
    })
    dropped <- colnames(x)[isConstant]
    if (!length(dropped)) return(formula)
    termLabels <- labels(terms(formula))
    droppedPattern <- paste0(dropped, collapse = "|")
    sel <- !sapply(termLabels, grepl, pattern = droppedPattern)
    termLabels <- termLabels[sel]
    reformulate(
        c("1", termLabels), intercept = TRUE,
        env = attr(formula, ".Environment")
    )
}

## Internal function that returns a list of contrast specification
## for variables encoded as factor in a table, as requested by the
## `contrasts.arg` argument in `model.matrix()`.
## @param x A DataFrame table.
.modelContrasts <- function(x) {
    factors <- colnames(x)[sapply(x, class) == "factor"]
    if (length(factors) > 0) {
        contrasts <- list()
        contrasts[factors] <- "contr.sum"
    } else {
        contrasts <- NULL
    }
    contrasts
}

## Internal function that returns a List of factor levels for all
## factors in a table.
## @param x A DataFrame table.
.modelLevels <- function(x) {
    factors <- colnames(x)[sapply(x, class) == "factor"]
    if (!length(factors)) return(List())
    levs <- lapply(factors, function(f) levels(x[[f]]))
    names(levs) <- factors
    as(levs, "List")
}


## ---- Model estimation ----

## Internal function that fits a design matrix to a response variable.
## The function return the estimates through a ScpModelFit object.
##
## @param y A numeric() containing the response variable
## @param x A matrix() corresponding to the model matrix where each
##     column is a parameter to fit.
## @param effectNames A character() providing the names of the model
##     variables associated to the . This is only used to know which
##     effects can be computed
##
.fitModel <- function(y, x, effectNames) {
    model <- ScpModelFit(nrow(x), ncol(x))
    if (!ncol(x)) return(model)
    res <- .fitRidge(y[rownames(x)], x)
    scpModelFitCoefficients(model) <- res$coefficients
    scpModelFitResiduals(model) <- res$residuals
    scpModelFitEffects(model) <-
        .computeModelEffects(res$coefficients, x, effectNames)
    scpModelFitDf(model) <- res$df
    scpModelFitVar(model) <- res$var
    scpModelFitUvcov(model) <- res$uvcov
    scpModelFitLevels(model) <- attr(x, "levels")
    model
}

## Internal function to fit a vector y to the design matrix x given a
## ridge penalty given by lambda. The function returns the
## estimated coefficients, the residuals, the effective degrees of
## freedom, the residual variance and the unscaled variance-covariance
## matrix.
## @param y A numeric() containing the response variable
## @param x A matrix() corresponding to the model matrix where each
##     column is a parameter to fit.
## @param lalmbda A numeric(1) providing the regularisation penalty
##
## Note: when performing ridge regression, the dof are computed using
## the hat matrix:
##> hat <- x %*% tcrossprod(xtxInvLam, x)
##> df <- length(y) - sum(diag(2 * hat - tcrossprod(hat)))
## However, the hat matrix can become huge for large number of
## cells! Also computing the effective degrees of freedom takes most
## of the time... Since the lambda used is small, the effective df
## is close to n - p, so we don't compute the edf for now, but
## this can be changed setting effectiveDf = TRUE.
## See refs for background:
## https://doi.org/10.1186/1471-2105-12-372
## https://stats.stackexchange.com/a/44841
##
## Note: the degrees of freedom are often given by n - p - 1,
## where - 1 is to account for the degree of freedom taken by the
## estimation of the intercept. However, because intercept is in x, we
## compute df = n - p = lenght(y) - ncol(x)
##
## Note: I don't expose 'lambda' to the users because it
## is solely implemented to avoid numerical instability issues due to
## correlated variables. I expect tuning lambda for ridge
## regression to be beneficial, but we then need to implement more
## advanced functionality for fitting and selecting the best lambda.
## The difficulty here is that the best lambda will depend on the
## model and the input data and hence should be different for each
## feature.
.fitRidge <- function(y, x, lambda = 1E-3) {
    xtxInvLam <- solve(crossprod(x) + lambda * diag(ncol(x)))
    beta <- xtxInvLam %*% crossprod(x, y)
    df <- max(length(y) - ncol(x), 1)
    residuals <- y - x %*% beta
    list(coefficients = beta[, 1],
         residuals = residuals[, 1],
         df = df,
         var = sum(residuals^2) / df,
         uvcov = crossprod(x %*% xtxInvLam))
}

## Internal function that computes an effect vector a given
## variables. An effect vector is the fitted data by a given
## variable, i.e. $\hat{y} = X_{sel} \beta{sel}$, where $sel$ are the
## subset of model parameters related to a model variable. The
## functions returns the effects vector in a named numeric vector.
##
## @param beta A numeric() of estimated coefficient
## @param x The model matrix
## @param effect A character(1) providing the name of the variable for
##     which to compute effect matrices.
##
## Note: this should be adapted in case of interaction
.computeModelEffect <- function(beta, x, effect) {
    effect <- gsub("[(]", "\\\\(", effect)
    effect <- gsub("[)]", "\\\\)", effect)
    sel <- grep(paste0("^", effect, "\\d*"), names(beta))
    if (!length(sel)){
        return(NULL)
    }
    out <- tcrossprod(beta[sel], x[, sel, drop = FALSE])
    out[1, ]
}

## Internal function that computes the effect vectors for all
## requested variables (see `computeModelEffect()`). The
## functions returns the effects vectors in a named list.
##
## @param beta A numeric() of estimated coefficient
## @param x The model matrix
## @param effects A character() of variable names for which to compute
##     effect matrices
##
## Note: this should be adapted in case of interaction, namely the
## function should automatically detect variables containing an
## interaction operator ("*", ":") and define effects accordingly
## (include main effects or not?).
.computeModelEffects <- function(beta, x, effects) {
    out <- lapply(effects, function(e) {
        .computeModelEffect(beta, x, e)
    })
    names(out) <- effects
    out <- out[!sapply(out, is.null)]
    as(out, "List")
}

## ---- ScpModel filter plot ----

##' @name ScpModel-Workflow
##'
##' @import ggplot2
##'
##' @export
scpModelFilterPlot <- function(object, name) {
    message(
        "To change the threshold, use:\n",
        "scpModelFilterThreshold(object, name) <- threshold"
    )
    npRatio <- scpModelFilterNPRatio(object, name, FALSE)
    npThreshold <- scpModelFilterThreshold(object, name)
    .filterPlot(npRatio, npThreshold) +
        ggtitle(
            "Distribution of the n/p ratio",
            subtitle = .filterSubtitle(npRatio, npThreshold)
        )
}

## Internal function that generates a description string for the
## scpModelFilterPlot. This description contains the number of total
## features in the data and the number of estimable features. If
## threshold is greater than 1 (ie the model is not saturated), the
## description also includes the number of feature that pass N/P ratio
## filtering.
##
## @param npRatio A numeric() containing the N/P ratios for each feature
## @param threshold A non-negative numeric(1) that provides the N/P
##     ratio threshold used to filter features.
##
.filterSubtitle <- function(npRatio, threshold) {
    out <- paste(
        "Total features:", length(npRatio),
        "\nEstimated features (n/p >= 1):", sum(npRatio >= 1)
    )
    if (threshold > 1) {
        out <- paste0(
            out, "\nSelected features (n/p >= ", threshold, "): ",
            sum(npRatio >= threshold)
        )
    }
    out
}

.filterPlot <- function(npRatio, threshold) {
    breaks <- if (threshold > 1) c(1, threshold) else 1
    type <- cut(npRatio, breaks = c(-Inf, breaks, Inf))
    cols <- c(
        "inestimable" = "grey70", "estimated" = "grey40",
        "selected" = "dodgerblue"
    )
    levels(type) <- names(cols)
    df <- data.frame(npRatio = npRatio, type = type)
    ggplot(df) +
        aes(x = npRatio, fill = type) +
        geom_histogram(bins = 50) +
        scale_fill_manual(values = cols, name = "") +
        theme_minimal() +
        theme(legend.position = "top")
}
