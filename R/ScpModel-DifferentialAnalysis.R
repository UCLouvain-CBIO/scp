## ---- SCP Differential Analysis ----

##' @name ScpModel-DifferentialAnalysis
##'
##' @title Differential abundance analysis for single-cell proteomics
##'
##' @description
##'
##' Differential abundance analysis assess the statistical
##' significance of the differences observed between group of samples
##' of interest.
##'
##' @section Running the differential abundance analysis:
##'
##' `scpDifferentialAnalysis()` performs statistical inference by
##' means of a t-test on the estimatated parameters. There are 2 use
##' cases:
##'
##' 1. **Statistical inference for differences between 2 groups**
##'
##' You can **contrast** 2 groups of interest through the `contrasts`
##' argument. Multiple contrasts, that is multiple pairwise group
##' comparisons, can be performed. Therefore, `contrasts` must be
##' provided as a list where each element describes the comparison to
##' perform as a three-element character vector (see examples). The
##' first element is the name of the annotation variable that contains
##' the two groups to compare. This variable must be **categorical**.
##' The second element is the name of the reference group. The third
##' element is the name of the other group to compare against the
##' reference.
##'
##' 2. **Statistical inference for numerical variables**
##'
##' Numerical variables can be tested by providing the `coefficient`
##' argument, that is the name of the numerical annotation variable.
##'
##' The statistical tests in both use cases are conducted for each
##' feature independently. The p-values are adjusted using
##' [IHW::ihw()], where each test is weighted using the feature
##' intercept (that is the average feature intensity). The function
##' returns a list of `DataFrame`s with one table for each test
##' contrast and/or coefficient. It provides the adjusted p-values and
##' the estimates. For contrast, the estimates represent the estimated
##' log fold changes between the groups. For coefficients, the
##' estimates are the estimated slopes. Results are only provided for
##' features for which contrasts or coefficients are estimable, that
##' are features for which there is sufficient observations for
##' inference.
##'
##' @section Differential abundance at the protein level:
##'
##' `scpDifferentialAggregate()` combines the differential abundance
##' analysis results for groups of features. This is useful, for
##' example, to return protein-level results when data is modelled at
##' the peptide level. The function heavily relies on the approaches
##' implemented in [metapod::combineGroupedPValues()]. The p-values
##' are combined into a single value using one of the following
##' methods: Simes' method
##' (default), Fisher's method, Berger's method, Pearson's method,
##' minimum Holm's approach, Stouffer's Z-score method, and
##' Wilkinson's method. We refer to the `metapod` documentation for
##' more details on the assumptions underlying each approach. The
##' estimates are combined using the representative estimate, as
##' defined by `metapod`. Which estimate is representative depends on
##' the selected combination method. The function takes the list of
##' tables generated by `scpDifferentialAnalysis()` and returns a new
##' list of `DataFrame`s with aggregated results. Note that we cannot
##' meaningfully aggregate degrees of freedom. Those are hence removed
##' from the aggregated result tables.
##'
##' @section Volcano plots:
##'
##' [scpAnnotateResults()] adds annotations to the differential abundance
##' analysis results. The annotations are added to all elements of the
##' list returned by `()`. See the associated
##' man page for more information.
##'
##' `scpVolcanoPlot()` takes the list of tables generated by
##' `scpDifferentialAnalysis()` and returns a `ggplot2` scatter plot.
##' The plots show the adjusted p-values with respect to the estimate.
##' A horizontal bar also highlights the significance threshold
##' (defaults to 5%, `fdrLine`). The top (default 10) features with lowest
##' p-values are labeled on the plot. You can control which features
##' are labelled using the `top`, `by` and `decreasing` arguments.
##' Finally, you can change the point and label aesthetics thanks to
##' the `pointParams` and the `labelParams` arguments, respectively.
##'
##' @seealso
##'
##' - [ScpModel-Workflow] to run a model on SCP data upstream of
##'   differential abundance analysis.
##' - [scpAnnotateResults()] to annotate analysis of variance results.
##'
##' @author Christophe Vanderaa, Laurent Gatto
##'
##' @example inst/examples/examples_ScpModel-DifferentialAnalysis.R
##'
NULL

## ---- Analysis functions ----

##' @name ScpModel-DifferentialAnalysis
##'
##' @param object An object that inherits from the
##'     `SummarizedExperiment` class. It must contain an estimated
##'     `ScpModel` in its metadata.
##'
##' @param coefficients A `character()` vector with coefficient names
##'     to test. `coefficients` and `contrasts` cannot be both NULL.
##'
##' @param contrasts A `list()` where each element is a contrast to
##'     test. Each element must be a vector with 3 strings: 1. The
##'     name of a categorical variable to test; 2. The name of the
##'     reference group: 3. The name of the second group to contrast
##'     against the reference group. `coefficients` and `contrasts`
##'     cannot be both NULL.
##'
##' @param name A `character(1)` providing the name to use to retrieve
##'     the model results. When retrieving a model and `name` is
##'     missing, the name of the first model found in `object` is used.
##'
##' @export
scpDifferentialAnalysis <- function(object,
                                    coefficients = NULL,
                                    contrasts = NULL,
                                    name) {
    if (is.null(contrasts) & is.null(coefficients)) {
        stop("'contrasts' and 'coefficients' cannot be both NULL.")
    }
    if (missing(name)) name <- .defaultModelName(object)
    out <- List()
    if (!is.null(contrasts)) {
        out <- c(out, .scpDifferentialAnalysisOnContrast(
            object, contrasts, name
        ))
    }
    if (!is.null(coefficients)) {
        out <- c(out, .scpDifferentialAnalysisOnCoefficient(
            object, coefficients, name
        ))
    }
    as(out, "DataFrameList")
}

## Internal function that performs statistical inferences between 2
## groups defined by a modelled categorical variable, which we define
## as a contrast. The contrast are defined in a list where each
## element identifies the two groups to compare. The function returns
## a list of data tables as returned by .computeTTest().
##
## p-values are adjusted using the IHW approach given there is enough
## data. The weights required by IHW are assigned based on the model
## intercept, with features with heigher intercept are given more
## weight.
##
## @param object An object that inherits from the
##     `SummarizedExperiment` class. It must contain an estimated
##     `ScpModel` in its metadata.
## @param contrasts A `list()` where each element is a contrast to
##     test. Each element must be a vector with 3 strings: 1. The
##     name of a categorical variable to test; 2. The name of the
##     reference group: 3. The name of the second group to contrast
##     against the reference group. `coefficients` and `contrasts`
##     cannot be both NULL.
## @param name A `character(1)` providing the name to use to retrieve
##     the model results. When retrieving a model and `name` is
##     missing, the name of the first model found in `object` is used.
##
.scpDifferentialAnalysisOnContrast <- function(object, contrasts,
                                               name) {
    contrasts <- .checkContrasts(object, contrasts, name)
    df <- scpModelDf(object, name)
    weights <- scpModelIntercept(object, name)
    daTables <- lapply(contrasts, function(contrast) {
        est <- .contrastToEstimates(object, contrast, name)
        daTable <- .computeTTest(est$logFc, est$se, df, weights)
        metadata(daTable)$contrast <- contrast
        daTable[!is.na(daTable$Estimate), ]
    })
    names(daTables) <- names(contrasts)
    daTables
}

## Internal function that checks whether the provided contrasts
## unambiguously identify two groups to compare. The function will
## automatically assign names to each contrast and return a set of
## unique and valid contrast. Invalid contrast will throw an error.
## Validity conditions are:
## - contrasts must be a list
## - contrast can only contain elements of length 3
## - the first element of each contrast must be a character corresponding
##   to a modelled variable which must be categorical
## - The second and third element of each contrast must contain two
##   level names present in the modelled variables.
##
## Empty contrast are valid, the function then returns NULL.
##
## @param object An object that inherits from the
##     `SummarizedExperiment` class. It must contain an estimated
##     `ScpModel` in its metadata.
## @param contrasts A `list()` where each element is a contrast to
##     test. Each element must be a vector with 3 strings: 1. The
##     name of a categorical variable to test; 2. The name of the
##     reference group: 3. The name of the second group to contrast
##     against the reference group. `coefficients` and `contrasts`
##     cannot be both NULL.
## @param name A `character(1)` providing the name to use to retrieve
##     the model results. When retrieving a model and `name` is
##     missing, the name of the first model found in `object` is used.
##
.checkContrasts <- function(object, contrasts, name) {
    if (is.null(contrasts)) return(NULL)
    if (!is.list(contrasts) & !inherits(contrasts, "List")) {
        stop("'contrasts' must be a list.")
    }
    if (any(sapply(contrasts, length) != 3)) {
        stop("'contrasts' must be a list with elements of length 3.")
    }
    effects <- unique(sapply(contrasts, "[[", 1))
    if (any(mis <- !effects %in% scpModelEffectNames(object, name))) {
        stop(
            "Effect(s) '", paste(effects[mis], collapse = "', '"),
            "' not found."
        )
    }
    effectClass <- sapply(effects, function(e) class(colData(object)[, e]))
    if (any(notCateg <- !effectClass %in% c("factor", "character"))) {
        stop(
            "Provide 'contrasts' only for categorical variables. ",
            "Problematic variable(s): ",
            paste(effects[notCateg], collapse = ", ", ".")
        )
    }
    for (i in seq_along(contrasts)) {
        x <- contrasts[[i]]
        levs <- x[2:3]
        if (any(mis <- !levs %in% colData(object)[, x[1]])) {
            stop(
                "Level(s) not found for effect '", x[1], "': ",
                paste(levs[mis], collapse = ", "), "."

            )
        }
        constrastn <- make.names(paste0(x[1], "_", x[2], "_vs_", x[3]))
        names(contrasts)[i] <- constrastn
    }
    contrasts[unique(names(contrasts))]
}

## Internal function that takes a contrast and retrieves the estimated
## log fold change and associated standard error between the two
## defined groups. The inference is performed on the model
## coefficients and the variance covariance matrix that are
## automatically retrieved from the SummarizedExperimeent object. The
## log fold change and the standard error are returned for all
## features in a list.
##
## Note that the function may return NA estimates for features where
## missing data leads to loss of one or two levels of interest
## supplied by `contrast`.
##
## @param object An object that inherits from the
##   `SummarizedExperiment` class. It must contain an estimated
##   `ScpModel` in its metadata.
## @param contrast A `character(3)` with the following elements: 1.
##   The name of a categorical variable to test; 2. The name of the
##   reference group: 3. The name of the second group to contrast
##   against the reference group. `coefficients` and `contrasts`
##   cannot be both NULL.
## @param name A `character(1)` providing the name to use to retrieve
##   the model results. When retrieving a model and `name` is missing,
##   the name of the first model found in `object` is used.
##
##   cf https://stats.stackexchange.com/a/446699
.contrastToEstimates <- function(object, contrast, name) {
    fits <- scpModelFitList(object, name, filtered = TRUE)
    out <- sapply(fits, function(fit) {
        coef <- scpModelFitCoefficients(fit)
        vcov <- scpModelFitVcov(fit)
        levs <- scpModelFitLevels(fit)[[contrast[[1]]]]
        contrastMat <- .levelsToContrastMatrix(contrast, levs)
        if (is.null(contrastMat)) return(c(logFc = NA, se = NA))
        sel <- grepl(contrast[[1]], names(coef))
        logFc <- contrastMat %*% coef[sel]
        se <- sqrt(contrastMat %*% vcov[sel, sel] %*% t(contrastMat))
        c(logFc = logFc, se = se)
    })
    list(logFc = out["logFc", ], se = out["se", ])
}

## Internal function that converts a set of levels and a desired
## contrast into a contrast matrix. The contrast matrix has 1 row
## where each column is named after the model coefficients that enable
## the quantification of the log fold change between the 2 groups of
## interest.
##
## If `levels` contains only a single level or if one of the levels
## requested in `contrast` is absent (happens when a level has been
## dropped during modelling), the function return NULL since no
## contrast matrix can be computed.
##
## @param contrasts A `character(3)` with the following elements: 1.
##   The name of a categorical variable to test; 2. The name of the
##   reference group: 3. The name of the second group to contrast
##   against the reference group. `coefficients` and `contrasts`
##   cannot be both NULL.
## @param levels A `character()` containing all possible levels
##   contained by the model variable identified by contrast[[1]].
##
.levelsToContrastMatrix <- function(contrast, levels) {
    if (length(levels) <= 1 || any(!contrast[2:3] %in% levels))
        return(NULL)
    df <- data.frame(group = factor(levels, levels = levels))
    mm <- model.matrix(
        ~ 1 + group, data = df,
        contrasts.arg = .modelContrasts(df)
    )[, -1, drop = FALSE] ## to do test drop = false
    colnames(mm) <- sub("group", contrast[[1]], colnames(mm))
    rownames(mm) <- levels
    l <- rep(0, length(levels))
    names(l) <- rownames(mm)
    l[[contrast[[2]]]] <- -1
    l[[contrast[[3]]]] <- 1
    t(l) %*% mm
}

## Internal function that performs statistical inferences for a
## coefficient of interest. The function returns
## a list of data tables as returned by .computeTTest().
##
## p-values are adjusted using the IHW approach given there is enough
## data. The weights required by IHW are assigned based on the model
## intercept, with features with heigher intercept are given more
## weight.
##
## @param object An object that inherits from the
##     `SummarizedExperiment` class. It must contain an estimated
##     `ScpModel` in its metadata.
## @param coefficients A `character()` vector with coefficient names
##     to test.
## @param name A `character(1)` providing the name to use to retrieve
##     the model results. When retrieving a model and `name` is
##     missing, the name of the first model found in `object` is used.
##
.scpDifferentialAnalysisOnCoefficient <- function(object, coefficients,
                                                  name) {
    allCoefs <- sapply(scpModelCoefficients(object, name), names)
    allCoefs <- unique(unlist(allCoefs))
    if (any(mis <- !coefficients %in% allCoefs)) {
        stop(
            "Some coefficients not found: ",
            paste0(coefficients[mis], collapse = ", "), ".")
    }
    df <- scpModelDf(object, name)
    weights <- scpModelIntercept(object, name)
    out <- lapply(coefficients, function(coefficient) {
        est <- .coefficientToEstimates(object, coefficient, name)
        .computeTTest(est$beta, est$se, df, weights)
    })
    names(out) <- coefficients
    out
}

## Internal function that takes a coefficient name and retrieves the
## corresponding coefficient estimation and associated standard error,
## based on the variance covariance matrix. These information are
## automatically retrieved from the SummarizedExperimeent object. The
## log fold change and the standard error are returned for all
## features in a list.
##
## @param object An object that inherits from the
##     `SummarizedExperiment` class. It must contain an estimated
##     `ScpModel` in its metadata.
## @param coefficient A `character(1)` with the coefficient name to
##     test.
## @param name A `character(1)` providing the name to use to retrieve
##     the model results. When retrieving a model and `name` is
##     missing, the name of the first model found in `object` is used.
##
.coefficientToEstimates <- function(object, coefficient, name) {
    fits <- scpModelFitList(object, name, filtered = TRUE)
    out <- sapply(fits, function(fit) {
        if (!coefficient %in% names(scpModelFitCoefficients(fit))) {
            return(c(beta = NA, se = NA))
        }
        beta <- scpModelFitCoefficients(fit)[coefficient]
        se <- sqrt(scpModelFitVcov(fit)[coefficient, coefficient])
        c(beta = unname(beta), se = unname(se))
    })
    list(beta = out["beta", ], se = out["se", ])
}

## Internal function that computes significance of estimated
## coefficients/fold changes using t-tests. The inptu and the
## computed p-values (raw and adjusted) are returned in a DataFrame.
##
##' @importFrom IHW ihw adj_pvalues
##' @importFrom stats p.adjust pt
##' @importFrom S4Vectors DataFrame
.computeTTest <- function(beta, se, df, weights = NULL) {
    tstat <- beta / se
    pval <- 2 * pt(abs(tstat), df, lower.tail = FALSE)
    padj <- .adjustPvalue(pval, weights)
    DataFrame(
        feature = names(beta), Estimate = beta, SE = se, Df = df,
        tstatistic = tstat, pvalue = pval, padj = padj
    )
}

## Internal function to perform multiple adjustment for p-values. By
## default, the function perfors Benjamini-Hochberg correction, but
## when weigths are provided with enough features (>=3000), the
## function adjust using a weighted BH approach from the IHW method.
##
## Note the 3000 cutoff for IHW was defined after looking at the
## source code of `IHW:::ihw.default()`. When `bins = "auto"`, the
## number of bins is set to 1 (= BH correction) when the number of
## p-values is < 3000.
.adjustPvalue <- function(pval, weights = NULL) {
    if (is.null(weights) || sum(!is.na(pval)) < 3000) {
        padj <- p.adjust(pval, method = "BH")
    } else {
        res <- ihw(
            pvalues = pval, covariates = weights, alpha = 0.05,
            adjustment_type = "BH"
        )
        padj <- adj_pvalues(res)
        names(padj) <- names(pval)
    }
    padj
}


## ---- Aggregation functions ----

##' @name ScpModel-DifferentialAnalysis
##'
##' @param differentialList A list of tables returned by
##'     `scpDifferentialAnalysis()`.
##'
##' @param fcol A `character(1)` indicating the column to use for
##'     grouping features. Typically, this would be protein or gene
##'     names for grouping proteins.
##'
##' @param ... Further arguments passed to
##'     [metapod::combineGroupedPValues()].
##'
##' @importFrom metapod combineGroupedPValues
##' @importFrom QFeatures reduceDataFrame
##'
##' @export
##'
scpDifferentialAggregate <- function(differentialList, fcol, ...) {
    stopifnot(fcol %in% colnames(differentialList[[1]]))
    endoapply(differentialList, function(x) {
        x <- x[!is.na(x[[fcol]]), ]
        groupedPvals <- combineGroupedPValues(x$pvalue, x[[fcol]], ...)
        DataFrame(
            feature = sort(unique(x[[fcol]])),
            Estimate = x$Estimate[groupedPvals$representative],
            pvalue = groupedPvals$p.value,
            padj = .adjustPvalue(groupedPvals$p.value, weights = NULL),
            reduceDataFrame(x, x[[fcol]], count = TRUE, drop = TRUE)
        )
    })
}

## ---- Plotting functions ----

##' @name ScpModel-DifferentialAnalysis
##'
##' @param differentialList A list of tables returned by
##'     `scpDifferentialAnalysis()`.
##'
##' @param fdrLine A `numeric(1)` indicating the FDR threshold bar to
##'     show on the plot.
##'
##' @param top A `numeric(1)` indicating how many features should be
##'     labelled on the plot.
##'
##' @param by A `character(1)` used to order the features It
##'     indicates which variable should be considered when sorting the
##'     results. Can be one of: "Estimate", "SE", "Df", "tstatistic",
##'     "pvalue", "padj" or any other annotation added by the user.
##'
##' @param decreasing A `logical(1)` indicating whether the features
##'     should be ordered decreasingly (`TRUE`, default) or
##'     increasingly (`FALSE`) depending on the value provided by
##'     `by`.
##'
##' @param textBy A `character(1)` indicating the name of the column
##'     to use to label points.
##'
##' @param pointParams A `list` where each element is an argument that
##'     is provided to [ggplot2::geom_point()]. This is useful to
##'     change point size, transparency, or assign colour based on an
##'     annotation (see [ggplot2::aes()]).
##'
##' @param labelParams A `list` where each element is an argument that
##'     is provided to [ggrepel::geom_label_repel()]. This is useful
##'     to change label size, transparency, or assign
##'     colour based on an annotation (see [ggplot2::aes()]).
##'
##' @importFrom ggrepel geom_text_repel
##'
##' @export
scpVolcanoPlot <- function(differentialList,
                           fdrLine = 0.05,
                           top = 10, by = "padj",
                           decreasing = FALSE,
                           textBy = "feature",
                           pointParams = list(),
                           labelParams = list()) {
    stopifnot(inherits(differentialList, "List"))
    if (!textBy %in% colnames(differentialList[[1]]))
        stop("'", textBy, "' not found in results. Use ",
             "scpAnnotateResults() to add custom annotations.")
    pl <- lapply(names(differentialList), function(i) {
        contrast <- metadata(differentialList[[i]])$contrast
        df <- data.frame(differentialList[[i]])
        labelParams$data <- .filterDifferentialData(df, top, by, decreasing)
        .plotVolcano(df, pointParams, labelParams, textBy, contrast, fdrLine) +
            ggtitle(i)
    })
    names(pl) <- names(differentialList)
    pl
}

## Internal function that filters the result table as returned in the
## List by scpDifferentialAnalysis().
##
## Filtering proceeds as follows:
## 1. Focus the filter based on the value of the 'by' column
## 2. Define the best values as the highest (decreasing = TRUE) or as
##    the lowest (decreasing = FALSE).
## 3. Keep the 'top' best values.
##
## The functions returns a data.frame with filtered rows.
## If top = 0, the function returns an empty data.frame. If top is
## greater or equal to the number of rows of df, the unfiltered table
## is returned.
##
## @param df A data.frame containing inference results that need to be
##     filtered.
## @param top A `numeric(1)` indicating how many features should be
##     labelled on the plot.
## @param by A `character(1)` used to order the features It
##     indicates which variable should be considered when sorting the
##     results. Can be one of: "Estimate", "SE", "Df", "tstatistic",
##     "pvalue", "padj" or any other annotation added by the user.
## @param decreasing A `logical(1)` indicating whether the features
##     should be ordered decreasingly (`TRUE`, default) or
##     increasingly (`FALSE`) depending on the value provided by
##     `by`.
##
.filterDifferentialData <- function(df, top, by, decreasing) {
    if (!by %in% colnames(df))
        stop("'", by, "' not found in differentialList tables.")
    if (top == 0) return(df[0, , drop = FALSE])
    if (top >= nrow(df)) return(df)
    topIndex <- order(df[[by]], decreasing = decreasing)[seq_len(top)]
    df[topIndex, ]
}

## Internal function that generates a ggplot based on an input
## data.frame. The plot visualized the -log10 adjusted p-value as a
## function of the contrast/coefficient estimate. An FDR line is shown
## as well. If contrast is provided, the function will adapt the x
## axis annotation to explicit the direction of the fld change. The
## function returns a ggplot object that creates the visualisation.
##
## @param df A data.frame containing inference results that need to be
##     plotted.
## @param pointParams A `list` where each element is an argument that
##     is provided to [ggplot2::geom_point()]. This is useful to
##     change point size, transparency, or assign colour based on an
##     annotation (see [ggplot2::aes()]).
## @param labelParams A `list` where each element is an argument that
##     is provided to [ggrepel::geom_label_repel()]. This is useful
##     to change label size, transparency, or assign
##     colour based on an annotation (see [ggplot2::aes()]).
## @param textBy A `character(1)` indicating the name of the column
##     to use to label points.
## @param contrast A `character(3)` with the following elements: 1. The
##     name of a categorical variable to test; 2. The name of the
##     reference group: 3. The name of the second group to contrast
##     against the reference group. `coefficients` and `contrasts`
##     cannot be both NULL.
## @param fdrLine A `numeric(1)` indicating the FDR threshold bar to
##     show on the plot.
##
.plotVolcano <- function(df, pointParams, labelParams,
                         textBy, contrast = NULL, fdrLine = 0.05) {
    pl <- ggplot(df) +
        aes(x = .data$Estimate,
            y = -log10(.data$padj),
            label = .data[[textBy]]) +
        do.call(geom_point, pointParams) +
        geom_hline(yintercept = -log10(fdrLine)) +
        do.call(geom_text_repel, labelParams) +
        ylab("-log10(Adjusted p-value)") +
        theme_minimal()
    if (!is.null(contrast)) {
        pl <- pl + .annotateDirection(df$Estimate, contrast)
    }
    pl
}

## Internal function that generates an intuitively annotates the log
## fold change axis. This is performed by highlighting the direction
## of change for the two groups in the contrast. Note however that the
## annotation is adapted in case that the changes are unidirectional.
## Then, only the contrast group where intensities are estimated to be
## higher is shown.
##
## The function returns an object of class "labels" that can be used
## as the x axis for a ggplot fig.
##
## @param logFoldChange A numeric() containing estimated log fold
##   changes.
## @param contrast A `character(3)` with the following elements: 1.
##   The name of a categorical variable to test; 2. The name of the
##   reference group: 3. The name of the second group to contrast
##   against the reference group. `coefficients` and `contrasts`
##   cannot be both NULL.
.annotateDirection <- function(logFoldChange, contrast) {
    xAnnotation <- "log2(Fold change)"
    if (any(logFoldChange > 0))
        xAnnotation <- paste(xAnnotation, "  ->", contrast[[3]])
    if (any(logFoldChange < 0))
        xAnnotation <- paste(contrast[[2]], "<-  ", xAnnotation)
    xlab(xAnnotation)
}
