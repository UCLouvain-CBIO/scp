
## ---- ScpModel class definition ----

##' @rdname ScpModel-class
##'
##' @title Class to store the results of single-cell proteomics
##' modelling
##'
##' @description
##'
##' An `ScpModel` object must be always stored in the `metadata()` of
##' an object that inherits from the `SingleCellExperiment` class. The
##' `ScpModel` object should **never be accessed directly** by the
##' user. Instead, we provide several setter function to retrieve
##' information that may be useful to the user.The `ScpModel` class
##' contains several slots:
##'
##' - `scpModelFormula`: a `formula` object controlling which
##'   variables are to be modelled.
##' - `scpModelInputIndex`: a `numeric(1)`, selecting the assay to use
##'   in the `SingleCellExperiment` object as input matrix. Note that
##'   this slot serves as a pointer, meaning that the quantitative
##'   data is not duplicated. Any change to the assay in the
##'   `SingleCellExperiment` will impact the estimation of the
##'   `ScpModel` object.
##' - `scpModelFilterThreshold`: A `numeric(1)` indicating the minimal
##'   n/p ratio required for a feature to be included in further model
##'   exploration. n/p is the number of measured values for a features
##'   divided by the number of coefficients to estimate. n/p cannot be
##'   smaller than 1 because this would lead to over-specified models.
##' - `scpModelFitList`: A `List` that contains the model results for
##'   each feature. Each element is a `ScpModelFit` object (see
##'   [`ScpModelFit-class`])
##'
##' @section Getters:
##'
##' Each slot has a getter function associated:
##'
##' - `scpModelNames()`: returns a vector of names of `ScpModel`
##'   objects stored in the `SingleCellExperiment` object.
##' - `scpModelFormula()`: returns the `formula` slot of the `ScpModel`
##'   within an object that inherits from the `SummarizedExperiment`
##'   class.
##' - `scpModelFilterThreshold()`: returns the n/p ration threshold
##'   used for feature filtering.
##' - `scpModelInput()`: returns a `matrix` with the quantitative
##'   values used as input of the model. Hence, the matrix contains
##'   the data before modelling. If `filtered = TRUE`, the feature of
##'   the matrix are restricted to the features that satisfy the n/p
##'   ratio threshold.
##' - `scpModelFilterNPRatio()`: returns the computed n/p ratio for
##'   each feature. If `filtered = TRUE`, the function returns only
##'   the n/p of the features that satisfy the n/p ratio threshold.
##' - `scpModelResiduals()`: when `join = FALSE`, the function returns
##'   a list where each element corresponds to a feature and contains
##'   the estimated residuals. When `join = TRUE` (default), the function
##'   combines the list into a matrix with features in rows and cells
##'   in columns, and filling the gaps with `NA`. If `filtered = TRUE`,
##'   the feature of the matrix are restricted to the features that
##'   satisfy the n/p ratio threshold.
##' - `scpModelEffects()`: when `join = FALSE`, the function return a
##'   list where each element of the list corresponds to a feature.
##'   Each element contains another list with as many elements as
##'   variable in the model and each element contains the data  effect vector
##'   for that vector. When `join = TRUE` (default), each element of the list is
##'   a matrix with features in rows and cells in columns where gaps
##'   are filled with `NA`. If `filtered = TRUE`, the feature of the
##'   matrix are restricted to the features that satisfy the n/p
##'   ratio threshold.
##'
##' Setter:
##'
##' - `scpModelFilterThreshold<-()`: the function changes the n/p
##'   ratio threshold used for filtering features.
##'
##' @seealso
##' - [ScpModelFit-class] for a description of the class that store
##'   modelling results
##' - [ScpModel-Workflow] that uses the class to store the estimated
##'   model.
##'
##' @author Christophe Vanderaa, Laurent Gatto
##'
##' @example inst/examples/examples_ScpModel-Class.R
##'
##' @importFrom methods new
##'
##' @name ScpModel
##' @aliases ScpModel ScpModel-class class:ScpModel
##' @exportClass ScpModel
setClass("ScpModel", slots = c(
    scpModelFormula = "formula",
    scpModelInputIndex = "numeric",
    scpModelFilterThreshold = "numeric",
    scpModelFitList = "List"
))

## Class constructors
## Keep the constructor hidden from the user. It is better they
## build this object within an S(C)E object through `scpModelWorkflow()`
ScpModel <- function() {
    ans <- new("ScpModel")
    ans@scpModelFitList <- List()
    ans
}

## ---- Exported getters ----

##' @rdname ScpModel-class
##'
##' @export
scpModelFormula <- function(object, name) {
    out <- scpModel(object, name)@scpModelFormula
    .checkModelElement(
        out, "scpModelFormula", .runWorkflowMessage,
        name = .checkModelName(object, name)
    )
    out
}

##' @rdname ScpModel-class
##'
##' @param filtered A `logical(1)` indicating whether the output
##'     should return all features (`FALSE`) or the features that
##'     comply to the n/p ratio threshold (`TRUE`).
##'
##' @importFrom SummarizedExperiment assay
##'
##' @export
scpModelInput <- function(object, name, filtered = TRUE) {
    out <- assay(object, scpModelInputIndex(object, name))
    if (filtered) {
        out <- out[scpModelFeatureNames(object, name), , drop = FALSE]
    }
    out
}

##' @rdname ScpModel-class
##'
##' @export
scpModelFilterThreshold <- function(object, name) {
    out <- scpModel(object, name)@scpModelFilterThreshold
    .checkModelElement(
        out, "scpModelFilterThreshold", .runWorkflowMessage,
        name = .checkModelName(object, name)
    )
    out
}

##' @rdname ScpModel-class
##'
##' @export
scpModelFilterNPRatio <- function(object, name, filtered = TRUE) {
    out <- scpModelN(object, name, filtered) /
        scpModelP(object, name, filtered)
    ## when all missing or no params
    out[is.na(out) | is.infinite(out)] <- 0
    out
}

##' @rdname ScpModel-class
##'
##' @param join A `logical(1)` indicating whether the output should be
##'     combined in a single matrix (`TRUE`) or it should be returned
##'     as a list with one element for each feature (`FALSE`). When
##'     `TRUE`, any gaps across features will be filled with NA's.
##' @export
scpModelResiduals <- function(object, name, join = TRUE, 
                              filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "Residuals", filtered, .runWorkflowMessage
    )
    if (join) out <- .joinScpModelOutput(out, object)
    out
}

##' @rdname ScpModel-class
##'
##' @export
scpModelEffects <- function(object, name, join = TRUE,
                            filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "Effects", filtered, .runWorkflowMessage
    )
    if (join) {
        out <- lapply(scpModelEffectNames(object, name), function(e) {
            effect <- endoapply(out, "[[", e)
            .joinScpModelOutput(effect, object)
        })
        names(out) <- scpModelEffectNames(object, name)
        out <- as(out, "List")
    }
    out
}

##' @rdname ScpModel-class
##'
##' @export
scpModelNames <- function(object) {
    errorMessage <- paste0(
        "No 'ScpModel' found in object. ", .runWorkflowMessage
    )
    if (!length(metadata(object))) {
        stop(errorMessage)
    }
    sel <- which(sapply(
        metadata(object),
        function(x) inherits(x, "ScpModel")
    ))
    if (!length(sel)) {
        stop(errorMessage)
    }
    names(metadata(object))[sel]
}

## A vector of strings with the available component analysis methods
##' @rdname ScpModel-ComponentAnalysis
##'
##' @export
scpModelComponentMethods <- c("APCA", "ASCA", "ASCA.E")

## ---- Internal getters ----

## Internal functions to get an ScpModel object or its assays from an
## SE object. This is not exported because there is no good reason for
## the users to extract an ScpModel from its container SE.
scpModel <- function(object, name) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    name <- .checkModelName(object, name)
    metadata(object)[[name]]
}

## Internal function that retrieves the scpModelInputIndex slot.
scpModelInputIndex <- function(object, name) {
    out <- scpModel(object, name)@scpModelInputIndex
    .checkModelElement(
        out, "scpModelInputIndex", .runWorkflowMessage,
        name = .checkModelName(object, name)
    )
    out
}

## Internal function that retrieves the scpModelFitList slot.
scpModelFitList <- function(object, name, filtered = FALSE) {
    out <- scpModel(object, name)@scpModelFitList
    .checkModelElement(
        out, "scpModelFitList", .runWorkflowMessage,
        name = .checkModelName(object, name)
    )
    if (filtered) out <- out[scpModelFeatureNames(object, name)]
    out
}

## Internal function that extracts a specific element from a
## scpModelFitList in an ScpModel object.
## @param object Object that inherits from SummarizedExperiment class
## @param name A `character(1)` providing the name to use to store or
##     retrieve the modelling results.
## @param what A `character(1)` providing the name of the element in
##     the ScpModelFit object to retrieve.
## @param filtered A `logical(1)` indicating whether the output should
##     be filtered based on the NP ratio.
## @param helpMessage A `character(1)` that provides additional
##     information in case the retrieved elements are empty in the
##     object.
scpModelFitElement <- function(object, name, what, filtered,
                               helpMessage = "") {
    out <- scpModelFitList(object, name)
    scpModelFitSlot <- try(get(paste0("scpModelFit", what)), TRUE)
    if (inherits(scpModelFitSlot, "try-error"))
        stop("'", what, "' is not a slot of an ScpModelFit object.")
    out <- endoapply(out, scpModelFitSlot)
    .checkModelElement(
        out, paste(what, "(in ScpModelFit)"), helpMessage,
        name = .checkModelName(object, name)
    )
    if (filtered) out <- out[scpModelFeatureNames(object, name)]
    out
}

scpModelN <- function(object, name, filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "N", filtered, .runWorkflowMessage
    )
    unlist(out)
}

scpModelP <- function(object, name, filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "P", filtered, .runWorkflowMessage
    )
    unlist(out)
}

scpModelCoefficients <- function(object, name, filtered = TRUE) {
    scpModelFitElement(
        object, name, "Coefficients", filtered, .runWorkflowMessage
    )
}

scpModelDf <- function(object, name, filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "Df", filtered, .runWorkflowMessage
    )
    unlist(out)
}

scpModelVar <- function(object, name, filtered = TRUE) {
    out <- scpModelFitElement(
        object, name, "Var", filtered, .runWorkflowMessage
    )
    unlist(out)
}

scpModelUvcov <- function(object, name, filtered = TRUE) {
    scpModelFitElement(
        object, name, "Uvcov", filtered, .runWorkflowMessage
    )
}

scpModelVcov <- function(object, name, filtered = TRUE) {
    scpModelFitElement(
        object, name, "Vcov", filtered, .runWorkflowMessage
    )
}

scpModelIntercept <- function(object, name, filtered = TRUE) {
    coefs <- scpModelCoefficients(object, name, filtered)
    sapply(coefs, function(x) x[["(Intercept)"]])
}

## Internal function that returns a vector of strings with the feature
## names (= row names) that have been kept for modelling. If no filter
## is found, it returns all feature names.
scpModelFeatureNames <- function(object, name) {
    threshold <- scpModelFilterThreshold(object, name)
    if (!length(threshold)) {
        return(rownames(object))
    }
    sel <- scpModelFilterNPRatio(object, name, FALSE) >= threshold
    rownames(object)[sel]
}

## Internal function that returns the effect names in the model
scpModelEffectNames <- function(object, name) {
    labels(terms(
        scpModelFormula(object, name),
        data = colData(object)
    ))
}

## ---- Exported setters ----

##' @rdname ScpModel-class
##'
##' @param object An object that inherits from the
##'     `SingleCellExperiment` class.
##'
##' @param value An `numeric(1)`, the new value for the n/p ratio
##'     threshold
##'
##' @param name A `character(1)` providing the name to use to store or
##'     retrieve the modelling results. When retrieving a model and
##'     `name` is missing, the name of the first model found in
##'     `object` is used.
##'
##' @export
`scpModelFilterThreshold<-` <- function(object, name, value) {
    stopifnot(length(value) == 1)
    stopifnot(value >= 1)
    scpModel(object, name)@scpModelFilterThreshold <- value
    object
}

## ---- Internal setters ----

## All setters are internal function and hence not visible to the
## user, otherwise it may lead to corrupt model results.
## We implement setters instead of initializing the slots when
## constructing the object because the validity depends on the SE
## object that is outside of the ScpModel object.

`scpModel<-` <- function(object, name, value) {
    stopifnot(inherits(value, "ScpModel") | is.null(value))
    if (missing(name)) name <- .defaultModelName(object)
    metadata(object)[[name]] <- value
    object
}

`scpModelFormula<-` <- function(object, name, value) {
    stopifnot(inherits(value, "formula"))
    value <- .checkScpModelFormula(value, object)
    scpModel(object, name)@scpModelFormula <- value
    object
}

##' @importFrom SummarizedExperiment assayNames
`scpModelInputIndex<-` <- function(object, name, value) {
    stopifnot(!is.null(colnames(object)))
    stopifnot(!is.null(rownames(object)))
    value <- .checkInputIndex(value, assayNames(object), "i")
    x <- assay(object, value)
    if (any(is.infinite(x)))
        stop("The selected assay ('assay(object, i)') contains infinite values.")
    scpModel(object, name)@scpModelInputIndex <- value
    object
}

`scpModelFitList<-` <- function(object, name, value) {
    stopifnot(inherits(value, "List"))
    stopifnot(all(sapply(value, inherits, "ScpModelFit")))
    stopifnot(identical(rownames(object), names(value)))
    scpModel(object, name)@scpModelFitList <- value
    object
}

## ---- Internal utility functions ----

## A string suggesting how to run the scplainer workflow the model
.runWorkflowMessage <-
    "Use 'scpModelWorkflow()' to run the scplainer modelling workflow."

## Internal function that returns the name of the default model.
## The default model is the first ScpModel object in the metadata
## Note: creating this function may seem overkill, but I add it
## in case we come up with another rationale for selecting the default
## model.
.defaultModelName <- function(object) {
    scpModelNames(object)[[1]]
}

## Internal function that checks whether `name` points to a valid
## ScpModel in `object`. If `name` is not a single value or if it does
## not point to a ScpModel object in `metadata(object)`, the function
## throws an error.
## @param object A SummarizedExperiment object
## @param name A character(1) to select a model. When missing, name is
##      assigned as the name of the default model
.checkModelName <- function(object, name) {
    if (missing(name)) name <- .defaultModelName(object)
    stopifnot(length(name) == 1)
    if (!name %in% scpModelNames(object)) {
        stop(
            "Model name '", name, "' not found in object. ",
            .runWorkflowMessage
        )
    }
    name
}

## Internal functions that checks whether the provided component is
## empty or not. This component is expected to be extracted from one
## of the slot of an ScpModel object contained in a
## SummarizedExperiment. If `x` is empty, it throws a meaningful error
## using `name`, `what` and `additionalMessage`.
## @param x An element to check, can be atomic or a list
## @param name The name of the model from which x is retrieved.
## @param what The name of the slot in the model from which x is
##     retrieved.
## @param additionalMessage A string that provides additional
##     information, typically indicating how the empty element should
##     be created.
.checkModelElement <- function(x, name, what, additionalMessage = "") {
    if (inherits(x, "List") || is.list(x)) {
        isEmpty <- length(x) == 0 || sum(sapply(x, length)) == 0
    } else {
        isEmpty <- length(x) == 0
    }
    if (isEmpty) {
        stop(
            "No available ", what, " for model '", name, "'. ",
            additionalMessage
        )
    }
    NULL
}

## Internal functions that checks the formula is valid for modelling.
## If a response variable is provided it is removed automatically and
## if an intercept term is missing, it is added automatically. When the
## formula contains `.`, it is replaced by all variables present in
## `colData(object)`, except for variable already present in the
## formula. The function returns the cleaned formula.
##
## @param formula A formula
## @param object An object that inherits from SummarizedExperiment
##
##' @importFrom stats as.formula
.checkScpModelFormula <- function(formula, object) {
    fterms <- terms(formula, data = colData(object))
    formula <- .removeResponseVariables(formula, fterms)
    formula <- .checkExplanatoryVariables(
        formula, fterms, colnames(colData(object))
    )
    formula
}

## Internal function that removes the response variable from formula
## terms. If the formula contains a response variable, the function
## throws a warning
## @param formula A formula.
## @param fterms A terms object derived from formula.
.removeResponseVariables <- function(formula, fterms) {
    if (!identical(attr(fterms, "response"), 0L)) {
        warning(
            "The formula contains a response variable and is ignored."
        )
        formula <- as.formula(
            paste("~", as.character(formula)[[3]]),
            env = attr(fterms, ".Environment")
        )
    }
    formula
}

## Internal functions that checks the variables in a formula are
## available from the sample annotations contained in the colData of
## an object that inherits from SummarizedExperiment. The function
## returns a cleaned formula.
## - The formula contains a no intercept = warning
## - "Residuals" cannot be a variable name = error
## - The formula has no variables (excluding intercept) = error
## - The colData is empty = error
## - The colData is missing some variables from the formula = error
## @param formula A formula.
## @param fterms A terms object derived from formula.
## @param availableVariables A vector of available variable names to
##     model.
.checkExplanatoryVariables <- function(formula, fterms,
                                       availableVariables) {
    modelVars <- all.vars(formula)
    modelVars <- .replaceDotVariable(modelVars, availableVariables)
    if (attr(fterms, "intercept") != 1)
        warning("No intercept in the formula. It is added automatically.")
    if ("Residuals" %in% modelVars)
        stop("'Residuals' is reserved. Please rename that variable.")
    if (!length(modelVars)) {
        stop("You provided a formula with no variable to model.")
    } else {
        if (!length(availableVariables))
            stop("colData(object) is empty.")
        if (any(mis <- !modelVars %in% availableVariables))
            stop("colData(object) is missing one or more variables ",
                 "from the formula: ",
                 paste(modelVars[mis], collapse = ", "), ".")
    }
    as.formula(
        paste(c("~ 1", modelVars), collapse = " + "),
        env = attr(formula, ".Environment")
    )
}

## Internal function that take model variables and replace the "."
## shorthand with all remaining available variables.
## @param modelVars A vector of model variable names extracted
##     from a formula.
## @param availableVariables A vector of available variable names. If
##     'modelVars' contains a '.', it will be replaced by all
##     names in 'availableVariables' not present in 'modelVars'.
.replaceDotVariable <- function(modelVars, availableVariables) {
    if (any(modelVars == ".")) {
        modelVars <- unique(c(
            modelVars, availableVariables
        ))
        modelVars <- modelVars[modelVars != "."]
    }
    modelVars
}


## Internal function that converts and checks a provided assay or model
## index given a set of choices. The function return a single numeric
## corresponding to a valid assay or model, otherwise it throws an error.
## @param index The index to convert and check. The index can be
##     logical, numeric, or character.
## @param choices The names of valid assays or models that the index
##     is pointing to.
## @param what A character string that indicates the name of the
##     index. It is used to provide meaningful error messages.
.checkInputIndex <- function(index, choices, what = "name") {
    if (is.logical(index) && !identical(index, NA)) {
        index <- which(index)
    }
    if (length(index) != 1) {
        stop("'", what, "' points to multiple input assays.")
    }
    if (is.character(index)) {
        if (!index %in% choices) stop("'", index, "' not found.")
        index <- which(choices == index)
    } else if (is.numeric(index) && !inherits(index, "matrix")) {
        if (length(choices) && index > length(choices)) {
            stop("'", what, "' is out of bounds")
        }
    } else {
        stop("'", what, "' must be a character, numeric or logical")
    }
    index
}

## Internal function that combines a list of model output elements
## into a matrix with as many rows as fitted models and as many
## columns a the number of samples in object. Any piece of information
## missing for a sample in a model is filled with NA.
## Typically, these elements are one of the following slots from
## multiple ScpModelList objects: residuals, effect matrices.
## @param x A list of model output to combine. Each element contains
##     a named vector. The names should relate to names in object
## @param object A SummarizedExperiment
.joinScpModelOutput <- function(x, object) {
    stopifnot(length(names(x)) > 0)
    out <- matrix(NA, nrow = length(x), ncol = ncol(object),
                  dimnames = list(names(x), colnames(object)))
    for (i in seq_along(x)) {
        cols <- names(x[[i]])
        out[i, cols] <- x[[i]]
    }
    out
}
