## ---- ScpModelFit class definition ----

##' @rdname ScpModelFit-class
##'
##' @title Class to store the components of an estimated model for a
##' feature
##'
##' @description
##'
##' An `ScpModelFit` object is expected to be stored as a list element
##' in the `scpModelFitList` of an `ScpModel` object. The
##' `ScpModelFit` object should **never be accessed directly** by the
##' user. Refer to the [ScpModel-class] for a list of function to
##' access the information in an `ScpModelFit`. The `ScpModelFit`
##' class contains several slots that contain the model output for a
##' feature:
##'
##' - `n`: an `integer`, the number of observations for the feature
##' - `p`: an `integer`, the number of coefficient to estimate
##' - `coefficients`: a `numeric` vector with the estimated
##'   coefficients
##' - `residuals`: a `numeric` vector with the estimated residuals
##' - `effects`: a `List` with the
##' - `df`: an `integer` providing the number of degrees of freedom
##'   of the model estimation
##' - `var`: a `numeric` vector with the residual variance of the
##'   model estimation
##' - `uvcov`: the unscaled variance covariance `matrix`
##' - `levels`: a named `List` where each elements corresponds to a
##'   categorical model variable and contains a vector with the
##'   possible categories.
##'
##' @seealso
##' [ScpModel-class] for a description of the class that relies on
##' `ScpModelFit`
##'
##' @author Christophe Vanderaa, Laurent Gatto
##'
##' @examples
##' new("ScpModelFit") ## this should never be used by the user
##'
##' @name ScpModelFit
##' @aliases ScpModelFit ScpModelFit-class class:ScpModelFit
##' @exportClass ScpModelFit
setClass("ScpModelFit", slots = c(
    n = "integer",
    p = "integer",
    coefficients = "numeric",
    residuals = "numeric",
    effects = "List",
    df = "numeric",
    var = "numeric",
    uvcov = "matrix",
    levels = "List"
))

## Class constructors
ScpModelFit <- function(n, p) {
    stopifnot(n >= 0)
    stopifnot(p >= 0)
    stopifnot(length(n) == 1)
    stopifnot(length(p) == 1)
    new(
        "ScpModelFit", n = n, p = p, effects = List(),
        levels = List()
    )
}

## ---- Getters ----

scpModelFitN <- function(object) {
    object@n
}

scpModelFitP <- function(object) {
    object@p
}

scpModelFitCoefficients <- function(object) {
    object@coefficients
}

scpModelFitResiduals <- function(object) {
    object@residuals
}

scpModelFitEffects <- function(object) {
    object@effects
}

scpModelFitDf <- function(object) {
    object@df
}

scpModelFitVar <- function(object) {
    object@var
}

scpModelFitUvcov <- function(object) {
    object@uvcov
}

scpModelFitVcov <- function(object) {
    scpModelFitVar(object) * scpModelFitUvcov(object)
}

scpModelFitLevels <- function(object) {
    object@levels
}

## ---- Setters ----

`scpModelFitCoefficients<-` <- function(object, value) {
    stopifnot(length(value) == scpModelFitP(object))
    object@coefficients <- value
    object
}

`scpModelFitResiduals<-` <- function(object, value) {
    stopifnot(length(value) == scpModelFitN(object))
    object@residuals <- value
    object
}

`scpModelFitEffects<-` <- function(object, value) {
    stopifnot(all(sapply(value, function(x) inherits(x, "numeric"))))
    enames <- lapply(value, names)
    if (length(enames) == 0) {
        refNames <- c()
    } else {
        refNames <- enames[[1]]
    }
    if (!all(sapply(enames, function(x) identical(x, refNames))))
        stop("Effect vectors do not share identical names.")
    if (!identical(names(scpModelFitResiduals(object)), refNames))
        stop("Effects and residuals do not share identical names.")
    object@effects <- value
    object
}

`scpModelFitDf<-` <- function(object, value) {
    stopifnot(length(value) == 1)
    stopifnot(value >= 0)
    object@df <- value
    object
}

`scpModelFitVar<-` <- function(object, value) {
    stopifnot(length(value) == 1)
    stopifnot(value >= 0)
    object@var <- value
    object
}

`scpModelFitUvcov<-` <- function(object, value) {
    coef <- scpModelFitCoefficients(object)
    stopifnot(identical(names(coef), rownames(value)))
    stopifnot(identical(names(coef), colnames(value)))
    object@uvcov <- value
    object
}

## Internal function to set the ScpModelFit levels. The value must be
## a List where each element must be a vector of character. An empty
## List is allowed. When value is not null, the names of the values
## must match some or all of the effect names.
`scpModelFitLevels<-` <- function(object, value) {
    stopifnot(all(sapply(value, function(x) inherits(x, "character"))))
    if (length(value) > 0) {
        if (is.null(names(value)))
            stop("List of levels must be named.")
        if (any(!names(value) %in% names(scpModelFitEffects(object))))
            stop("Some levels are not matched to effects.")
    }
    object@levels <- value
    object
}
