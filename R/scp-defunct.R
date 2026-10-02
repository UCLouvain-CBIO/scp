##' @name scp-defunct
##'
##' @title Defunct Functions
##'
##' @rdname scp-deprecated
NULL

##' @rdname scp-deprecated
##' @details
##'
##' The `aggregateFeaturesOverAssays` function is defunct. Please use the
##' `QFeatures::aggregateFeatures` method from the `QFeatures` package instead.
##'
##' @usage NULL
##'
##' @export
aggregateFeaturesOverAssays <-
    function(object, i, fcol, name, fun, ...)
        .Defunct("QFeatures::aggregateFeatures", package = "scp")