##' @name aggregateFeaturesOverAssays-deprecated
##'
##' @title Aggregate features over multiple assays
##'
##' @description
##' The `aggregateFeaturesOverAssays` function is deprecated and will be
##' removed in a future release. Please use the `aggregateFeatures` method
##' from the `QFeatures` package instead.
##'
##' This function is a wrapper function around
##' [QFeatures::aggregateFeatures].
##' It allows the user to provide multiple assays for which
##' `aggregateFeatures` will be applied sequentially.
##'
##' @param object A `QFeatures` object
##'
##' @param i A `numeric(1)` or `character(1)` indicating which assay
##'     to transfer the `colData` to.
##'
##' @param fcol The feature variables for each assays `i` defining how
##'     to summarise the QFeatures. If `fcol` has length 1, the
##'     variable name is assumed to be the same for all assays
##'
##' @param name A `character()` naming the new assay. `name` must have
##'     the same length as `i`. Note that the function will fail if of
##'     the names in `name` is already present.
##'
##' @param fun A function used for quantitative feature aggregation.
##'
##' @param ... Additional parameters passed the `fun`.
##'
##' @return A `QFeatures` object
##'
##' @export
##'
##' @importFrom utils flush.console
##' @importFrom methods new
##' @importFrom S4Vectors metadata
##' @importFrom MultiAssayExperiment experiments
##'
##' @seealso [QFeatures::aggregateFeatures]
##' @aliases aggregateFeaturesOverAssays
##'
##' @examples
##' data("scp1")
##' scp1 <- aggregateFeaturesOverAssays(scp1,
##'                                     i = 1:3,
##'                                     fcol = "peptide",
##'                                     name = paste0("peptides", 1:3),
##'                                     fun = colMeans,
##'                                     na.rm = TRUE)
##' scp1
##'
aggregateFeaturesOverAssays <- function(object, i, fcol, name, fun, ...) {
    .Deprecated("aggregateFeatures", package = "QFeatures")
    if (length(i) != length(name)) stop("'i' and 'name' must have same length")
    if (length(fcol) == 1) fcol <- rep(fcol, length(i))
    if (length(i) != length(fcol)) stop("'i' and 'fcol' must have same length")
    if (is.numeric(i)) i <- names(object)[i]

    ## Compute the aggregated assays
    el <- experiments(object)[i]
    for (j in seq_along(el)) {
        suppressMessages(
            el[[j]] <- aggregateFeatures(el[[j]], fcol = fcol[j],
                                         fun = fun, ...)
        )
        ## Print progress
        message("\rAggregated: ", j, "/", length(el), "\n")
    }
    names(el) <- name
    ## Get the AssayLinks for the aggregated assays
    alnks <- lapply(seq_along(i), function(j) {
        hits <- QFeatures:::.get_Hits(rdFrom = rowData(object[[i[j]]]),
                                      rdTo = rowData(el[[j]]),
                                      varFrom = fcol[[j]],
                                      varTo = fcol[[j]])
        AssayLink(name = name[j], from = i[j], fcol = fcol[j], hits = hits)
    })
    ## Append the aggregated assays and AssayLinks to the previous assays
    el <- c(object@ExperimentList, el)
    alnks <- append(object@assayLinks, AssayLinks(alnks))
    ## Update the sampleMapfrom the data
    smap <- MultiAssayExperiment:::.sampleMapFromData(colData(object), el)
    ## Create the new QFeatures object
    new("QFeatures",
        ExperimentList = el,
        colData = colData(object),
        sampleMap = smap,
        metadata = metadata(object),
        assayLinks = alnks)
}
