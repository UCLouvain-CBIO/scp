
##' Aggregate features over multiple assays
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
##' @param uniformRowData A `logical()` that specify if the new assays should
##'     have the same rowData invariant columns
##'     (Note that turning uniformRowData to FALSE can have a negative impact
##'     on perfomances)
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
aggregateFeaturesOverAssays <- function(object,
                                        i,
                                        fcol,
                                        name,
                                        fun,
                                        uniformRowData = TRUE,
                                        ...) {
    if (length(i) != length(name)) stop("'i' and 'name' must have same length")
    if (length(fcol) == 1) fcol <- rep(fcol, length(i))
    if (length(i) != length(fcol)) stop("'i' and 'fcol' must have same length")
    if (is.numeric(i)) i <- names(object)[i]

    ## Compute the aggregated assays
    el <- experiments(object)[i]
    if (uniformRowData) rowDataColsKept <- colnames(rowData(el[[1]]))
    for (j in seq_along(el)) {
        if (uniformRowData) rowData(el[[j]]) <- rowData(el[[j]])[, rowDataColsKept]
        suppressMessages(
            el[[j]] <- aggregateFeatures(el[[j]], fcol = fcol[j],
                                         fun = fun, ...)
        )
        if (uniformRowData) rowDataColsKept <- intersect(rowDataColsKept, colnames(rowData(el[[j]])))
        ## Print progress
        message("\rAggregated: ", j, "/", length(el), "\n")
    }
    if (uniformRowData) {
        for (j in seq_along(el)) {
        rowData(el[[j]]) <- rowData(el[[j]])[, rowDataColsKept]
        }
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

