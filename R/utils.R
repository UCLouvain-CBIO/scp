##' Extract the `rowData` of a `QFeatures` object to a 
##' `DataFrame`
##'
##' The methods takes the `rowData` of one or more given assay in a
##' `QFeatures` object and combines the data in a single `DataFrame`.
##' 
##' Along with the reuired `rowData` an additional `.assay` variable
##' is created and holds the assay name from which the metadata was
##' taken.
##' 
##' @param obj A `QFeatures` object
##'
##' @param i A `numeric()` or `character()` vector indicating from
##'     which assays the `rowData` should be taken.
##' 
##' @param vars A `character()` vector indicating which variables from
##'     the `rowData` should be extracted.
##'
##' @return A `DataFrame` object with the `rowData` row-binded over
##'     the required assays.
##'     
##' @export
##'
##' @importFrom MultiAssayExperiment ExperimentList experiments
##'
##' @examples
##' ## Extract the peptide length and sequence from the first 3 assays
##' data("scp1")
##' rowDataToDF(scp1, i = 1:3, c("Length", "Sequence"))
##' 
rowDataToDF <- function(obj, 
                        i, 
                        vars) {
    if (!inherits(obj, "QFeatures")) stop("'obj' must be a QFeatures object")
    if (is.numeric(i)) i <- names(obj)[i]
    ## Make sure that the variables to extract are present in the rowData
    mis <- vapply(experiments(obj)[, , i], 
                  function(x) any(!vars %in% colnames(rowData(x))),
                  logical(1))
    if (any(mis)) 
        stop("rowData variable(s) not found in:\n", 
             paste(i[mis], collapse = ", "))
    ## Extract the rowData and add from which assay it was extracted
    out <- lapply(i, function(ii) {
        x <- rowData(obj[[ii]])[, vars, drop = FALSE]
        cbind(x, .assay = ii, .rowname = rownames(x))
    })
    do.call(rbind, out)
}


##' Transfer the `colData` to an Assay
##' 
##' The function transfers the `colData` from a `QFeatures` object to
##' one of the assays it contains. The transfered data is bound to the
##' existing `colData` of the target assay.
##'
##' @param obj A `QFeatures` object
##' 
##' @param i A `numeric(1)` or `character(1)` indicating which assay
##'     to transfer the `colData` to.
##'
##' @return A `QFeatures` object
##' 
##' @export
##'
##' @examples
##' data("scp1")
##' colData(scp1[["peptides"]])
##' scp1 <- transferColDataToAssay(scp1, i = "peptides")
##' colData(scp1[["peptides"]])
##' 
transferColDataToAssay <- function (obj, 
                                    i) {
    cd <- colData(obj)[colnames(obj[[i]]),  ]
    if (all(colnames(cd) %in% colnames(colData(obj@ExperimentList[[i]])))) {
        message("The colData is already present in assay '", i, "'.")
        return(obj)
    }
    colData(obj@ExperimentList[[i]]) <- cbind(colData(obj[[i]]), cd)
    return(obj)
}

##' Aggregate features over multiple assays
##' 
##' This function is a wrapper function around
##' [QFeatures::aggregateFeatures]. 
##' It allows the user to provide multiple assays for which 
##' `aggregateFeatures` will be applied sequentially.
##' 
##' @param obj A `QFeatures` object
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
aggregateFeaturesOverAssays <- function(obj, 
                                        i, 
                                        fcol, 
                                        name, 
                                        fun, 
                                        ...) {
    if (length(i) != length(name)) stop("'i' and 'name' must have same length")
    if (length(fcol) == 1) fcol <- rep(fcol, length(i))
    if (length(i) != length(fcol)) stop("'i' and 'fcol' must have same length")
    if (is.numeric(i)) i <- names(obj)[i]
    
    ## Compute the aggregated assays
    el <- experiments(obj)[, , i]
    for (j in seq_along(el)) {
        suppressMessages(
            el[[j]] <- aggregateFeatures(el[[j]], 
                                         fcol = fcol[j], 
                                         fun = fun, 
                                         ...)
        )
        ## Print progress
        message("\rAggregated: ", j, "/", length(el), "\n")
    }
    names(el) <- name
    ## Get the AssayLinks for the aggregated assays 
    alnks <- lapply(seq_along(i), function(j) {
        hits <- QFeatures:::.get_Hits(rdFrom = rowData(obj[[i[j]]]),
                                      rdTo = rowData(el[[j]]), 
                                      varFrom = fcol[[j]], 
                                      varTo = fcol[[j]])
        AssayLink(name = name[j], from = i[j], fcol = fcol[j], hits = hits)
    })
    ## Append the aggregated assays and AssayLinks to the previous assays
    el <- c(obj@ExperimentList, el)
    alnks <- append(obj@assayLinks, AssayLinks(alnks))
    ## Update the sampleMapfrom the data 
    smap <- MultiAssayExperiment:::.sampleMapFromData(colData(obj), el)
    
    ## Create the new QFeatures object
    new("QFeatures",
        ExperimentList = el,
        colData = colData(obj),
        sampleMap = smap,
        metadata = metadata(obj),
        assayLinks = alnks)
}

