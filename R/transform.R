
####---- Internal functions ----####

## Internal function that normalizes SummarizedExperiment object using 
## proteomics normalization methods available in MsCoreUtils::normalize_matrix
##' @importFrom MsCoreUtils normalize_matrix
.normalizeSCP <- function(x, method, ...) {
  if(!inherits(x, "SummarizedExperiment")) 
    stop("The assay must be a 'SummarizedExperiment' object.")
  e <- normalize_matrix(assay(x), method, ...)
  rownames(e) <- rownames(assay(x))
  colnames(e) <- colnames(assay(x))
  assay(x) <- e
  x
}


####---- Exported functions ----####


##' Divide assay columns by a reference column
##'
##' The function divides the sample columns by a reference column. The sample 
##' and reference columns are defined based on the provided `colvar` 
##' variable and on regular expression matching.
##' 
##' The supplied assay(s) are replaced with the values computed after reference
##' division.
##' 
##' @param object A `QFeatures` object
##' 
##' @param i A `numeric()` or `character()` vector indicating from which 
##'     assays the `rowData` should be taken.
##' 
##' @param colvar A `character(1)` indicating the variable to take from 
##'     `colData(object)` that gives the sample annotation.
##' 
##' @param samplePattern A `character(1)` pattern that matches the sample 
##'     encoding in `colvar`. By default all samples are devided (using the
##'     regex wildcard `.`).
##' 
##' @param refPattern A `character(1)` pattern that matches the carrier 
##'     encoding in `colvar`. Only one match per assay is allowed, otherwise
##'     only the first match is taken
##'
##' @return A `QFeatures` object
##' 
##' @export
##'
##' @examples
##' data("scp1")
##' scp1 <- divideByReference(scp1, 
##'                           i = 1, 
##'                           colvar = "SampleType",
##'                           samplePattern = "Macrophage",
##'                           refPattern = "Ref")
##'                           
divideByReference <- function(object, 
                              i, 
                              colvar, 
                              samplePattern = ".", 
                              refPattern) {
    ## Check arguments
    if (!inherits(object, "QFeatures")) 
        stop("'object' must be a QFeatures object")
    for (ii in i){
        ## Get the reference channel 
        annot <- colData(object)[colnames(object[[ii]]), ][, colvar]
        refIdx <- grep(refPattern, annot)
        sampIdx <- grep(samplePattern, annot)
        if (!length(refIdx)) 
            stop("The reference pattern '", refPattern, 
                 "' did not match any column in '", names(object)[ii], "'")
        if (!length(sampIdx)) 
            stop("The sample pattern '", samplePattern, 
                 "' did not match any column in '", names(object)[ii], "'")
        if (length(refIdx) != 1) 
            warning("Multiple references found in assay '", names(object)[ii], 
                    "'. Only the first match will be used")
        ## Divide all channels by the reference channel
        y <- object[[ii]]
        ref <- assay(y)[, refIdx]
        assay(y)[, sampIdx] <- assay(y)[, sampIdx, drop = FALSE] / ref
        ## Store the normalized assay
        object <- replaceAssay(object, y, ii)
    }
    return(object)
}

##' Normalize single-cell proteomics (SCP) data
##'
##' This function normalises an assay in a `QFeatures` according to 
##' the supplied method (see Details). The normalized data is added as
##' a new assay 
##' 
##' @param object An object of class `QFeatures`.
##' 
##' @param i A numeric vector or a character vector giving the index 
##'     or the name, respectively, of the assay(s) to be processed.
##'
##' @param name A `character(1)` naming the new assay name. Defaults is
##'     are `normAssay`.
##' 
##' @param method `character(1)` defining the normalisation method to
##'     apply. See Details.`
##' 
##' @param ... Additional parameters passed to 
##'     [MsCoreUtils::normalizeMethods()].
##' 
##' @return A `QFeatures` object with an additional assay containing the 
##'     normalized data.
##' 
##' @details 
##' 
##' The `method` parameter in `normalize` can be one of `"sum"`,
##' `"max"`, `"center.mean"`, `"center.median"`, `"div.mean"`,
##' `"div.median"`, `"diff.meda"`, `"quantiles`", `"quantiles.robust`"
##' or `"vsn"`. The [MsCoreUtils::normalizeMethods()] function returns
##' a vector of available normalisation methods.
##'
##' - For `"sum"` and `"max"`, each feature's intensity is divided by
##'   the maximum or the sum of the feature respectively. These two
##'   methods are applied along the features (rows).
##'
##' - `"center.mean"` and `"center.median"` center the respective
##'   sample (column) intensities by subtracting the respective column
##'   means or medians. `"div.mean"` and `"div.median"` divide by the
##'   column means or medians. These are equivalent to `sweep`ing the
##'   column means (medians) along `MARGIN = 2` with `FUN = "-"` (for
##'   `"center.*"`) or `FUN = "/"` (for `"div.*"`).
##'
##' - `"diff.median"` centers all samples (columns) so that they all
##'   match the grand median by subtracting the respective columns
##'   medians differences to the grand median.
##'
##' - Using `"quantiles"` or `"quantiles.robust"` applies (robust) 
##'   quantile normalisation, as implemented in 
##'   [preprocessCore::normalize.quantiles()] and 
##'   [preprocessCore::normalize.quantiles.robust()]. `"vsn"` uses the
##'   [vsn::vsn2()] function. Note that the latter also glog-transforms
##'   the intensities. See respective manuals for more details and 
##'   function arguments.
##'
##' For further details and examples about normalisation, see
##' [MsCoreUtils::normalize_matrix()].
##'
##' @seealso [QFeatures::normalize] for more details about `normalize`
##'
##' @export
##' 
##' @examples 
##' 
##' data("scp1")
##' scp1
##' normalizeSCP(scp1, i = "proteins", name = "normproteins", 
##'              method = "center.mean")
##' 
normalizeSCP <- normaliseSCP <- function(object, i, name = "normAssay", 
                                         method, ...) {
    if(!inherits(object, "QFeatures"))
      stop("'object' must be a 'QFeatures' object.")
    if (length(i) != 1)
      stop("Only one assay to be processed at a time")
    if (is.numeric(i)) i <- names(object)[[i]]
    object <- addAssay(object,
                       .normalizeSCP(object[[i]], method, ...),
                       name)
    addAssayLinkOneToOne(object, from = i, to = name)
}
