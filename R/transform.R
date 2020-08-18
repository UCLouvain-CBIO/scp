##' Internal function to efficiently replace a `SummarizedExperiment`
##' (and subclasses of it) assay in a `QFeatures` object. Watch out
##' this function should only be used when the old assay and the assay
##' to replace have same size. Otherwise, this can lead to
##' inconsistencies in the colData and functionality from
##' `MultiAssayExperiment` should be used.
##' 
##' @param obj A QFeatures object
##' 
##' @param y A SummarizedExperiment object, or an object that inherits from it
##' 
##' @param i a characther(1) or logical(1) that indicates which assay must be replaced
##'
##' @noRd
.replaceAssay <- function(obj, y, i) {
  if (length(i) > 1) stop("Only 1 assay can be replaced at a time.")
  if (!inherits(obj, "QFeatures")) stop("'obj' must be a 'QFeatures' object")
  if (!inherits(y, "SummarizedExperiment")) 
    stop("'y' must inherits from a 'SummarizedExperiment' object")
  if (!identical(colnames(experiments(obj)[[i]]), colnames(y)))
    stop("Colnames of old and new assays must match. Otherwise, consider ",
         "using 'experiments(obj)[[i]] <- y' to avoid bad surprises.")
  obj@ExperimentList@listData[[i]] <- y
  return(obj)
}

##' Divide assay columns by a reference column
##'
##' The function divides the sample columns by a reference column. The sample 
##' and reference columns are defined based on the provided `colDataCol` 
##' variable and on regular expression matching.
##' 
##' The supplied assay(s) are replaced with the values computed after reference
##' division.
##' 
##' @param obj A `QFeatures` object
##' 
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' 
##' @param colDataCol A `character(1)` indicating the variable to take from 
##'     `colData(obj)` that gives the sample annotation.
##' 
##' @param samplePattern A `character(1)` pattern that matches the sample 
##'     encoding in `colDataCol`. By default all samples are devided (using the
##'     regex wildcard `.`).
##' 
##' @param refPattern A `character(1)` pattern that matches the carrier 
##'     encoding in `colDataCol`. Only one match per assay is allowed, otherwise
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
##'                           colDataCol = "SampleType",
##'                           samplePattern = "Macrophage",
##'                           refPattern = "Ref")
divideByReference <- function(obj, 
                              i,
                              colDataCol,
                              samplePattern = ".",
                              refPattern) {
  ## Check arguments
  if (!inherits(obj, "QFeatures")) stop("'obj' must be a QFeatures object")
  for (ii in i){
    ## Get the reference channel 
    annot <- colData(obj)[colnames(obj[[ii]]), ][, colDataCol]
    refIdx <- grep(refPattern, annot)
    sampIdx <- grep(samplePattern, annot)
    if (!length(refIdx)) 
      stop("The reference pattern '", refPattern, 
           "' did not match any column in '", names(obj)[ii], "'")
    if (!length(sampIdx)) 
      stop("The sample pattern '", samplePattern, 
           "' did not match any column in '", names(obj)[ii], "'")
    if (length(refIdx) != 1) 
      warning("Multiple references found in assay '", names(obj)[ii], 
              "'. Only the first match will be used")
    ## Divide all channels by the reference channel
    y <- obj[[ii]]
    assay(y)[, sampIdx] <- assay(y)[, sampIdx, drop = FALSE] / assay(y)[, refIdx]
    ## Store the normalized assay
    obj <- .replaceAssay(obj, y, ii)
  }
  return(obj)
}


##' Remove infinite data 
##' 
##' This function coinsiders any infinite value as missing data. So,
##' any value in assay `i` that return `TRUE` to `is.infinite` will be
##' replaced by `NA`.
##' 
##' @param obj A `QFeatures` object.
##' 
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' 
##' @return A `QFeatures` object
##' 
##' @export
##' 
##' @examples
##' data("scp1")
##' scp1 <- infIsNA(scp1, i = "peptides")
infIsNA <- function(obj, i) {
  for (ii in i) {
    sel <- is.infinite(assay(obj[[ii]])) 
    assay(obj[[ii]])[sel] <- NA
  }
  obj
}
