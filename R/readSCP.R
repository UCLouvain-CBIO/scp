##' @title Read single-cell proteomics tabular data
##'
##' @description
##'
##' Function to import and convert tabular data from a spreadsheet or
##' a `data.frame` into a `SingleCellExperiment` and `QFeatures`
##' object.
##'
##' @param ... Parameters passed to [readSummarizedExperiment()],
##'     [readQFeatures()] or [readQFeaturesFromDIANN()]. See these
##'     respective manual pages for details.
##'
##' @param experimentsAsSce A `logical(1)` indicating whether the
##'     `QFeatures` object should be composed of `SingleCellExperiment`
##'     objects. By default: `FALSE`, the `QFeatures` object will be composed
##'     of `SummarizedExperiment` objects. Note that using
##'     `SingleCellExperiment` can impact the performance.
##'
##' @return An instance of class `SingleCellExperiment` or a
##'     `QFeatures`, composed of `SingleCellExperiment` or
##'     `SummarizedExperiment` objects.
##'
##' @note The `SingleCellExperiment` class is built on top of the
##'     `RangedSummarizedExperiment` class. This means that some column names
##'     are forbidden in the `rowData`. Avoid using the following names:
##'     `seqnames`, `ranges`, `strand`, `start`, `end`,
##'     `width`,  `element`
##'
##' @importFrom methods as
##' @import SingleCellExperiment
##' @import SummarizedExperiment
##' @importFrom MultiAssayExperiment ExperimentList experiments<-
##' @import QFeatures
##'
##' @seealso
##'
##' - The more general [QFeatures::readQFeatures()] function, which
##'   this function depends on.
##'
##' - The more general [QFeatures::readQFeaturesFromDIANN()] function,
##'   for details and an example on how to read label-free and plexDIA
##'   (mTRAQ) data processed with DIA-NN.
##'
##' - The [QFeatures::readSummarizedExperiment()] function, which
##' `readSingleCellExperiment()` depends on.
##'
##' - The [SingleCellExperiment::SingleCellExperiment()] class.
##'
##' @export
##'
##' @examples
##'
##' ######################################################
##' ## Load a single acquisition as a SingleCellExperiment
##'
##' ## Load a data.frame with PSM-level data
##' data("mqScpData")
##'
##' ## Create the QFeatures object
##' quantCols <- grep("Reporter.intensity.\\d", colnames(mqScpData))
##' sce <- readSingleCellExperiment(mqScpData, quantCols)
##' sce
##'
##' ######################################################
##' ## Load multiple acquisitions as a QFeatures
##'
##' ## Load an example table containing MaxQuant output
##' data("mqScpData")
##'
##' ## Load the (user-generated) annotation table
##' data("sampleAnnotation")
##'
##' ## Format the tables into a QFeatures object
##' readSCP(assayData = mqScpData,
##'         colData = sampleAnnotation,
##'         runCol = "Raw.file")
readSCP <- function(..., experimentsAsSce = FALSE) {
    ans <- readQFeatures(...)
    if (experimentsAsSce) {
        el <- ExperimentList(lapply(experiments(ans),
                                as, "SingleCellExperiment"))
        experiments(ans) <- el
    }
    setQFeaturesType(ans, "scp")
}


##' @export
##'
##' @rdname readSCP
readSCPfromDIANN <- function(..., experimentsAsSce = FALSE) {
    ans <- readQFeaturesFromDIANN(...)
    if (experimentsAsSce) {
        el <- ExperimentList(lapply(experiments(ans),
                                as, "SingleCellExperiment"))
        experiments(ans) <- el
    }
    setQFeaturesType(ans, "scp")
}

##' @export
##'
##' @rdname readSCP
readSingleCellExperiment <- function(...) {
    ## Read data as SummarizedExperiment
    sce <- readSummarizedExperiment(...)
    as(sce, "SingleCellExperiment")
}
