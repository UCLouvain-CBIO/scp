##' @import SingleCellExperiment
##' @import QFeatures

##' @title Read SingleCellExperiment from tabular data
##'
##' @description
##'
##' Convert tabular data from a spreadsheet or a `data.frame` into a
##' `SingleCellExperiment` object.
##'
##' @param table File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`. It is
##'     advised not to encode characters as factors.
##'
##' @param ecol A `numeric` indicating the indices of the columns to
##'     be used as assay values. Can also be a `character`
##'     indicating the names of the columns. Caution must be taken if
##'     the column names are composed of special characters like `(`
##'     or `-` that will be converted to a `.` by the `read.csv`
##'     function. If `ecol` does not match, the error message will
##'     display the column names as seen by the `read.csv` function.
##'
##' @param fnames An optional `character(1)` or `numeric(1)`
##'     indicating the column to be used as row names.
##'
##' @param ... Further arguments that can be passed on to [read.csv]
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @return An instance of class [SingleCellExperiment].
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @note The `SingleCellExperiment` class is built on top of the
##'     `RangedSummarizedExperiment` class. This means that some
##'     column names are forbidden in the `rowData`. Avoid using the
##'     following names: `seqnames`, `ranges`, `strand`, `start`,
##'     `end`, `width`, `element`
##'
##' @seealso The code relies on
##'     [QFeatures::readSummarizedExperiment()].
##'
##' @md
##'
##' @export
##'
##' @importFrom methods as
##'
##' @examples
##' ## Load a data.frame with PSM-level data
##' data("mqScpData")
##'
##' ## Create the QFeatures object
##' sce <- readSingleCellExperiment(mqScpData,
##'                                 grep("RI", colnames(mqScpData)))
readSingleCellExperiment <- function(table, ecol,
                                     fnames, ...) {
    ## Read data as SummarizedExperiment
    sce <- readSummarizedExperiment(table, ecol, fnames, ...)
    as(sce, "SingleCellExperiment")
}

##' @title Single-cell proteomics data from tabular data
##'
##' @description
##'
##' Converts tabular quantitative MS data and metadata from a
##' `data.frame` into a `QFeatures` object containing
##' `SingleCellExperiment` sets.
##'
##' @inheritParams QFeatures::readQFeatures
##'
##' @param colAnnotation A `data.frame` or any object that can be
##'     coerced to a `data.frame`. It is expected to contain all the
##'     sample meta information. Required fields are the acquisition
##'     batch (given by `runCol`) and the acquisition channel within
##'     the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     meta data.
##'
##' @return An instance of class [QFeatures::QFeatures()]. The
##'     expression data of each batch is stored in a separate assay as
##'     a [SingleCellExperiment::SingleCellExperiment()] object.
##'
##' @note The `SingleCellExperiment` class is built on top of the
##'     `RangedSummarizedExperiment` class. This means that some
##'     column names are forbidden in the `rowData`. Avoid using the
##'     following names: `seqnames`, `ranges`, `strand`, `start`,
##'     `end`, `width`, `element`
##'
##' @seealso
##'
##' - The more general [QFeatures::readQFeatures()] function, which
##'   this function depends on.
##'
##' - The [SingleCellExperiment::SingleCellExperiment()] class.
##'
##' - The [readSCPfromDIANN()] function to import DIA-NN quantitative
##'   single-cell data.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @importFrom S4Vectors DataFrame
##' @importFrom MultiAssayExperiment ExperimentList experiments experiments<-
##' @importFrom SummarizedExperiment colData rowData assay
##' @importFrom SummarizedExperiment rowData<- colData<- assay<-
##'
##' @md
##' @export
##'
##' @examples
##'
##' ## Load an example table containing MaxQuant output
##' data("mqScpData")
##'
##' ## Load the (user-generated) annotation table
##' data("sampleAnnotation")
##'
##' ## Format the tables into a QFeatures object
##' readSCP(assayData = mqScpData,
##'         colAnnotation = sampleAnnotation,
##'         runCol = "Raw.file",
##'         channelCol = "Channel")
readSCP <- function(assayData, colAnnotation, runCol, channelCol,
                    suffix = NULL, sep = "", removeEmptyCols = FALSE,
                    verbose = TRUE) {
    ans <- readQFeatures(assayData = assayData,
                         colAnnotation = colAnnotation,
                         runCol = runCol,
                         channelCol = channelCol,
                         suffix = suffix,
                         sep = sep,
                         removeEmptyCols = removeEmptyCols,
                         verbose = verbose)
    el <- ExperimentList(lapply(experiments(ans),
                                as, "SingleCellExperiment"))
    experiments(ans) <- el
    ans
}


##' @title  Read single-cell DIA-NN output as a QFeatures objects
##'
##' @description
##'
##' This function takes the output tables from DIA-NN and converts
##' them into a QFeatures object using the scp framework.
##'
##' @inheritParams QFeatures::readQFeaturesFromDIANN
##'
##' @return An instance of class `QFeatures`. The expression data of
##'     each acquisition run is stored in a separate set as a
##'     `SingleCellExperiment` object.
##'
##' @export
readSCPfromDIANN <- function(colAnnotation, assayData,
                             extractedData = NULL, ecol = "MS1.Area",
                             multiplexing = "none", # "none" or "mTRAQ"
                             ...) {
    ans <- readQFeaturesFromDIANN(colAnnotation, assayData,
                                  extractedData,
                                  ecol, ...)
    el <- ExperimentList(lapply(experiments(ans),
                                as, "SingleCellExperiment"))
    experiments(ans) <- el
    ans
}
