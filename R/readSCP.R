##' @import SingleCellExperiment
##' @import QFeatures
##' @import tidyverse


##' @title Read single-cell proteomics data as a QFeatures object from
##'     tabular data and metadata
##'
##' @description
##'
##' Convert tabular quantitative MS data and metadata from a
##' spreadsheet or a `data.frame` into a [QFeatures] object containing
##' [SingleCellExperiment] objects.
##'
##' @param quantTable File or object holding the quantitative
##'     data. Can be either a `character(1)` with the path to a
##'     text-based spreadsheet (comma-separated values by default, but
##'     see `...`) or an object that can be coerced to a
##'     `data.frame`. It is advised not to encode characters as
##'     factors.
##' 
##' @param metaTable A `data.frame` or any object that can be coerced
##'     to a `data.frame`. `metaTable` is expected to contains all the
##'     sample meta information. Required fields are the acquisition
##'     batch (given by `batchCol`) and the acquisition channel within
##'     the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     meta data.
##'
##' @param batchCol A `numeric(1)` or `character(1)` pointing to the
##'     column of `quantTable` and `metaTable` that contain the batch
##'     names. Make sure that the column name in both table are either
##'     identical (if you supply a `character`) or have the same index
##'     (if you supply a `numeric`).
##' 
##' @param channelCol A `numeric(1)` or `character(1)` pointing to the
##'     column of `metaTable` that contains the column names of the
##'     quantitive data in `quantTable` (see Example).
##' 
##' @param verbose A `logical(1)` indicating whether the progress of
##'     the data reading and formatting should be printed to the
##'     console. Default is `TRUE`.
##' 
##' @param ... Further arguments that can be passed on to [read.csv]
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @return An instance of class [QFeatures]. The expression data of
##'     each batch is stored in a separate assay as a
##'     [SingleCellExperiment] object.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##' 
##' @importFrom utils read.csv
##'
##' @md
##' @export
##' 
##' @examples 
##' 
##' ## Load an example table containing MaxQuant output
##' data("mqFile")
##' 
##' ## Load the (user-generated) annotation table
##' data("sampleAnnotation")
##' 
##' ## Format the tables into a QFeatures object
##' readSCP(quantTable = mqFile,
##'         metaTable = sampleAnnotation,
##'         batchCol = "Set",
##'         channelCol = "Channel")
##' 
readSCP <- function(quantTable, 
                    metaTable,
                    batchCol,
                    channelCol, 
                    verbose = TRUE,
                    ...) {
    metaTable <- as.data.frame(metaTable)
    ## Create the SingleCellExperiment object
    if (verbose) message("Loading data as a 'SingleCellExperiment' object")
    ecol <- unique(metaTable[, channelCol])
    scp <- readSingleCellExperiment(table = quantTable, 
                                    ecol = ecol, 
                                    ...)
    if (is.null(list(...)$row.names))
        rownames(scp@assays@data@listData[[1]]) <- paste0("PSM", 1:nrow(scp))
        ## Note the row names 
    
    ## Check the link between metaTable and scp
    mis <- !rowData(scp)[, batchCol] %in% metaTable[, batchCol]
    if (any(mis)) {
        warning("Missing metadata. The features are removed for ", 
                paste0(unique(rowData(scp)[mis, batchCol]), collapse = ", "))
      scp <- scp[!mis, ]
    }
    
    ## Split the SingleCellExperiment object by batch column
    if (verbose) message(paste0("Splitting data based on '", batchCol))
    scp <- .splitSCE(scp, f = batchCol)
    
    ## Add unique sample identifiers
    if (verbose) message(paste0("Formating sample metadata (colData)"))
    for (i in 1:length(scp)) {
        colnames(scp[[i]]) <- paste0(names(scp)[[i]], "_", colnames(scp[[i]]))
    }
    ## Create the colData 
    cd <- DataFrame(row.names = unlist(lapply(scp, colnames)))
    rownames(metaTable) <- paste0(metaTable[, batchCol], "_", 
                                  metaTable[, channelCol])
    cd <- cbind(cd, metaTable[rownames(cd), ])
    
    ## Store the data as a QFeatures object and add the experimental
    ## information
    if (verbose) message("Formating data as a 'QFeatures' object")
    QFeatures(experiments = scp, 
              colData = cd)
}

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
##'     dislpay the column names as seen by the `read.csv` function.
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
##' @seealso The code relies on [QFeatures::readSummarizedExperiment].
##'
##'
##' @md
##' 
##' @export
##'
##' @importFrom methods as
##'
##' @examples 
##' ## Load a data.frame with PSM-level data
##' data(mqFile)
##' 
##' ## Create the QFeatures object
##' sce <- readSingleCellExperiment(mqFile, grep("RI", colnames(mqFile)))
readSingleCellExperiment <- function(table, 
                                     ecol, 
                                     fnames,
                                     ...) {
    ## Read data as SummarizedExperiment
    sce <- readSummarizedExperiment(table, ecol, fnames, ...)
    sce <- as(sce, "SingleCellExperiment")
    return(sce)
}


##' Split SingleCellExperiment into an ExperimentList
##'
##' The fonction creates an [ExperimentList] containing
##' [SingleCellExperiment] objects from a [SingleCellExperiment]
##' object. `f` is used to split `x`` along the rows (`f`` was a feature
##' variable name) or samples/columns (f was a phenotypic variable
##' name). If f is passed as a factor, its length will be matched to
##' nrow(x) or ncol(x) (in that order) to determine if x will be split
##' along the features (rows) or sample (columns). Hence, the length of
##' f must match exactly to either dimension.
##' 
##' This function is not exported. If this is needed, create a pull
##' request to `rformassspectrometry/QFeatures`.
##' 
##' @param x a single [SingleCellExperiment] object
##' 
##' @param f a factor or a character of length 1. In the latter case,
##'     `f` will be matched to the row and column data variable names
##'     (in that order). If a match is found, the respective variable
##'     is extracted, converted to a factor if needed
##' @noRd
.splitSCE <- function(x, f) {
    ## Check that f is a factor
    if (is.character(f)) {
      if (length(f) != 1) 
          stop("Character must be of lenght one.")
      if (f %in% colnames(rowData(x))) {
          f <- rowData(x)[, f]
      }
      else if (f %in% colnames(colData(x))) {
          f <- colData(x)[, f]
      }
      else {
          stop(f, " not found in any feature/phenodata variables.")
      }
      if (!is.factor(f)) 
          f <- factor(f)
    }
    ## Check that the factor matches one of the dimensions
    if (!length(f) %in% dim(x)) 
        stop("length(f) not compatible with dim(x).")
    if (length(f) == nrow(x)) { ## Split along rows
        xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
    } else { ## Split along columns
        xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
    }
    ## Convert list to an ExperimentList
    do.call(ExperimentList, xl)
}
