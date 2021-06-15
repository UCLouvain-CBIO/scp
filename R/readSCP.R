##' @import SingleCellExperiment
##' @import QFeatures
##' @import dplyr
##' @import magrittr

##' @title Read single-cell proteomics data as a QFeatures object from
##'     tabular data and metadata
##'
##' @description
##'
##' Convert tabular quantitative MS data and metadata from a
##' spreadsheet or a `data.frame` into a [QFeatures] object containing
##' [SingleCellExperiment] objects.
##'
##' @param featureData File or object holding the quantitative
##'     data. Can be either a `character(1)` with the path to a
##'     text-based spreadsheet (comma-separated values by default, but
##'     see `...`) or an object that can be coerced to a
##'     `data.frame`. It is advised not to encode characters as
##'     factors.
##' 
##' @param colData A `data.frame` or any object that can be coerced
##'     to a `data.frame`. `colData` is expected to contains all the
##'     sample meta information. Required fields are the acquisition
##'     batch (given by `batchCol`) and the acquisition channel within
##'     the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     meta data.
##'
##' @param batchCol A `numeric(1)` or `character(1)` pointing to the
##'     column of `featureData` and `colData` that contain the batch
##'     names. Make sure that the column name in both table are either
##'     identical and syntactically valid (if you supply a `character`)
##'     or have the same index (if you supply a `numeric`). Note that
##'     characters can be converted to syntactically valid names using
##'     `make.names`
##' 
##' @param channelCol A `numeric(1)` or `character(1)` pointing to the
##'     column of `colData` that contains the column names of the
##'     quantitative data in `featureData` (see Example).
##' 
##' @param suffix A `character()` giving the suffix of the column 
##'     names in each assay. The length of the vector should equal the
##'     number of quantification channels and should contain unique 
##'     character elements. If NULL, the names of the quantification 
##'     columns in `featureData` are taken as suffix. 
##' 
##' @param removeEmptyCols A `logical(1)`. If true, the function will
##'     remove in each batch the columns that contain only missing 
##'     values.
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
##' @note The `SingleCellExperiment` class is built on top of the 
##'     `RangedSummarizedExperiment` class. This means that some column names 
##'     are forbidden in the `rowData`. Avoid using the following names:
##'     `seqnames`, `ranges`, `strand`, `start`, `end`, 
##'     `width`,  `element`
##'     
##' @author Laurent Gatto, Christophe Vanderaa
##' 
##' @importFrom utils read.csv
##' @importFrom S4Vectors DataFrame
##' @importFrom MultiAssayExperiment ExperimentList
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
##' readSCP(featureData = mqScpData,
##'         colData = sampleAnnotation,
##'         batchCol = "Raw.file",
##'         channelCol = "Channel")
##' 
readSCP <- function(featureData, 
                    colData, 
                    batchCol, 
                    channelCol,
                    suffix = NULL,
                    removeEmptyCols = FALSE,
                    verbose = TRUE,
                    ...) {
    ## Check the batch column name
    if (!identical(make.names(batchCol), batchCol)) 
        stop("'batchCol' is not a syntactically valid column name. ",
             "See '?make.names' for converting the column names to ",
             "valid names, e.g. '", batchCol, "' -> '", 
             make.names(batchCol), "'")
    
    colData <- as.data.frame(colData)
    
    ## Get the column contain the expression data
    ecol <- unique(colData[, channelCol])
    ## Get the sample suffix
    if (is.null(suffix))
        suffix <- ecol
    
    ## Create the SingleCellExperiment object
    if (verbose) message("Loading data as a 'SingleCellExperiment' object")
    scp <- readSingleCellExperiment(table = featureData, 
                                    ecol = ecol, 
                                    ...)
    if (is.null(list(...)$row.names))
        rownames(scp) <- paste0("PSM", seq_len(nrow(scp)))
    
    ## Check the link between colData and scp
    mis <- !rowData(scp)[, batchCol] %in% colData[, batchCol]
    if (any(mis)) {
        warning("Missing metadata. The features are removed for ", 
                paste0(unique(rowData(scp)[mis, batchCol]), collapse = ", "))
        scp <- scp[!mis, ]
    }
    
    ## Split the SingleCellExperiment object by batch column
    if (verbose) message(paste0("Splitting data based on '", batchCol, "'"))
    scp <- .splitSCE(scp, f = batchCol)
    
    ## Clean each element in the data list
    for (i in seq_along(scp)) {
        ## Add unique sample identifiers
        colnames(scp[[i]]) <- paste0(names(scp)[[i]], suffix)
        ## Remove the columns that are all NA
        if (removeEmptyCols) {
            sel <- colSums(is.na(assay(scp[[i]]))) != nrow(scp[[i]])
            scp[[i]] <- scp[[i]][, sel]
        } 
    }
    
    if (verbose) message(paste0("Formatting sample metadata (colData)"))
    ## Create the colData 
    cd <- DataFrame(row.names = unlist(lapply(scp, colnames)))
    rownames(colData) <- paste0(colData[, batchCol], suffix)
    cd <- cbind(cd, colData[rownames(cd), ])
    
    ## Store the data as a QFeatures object and add the experimental
    ## information
    if (verbose) message("Formatting data as a 'QFeatures' object")
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
##' @note The `SingleCellExperiment` class is built on top of the 
##'     `RangedSummarizedExperiment` class. This means that some column names 
##'     are forbidden in the `rowData`. Avoid using the following names:
##'     `seqnames`, `ranges`, `strand`, `start`, `end`, 
##'     `width`,  `element`
##'     
##' 
##' @seealso The code relies on 
##'     [QFeatures::readSummarizedExperiment].
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
##' data("mqScpData")
##' 
##' ## Create the QFeatures object
##' sce <- readSingleCellExperiment(mqScpData, 
##'                                 grep("RI", colnames(mqScpData)))
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
.splitSCE <- function(x, 
                      f) {
    ## Check that f is a factor
    if (is.character(f)) {
        if (length(f) != 1) 
            stop("'f' must be of lenght one")
        if (f %in% colnames(rowData(x))) {
            f <- rowData(x)[, f]
        }
        else if (f %in% colnames(colData(x))) {
            f <- colData(x)[, f]
        }
        else {
            stop("'", f, "' not found in rowData or colData")
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
