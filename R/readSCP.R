##' @import SingleCellExperiment
##' @import QFeatures
##' @import dplyr

##' @title Read single-cell proteomics data as a QFeatures object from
##'     tabular data and metadata
##'
##' @description
##'
##' Convert tabular quantitative MS data and metadata from a
##' spreadsheet or a `data.frame` into a [QFeatures] object containing
##' [SingleCellExperiment] objects.
##'
##' @param featureData File or object holding the identification and
##'     quantitative data. Can be either a `character(1)` with the
##'     path to a text-based spreadsheet (comma-separated values by
##'     default, but see `...`) or an object that can be coerced to a
##'     `data.frame`. It is advised not to encode characters as
##'     factors.
##'
##' @param colData A `data.frame` or any object that can be coerced
##'     to a `data.frame`. `colData` is expected to contain all the
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
##'     names in each assay. Sample/single-cell (column) names are
##'     automatically generated using: batch name + sep + suffix. Make
##'     sure suffix contains unique character elements. The length of
##'     the vector should equal the number of quantification channels.
##'     If NULL (default), the suffix is derived from the the names of
##'     the quantification columns in `featureData`.
##'
##' @param sep A `character(1)` that is inserted between the assay
##'     name and the `suffix` (see `suffix` argument for more details).
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
readSCP <- function(featureData, colData, batchCol, channelCol,
                    suffix = NULL, sep = "", removeEmptyCols = FALSE,
                    verbose = TRUE, ...) {
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
                                    ecol = ecol, ...)
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
    if (verbose) message("Splitting data based on '", batchCol, "'")
    scp <- .splitSCE(scp, f = batchCol)
    ## Clean each element in the data list
    for (i in seq_along(scp)) {
        ## Add unique sample identifiers
        colnames(scp[[i]]) <- paste0(names(scp)[[i]], sep, suffix)
        ## Remove the columns that are all NA
        if (removeEmptyCols) {
            sel <- colSums(is.na(assay(scp[[i]]))) != nrow(scp[[i]])
            scp[[i]] <- scp[[i]][, sel]
        }
    }
    if (verbose) message("Formatting sample metadata (colData)")
    ## Create the colData
    cd <- DataFrame(row.names = unlist(lapply(scp, colnames)))
    rownames(colData) <- paste0(colData[, batchCol], sep, suffix)
    cd <- cbind(cd, colData[rownames(cd), ])
    ## Store the data as a QFeatures object
    if (verbose) message("Formatting data as a 'QFeatures' object")
    QFeatures(experiments = scp, colData = cd)
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


##' @title  Read DIA-NN output as a QFeatures objects for single-cell 
##' proteomics data
##'
##' @description
##' 
##' This function takes the output tables from DIA-NN and converts them
##' into a QFeatures object using the scp framework. 
##'
##' @param colData A data.frame or any object that can be coerced to a
##'     data.frame. colData is expected to contains all the sample 
##'     annotations. We require the table to contain a column called
##'     `File.Name` that links to the `File.Name` in the DIA-NN report
##'     table. If `multiplexing = "mTRAQ"`, we require a second column
##'     called `Label` that links the label to the sample (the labels
##'     identified by DIA-NN can be retrieved from `Modified Sequence`
##'     column in the report table).
##'
##' @param reportData A data.frame or any object that can be coerced 
##'     to a data.frame that contains the data from the `Report.tsv`
##'     file generated by DIA-NN. 
##'
##' @param extractedData A data.frame or any object that can be coerced 
##'     to a data.frame that contains the data from the `*_ms1_extracted.tsv`
##'     file generated by DIA-NN. This argument is optional and is 
##'     currently only applicable for mTRAQ multiplexed experiments 
##'     where DIA-NN was run using the `plexdia` module 
##'
##' @param ecol A `character(1)` indicating which column in 
##'     `reportData` contains the quantitative information.
##'
##' @param multiplexing A `character(1)` indicating the type of 
##'     multiplexing used in the experiment. Provide `"none"` if the
##'     experiment is label-free (default). Alternative options are:
##'     `"mTRAQ"`. 
##'
##' @param ... Further arguments passed to `readSCP()`
##'
##' @return An instance of class QFeatures. The expression data of
##'     each acquisition run is stored in a separate assay as a
##'     SingleCellExperiment object.
##' 
##' @importFrom tidyr pivot_wider
##' 
##' @export
readSCPfromDIANN <- function(colData, reportData, extractedData = NULL,
                             ecol = "Ms1.Area", multiplexing = "none",
                             ...) {
    suppArgs <- .checkDiannArguments(
        colData, reportData, extractedData, ecol, multiplexing, ...
    )
    if (multiplexing == "mTRAQ") {
        reportData <- .formatMtraqReportData(reportData, colData, ecol)
    } else if (multiplexing == "none") {
        colData$Label <- ecol
    } 
    allArgs <- c(suppArgs, list(
        featureData = reportData, colData = colData,
        batchCol = "File.Name", channelCol = "Label"
    ))
    out <- do.call(readSCP, allArgs)
    if (!is.null(extractedData)) {
        out <- .addExtractedData(out, extractedData)
    }
    out
}

## Internal function that checks whether the provided arguments match
## the expected input that is generated by DIA-NN. 
## Parameter description is the same as for `readSCPfromDIANN()` 
.checkDiannArguments <- function(colData, reportData, extractedData, 
                                 ecol, multiplexing, ...) {
    diannReportCols <- c("File.Name", "Precursor.Id", "Modified.Sequence")
    if (!all(diannReportCols %in% colnames(reportData)))
        stop("'reportData' is not an expected DIA-NN report table ",
             "output. This function expects the main output file as ",
             "described here: https://github.com/vdemichev/DiaNN#main-output-reference")
    if (!ecol %in% colnames(reportData))
        stop("'", ecol, "' not found in 'reportData'")
    if (!"File.Name" %in% colnames(colData)) 
        stop("'colData' must contain a column named 'File.Name' that provides ",
             "a link to the 'File.Name' column in 'reportData'")
    if (multiplexing == "none" && !is.null(extractedData)) 
        stop("Providing 'extractedData' for label-free experiments ",
             "('multiplexed == \"none\"') is not expected. Raise an ",
             "issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")
    .checkDiannArgumentsDots(multiplexing, ...)
}

## Internal function that adapts the dots arguments (that will be used
## by `readSCP()`) depending on the multiplexing approach used.
.checkDiannArgumentsDots <- function(multiplexing, ...) {
    suppArgs <- list(...)
    if (multiplexing == "mTRAQ") {
        suppArgs$sep <- "."
    } else if (multiplexing == "none") {
        suppArgs$sep <- ""
        suppArgs$suffix <- ""
    } else {
        stop("The '", multiplexing, "' multiplexing strategy is not ",
             "implemented. Raise an issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")
    }
    suppArgs
}

## (Only for mTRAQ multiplexing!) Internal function that extracts the
## mTRAQlabels from the peptide sequence, removes the mTRAQ annotation
## from the precursor ID, identifies constant columns within precursor
## and puts the quantification data for different mTRAQ labels in 
## separate columns (wide format).
.formatMtraqReportData <- function(reportData, colData, ecol) {
    reportData$Label <-
        sub("^.*[Q-](\\d).*$", "\\1", reportData$Modified.Sequence)
    reportData$Precursor.Id <-
        gsub("\\(mTRAQ.*?\\)", "(mTRAQ)", reportData$Precursor.Id)
    .checkLabelsInColData(colData, reportData)
    idCols <- .findPrecursorVariables(reportData)
    pivot_wider(
        reportData, id_cols = all_of(idCols), 
        names_from = "Label", values_from = all_of(ecol)
    )
}

## Internal function that identifies which variables in the report 
## data are constant within each precursor (with each run).
.findPrecursorVariables <- function(reportData) {
    precIds <- paste0(reportData$Precursor.Id, reportData$File.Name)
    nUniqueIds <- length(unique(precIds))
    nLevels <- sapply(colnames(reportData), function(x) {
        nrow(unique(reportData[, c("Precursor.Id", "File.Name", x)]))
    })
    names(nLevels)[nLevels == nUniqueIds]
}

## Internal function that ensures that the reportData and the colData
## are correctly linked.
.checkLabelsInColData <- function(colData, reportData) {
    if (!"Label" %in% colnames(colData))
        stop("'colData' must contain a column named 'Label' that ",
             "provides the mTRAQ reagent used to label the ",
             "samples and/or single cells.")
    if (any(mis <- !colData$Label %in% reportData$Label)) {
        stop("Some labels from 'colData$Label' were not found as",
             "part of the mTRAQ labels found in ",
             "'reportData$Modified.Sequence': ",
             paste0(unique(colData$Label[mis]), collapse = ", "))
    }
    NULL
}

## Internal function that adds the extractedData to a QFeatures 
## object. The functions first converts the extractedData to a 
## SingleCellExperiment objects and subsets the SCE for the set of 
## shared samples. The added assay is automatically linked (using 
## AssayLinks) to the reportData.
## Developer's note: the function assumes that DIA-NN creates sample
## names in the extracted data by appending the labels to the run 
## names
## @param object A QFeatures object containing DIA-NN report data, as
##     generated by readSCP 
## @param extractedData A data.frame or any object that can be coerced 
##     to a data.frame that contains the data from the `*_ms1_extracted.tsv`
##     file generated by DIA-NN.
.addExtractedData <- function(object, extractedData) {
    quantColPattern <- paste0(unique(object$Label), "$", collapse = "|")
    quantCols <- grep(quantColPattern, colnames(extractedData))
    extractedData <- readSingleCellExperiment(
        extractedData, ecol = quantCols, fnames = "Precursor.Id"
    )
    extractedData <- .keepSharedSamples(extractedData, object)
    object <- addAssay(object, extractedData, name = "Ms1Extracted")
    addAssayLink(
        object, 
        from = grep("Ms1Extracted", names(object), invert = TRUE),
        to = "Ms1Extracted",
        varFrom = rep("Precursor.Id", length(names(object)) - 1),
        varTo = "Precursor.Id"
    )
}

## Internal functions that subsets the extractedData to keep only
.keepSharedSamples <- function(extractedData, object) {
    cnames <- unique(unlist(colnames(object)))
    if (any(mis <- !cnames %in% colnames(extractedData)))
        stop("Some columns present in reportData are not found in ",
             "extracted data", paste0(cnames[mis], collapse = ", "),
             "\nAre you sure the two tables were generated from ",
             "the same experiment?")
    extractedData[, cnames]
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
