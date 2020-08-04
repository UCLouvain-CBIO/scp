##' @import SingleCellExperiment
##' @import QFeatures
##' @importFrom utils read.csv
##' @import tidyverse


##' @title Read single-cell proteomics data as a QFeatures object from tabular 
##' data and metadata
##'
##' @description
##'
##' Convert tabular quantitative MS data and metadata from a spreadsheet or a 
##' `data.frame` into a `QFeatures` object containing `SingleCellExperiment` 
##' objects. 
##'
##' @param quantTable File or object holding the quantitative data. Can be
##'     either a `character(1)` with the path to a text-based
##'     spreadsheet (comma-separated values by default, but see `...`)
##'     or an object that can be coerced to a `data.frame`. It is
##'     advised not to encode characters as factors.
##' 
##' @param metaTable A `data.frame` or any object that can be coerced to a 
##'     `data.frame`. `metaTable` is expected to contains all the sample meta 
##'     information. Required fields are the acquisition batch (given by 
##'     `batchCol`) and the acquisition channel within the batch (e.g. TMT 
##'     channel, given by `channelCol`). Additional fields (e.g. sample type, 
##'     acquisition date,...) are allowed and will be stored as sample meta data.
##'
##' @param batchCol A `numeric(1)` or `character(1)` pointing to the column of 
##'     `quantTable` and `metaTable` that contain the batch names. Make sure 
##'     that the column name in both table are either identical (if you supply a
##'     `character`) or have the same index (if you supply a `numeric`).
##' 
##' @param channelCol A `numeric(1)` or `character(1)` pointing to the column of
##'     `metaTable` that contains the column names of the quantitive data in 
##'     `quantTable` (see Example).
##' 
##' 
##' @param verbose A `logical(1)` indicating whether the progress of the data 
##'     reading and formatting should be printed to the console. 
##' 
##' @param ... Further arguments that can be passed on to [read.csv]
##'     except `stringsAsFactors`, which is always `FALSE`.
##'
##' @details 
##' 
##'
##' @return An instance of class [QFeatures]. The expression data of each batch 
##'     is stored in a separate assay as a `SingleCellExperiment` object.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##' 
##' @importFrom utils read.csv
##'
##' @md
##' @export
##' 
##' @examples 
##' TODO
readSCP <- function(quantTable, 
                    metaTable,
                    batchCol,
                    channelCol, 
                    verbose = TRUE,
                    ...) {
  metaTable <- as.data.frame(metaTable)
  ## Create the SingleCellExperiment object
  if (verbose) cat("Loading data as a 'SingleCellExperiment' object\n")
  ecol <- unique(metaTable[, channelCol])
  scp <- readSingleCellExperiment(table = quantTable, 
                                  ecol = ecol, 
                                  ...)
  if (is.null(list(...)$row.names)) 
    rownames(scp) <- paste0("PSM", 1:nrow(scp))
  ## Check the link between metaTable and scp
  mis <- !rowData(scp)[, batchCol] %in% metaTable[, batchCol]
  if (any(mis)) {
    warning("Missing metadata. The features are removed for ", 
            paste0(unique(rowData(scp)[mis, batchCol]), collapse = ", "))
    scp <- scp[!mis, ]
  }
  
  ## Split the SingleCellExperiment object by batch column
  if (verbose) cat(paste0("Splitting data based on '", batchCol, "'\n"))
  scp <- .splitSCE(scp, f = batchCol)
  
  ## Add unique sample identifiers
  if (verbose) cat(paste0("Formating sample metadata (colData)\n"))
  for (i in 1:length(scp)) {
    colnames(scp[[i]]) <- paste0(names(scp)[[i]], "_", colnames(scp[[i]]))
  }
  ## Create the colData 
  cd <- DataFrame(row.names = unlist(lapply(scp, colnames)))
  rownames(metaTable) <- paste0(metaTable[, batchCol], "_", 
                                metaTable[, channelCol])
  cd <- cbind(cd, metaTable[rownames(cd), ])
  
  ## Store the data as a QFeatures object and add the experimental information
  if (verbose) cat("Formating data as a 'QFeatures' object\n")
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
##' @export
readSingleCellExperiment <- function(table, 
                                     ecol, 
                                     fnames,
                                     ...){
  ## Read data as SummarizedExperiment
  sce <- readSummarizedExperiment(table, ecol, fnames, ...)
  sce <- as(sce, "SingleCellExperiment")
  return(sce)
}


# Split SingleCellExperiment into an ExperimentList
#
# The fonction creates an [ExperimentList] containing [SingleCellExperiment] 
# objects from a [SingleCellExperiment] object. `f` is used to split `x`` along 
# the rows (`f`` was a feature variable name) or samples/columns (f was a 
# phenotypic variable name). If f is passed as a factor, its length will be 
# matched to nrow(x) or ncol(x) (in that order) to determine if x will be split 
# along the features (rows) or sample (columns). Hence, the length of f must 
# match exactly to either dimension.
# 
# This function is not exported. if this is needed, create a pull request to 
# `rformassspectrometry/QFeatures`.
# 
# @param x a single [SingleCellExperiment] object 
# @param f a factor or a character of length 1. In the latter case, `f`` will 
#     be matched to the row and column data variable names (in that order). If 
#     a match is found, the respective variable is extracted, converted to a 
#     factor if needed
#
.splitSCE <- function(x, f){
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



##' Compute the sample over carrier ratio (SCR)
##' 
##' The function computes the ratio of the intensities of sample channels over 
##' the intentisty of the carrier channel for each feature. The ratios are 
##' averaged within the assay.
##'
##' @param obj A `QFeatures` object.
##' @param i A `character()` or `integer()` indicating for which assay(s) the
##'     SCR needs to be computed. 
##' @param colDataCol A `character(1)` indicating the variable to take from 
##'     `colData(obj)` that gives the sample annotation. 
##' @param samplePattern A `character(1)` pattern that matches the sample 
##'     encoding in `colDataCol`.
##' @param carrierPattern A `character(1)` pattern that matches the carrier 
##'     encoding in `colDataCol`. Only one match per assay is allowed, otherwise
##'     only the first match is taken
##' 
##' @return A `QFeatures` object for which the `rowData` of the given assay(s) is 
##'     augmented with the mean SCR (`meanSCR` variable).
##'
##' @export
##'
##' @examples
##' TODO
computeSCR <- function(obj, 
                       i, 
                       colDataCol, 
                       samplePattern, 
                       carrierPattern,
                       verbose = TRUE){
  if (!inherits(obj, "QFeatures")) stop("'obj' must be a QFeatures object")
  if (is.numeric(samplePattern)) 
    warning("The pattern is numeric. This is only allowed for replicating the SCoPE2 analysis and will later get defunct.\n")
    
  ## Iterate over the different assay indices
  for (ii in i) {
    annot <- colData(obj)[colnames(obj[[ii]]), ][, colDataCol]
    ## Get the corresponding indices
    if (is.numeric(samplePattern)) {
      sampIdx <- samplePattern
    } else {
      sampIdx <- grep(samplePattern, annot)
    }
    carrIdx <- grep(carrierPattern, annot)
    if (length(carrIdx) > 1) {
      warning("Multiple carriers found in assay '", names(obj)[ii], 
              "'. Only the first match will be used")
      carrIdx <- carrIdx[1]
    }
    ## Compute ratios
    ratio <- assay(obj[[ii]])[, sampIdx, drop = FALSE] / assay(obj[[ii]])[, carrIdx]
    ## Compute mean sample to carrier ratios
    rowData(obj@ExperimentList@listData[[ii]])$meanSCR <- 
      rowMeans(ratio, na.rm = TRUE)
    ## more efficient than rowData(obj[[ii]])$meanSCR <- rowMeans(ratio, na.rm = TRUE)
  }
  obj
}


##' Extract the `rowData` of a `QFeatures` object to a `DataFrame`
##'
##' The methods takes the `rowData` of one or more given assay in a `QFeatures`
##' object and combines the data in a single `DataFrame`. 
##' 
##' Along with the reuired `rowData` an additional `.assay` variable is created 
##' and holds the assay name from which the metadata was taken. 
##' 
##' @param obj A `QFeatures` object
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' @param vars A `character()` vector indicating which variables from the 
##'     `rowData` should be extracted.
##'
##' @return A `DataFrame` object with the `rowData` row-binded over the required
##'     assays.
##'     
##' @export
##'
##' @examples
##' TODO
rowDataToDF <- function(obj, i, vars) {
  if (!inherits(obj, "QFeatures")) stop("'obj' must be a QFeatures object")
  if (is.numeric(i)) i <- names(obj)[i]
  ## Make sure that the variables to extract are present in the rowData
  mis <- sapply(experiments(obj)[i], 
                function(x) any(!vars %in% colnames(rowData(x))))
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

##' Compute false discovery rates (FDRs) from posterior error probabilities (PEPs)
##' 
##' The functions takes the PEPs from the given assay's `rowData` and creates a 
##' new variable (`.FDR`) to it that contains the computed FDRs. 
##'
##' @param object A `QFeatures` object
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' @param groupCol A `character(1)` indicating the variable names in the 
##'     `rowData` that contains the grouping variable. The FDR are usually 
##'     computed for PSMs grouped by peptide ID. 
##' @param pepCol A `character(1)` indicating the variable names in the 
##'     `rowData` that contains the PEPs.
##'
##' @return A `QFeatures` object.
##' 
##' @export
##'
##' @examples
##' TODO
computeFDR <- function(object, 
                       i, 
                       groupCol, 
                       pepCol) {
  if (!inherits(object, "QFeatures")) stop("'object' must be a QFeatures object")
  if (is.numeric(i)) i <- names(object)[i]
  
  ## Function to compute FDRs from PEPs
  fdrFromPEP <- function(x) ## this is calc_fdr from SCoPE2
    return((cumsum(x[order(x)]) / seq_along(x))[order(order(x))])
  
  ## Get the PEP from all assays
  peps <- rowDataToDF(object, i, vars = c(groupCol, pepCol))
  colnames(peps)[1:2] <- c("groupCol", "pepCol")

  ## Order the features to replicate SCoPE2
  warning("Features are sorted. This is only needed when replicating the SCoPE2 analysis\nSee also https://github.com/SlavovLab/SCoPE2/issues/3")
  peps <- arrange(data.frame(peps), .assay, groupCol, pepCol)
  ## Compute the FDR for every peptide ID separately
  peps <- group_by(data.frame(peps), groupCol)
  peps <- mutate(peps, FDR = fdrFromPEP(pepCol))

  ## Insert the FDR inside every assay
  pepID <- paste0(peps$.assay, peps$.rowname)
  for (ii in i) {
    rdID <- paste0(ii, rownames(object[[ii]]))
    .FDR <- peps[match(rdID, pepID), ]$FDR
    rowData(object@ExperimentList@listData[[ii]])$.FDR <- .FDR
  }
  return(object)
}

##' Compute cell median coefficient of variation (CV)
##' 
##' The function computes for each cell the median CV. The expression data is 
##' normalized twice. First, cell median expression is used as normalization 
##' factor, then, the mean for each batch and peptide. The CV is then computed 
##' for each protein in each cell. CV is the standard deviation divided by the 
##' mean expression. The CV is computed only if there are more than 5 
##' observations per protein per cell. 
##' 
##' A new columns, `medianCV`, is added to the `colData` of the assay `i` and
##' contains the computed median CVs.
##' 
##' *Watch out* that `peptideCol` and `proteinCol` are feature variables and 
##' hence taken from the `rowData`. `batchCol` is a sample variable and is taken
##' from the `colData` of the `QFeatures` object.
##'
##' @param object A `QFeatures` object
##' @param i  A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' @param peptideCol  A `character(1)` indicating the variable name in the 
##'     `rowData` that contains the peptide grouping.
##' @param proteinCol A `character(1)` indicating the variable name in the 
##'     `rowData` that contains the protein grouping.
##' @param batchCol A `character(1)` indicating the variable name in the 
##'     `colData` of `object` that contains the batch names.
##'     
##' @return A `QFeatures` object. 
##' 
##' @export
##'
##' @examples
##' TODO
computeMedianCV <- function(object, 
                            i, 
                            peptideCol, 
                            proteinCol,
                            batchCol) {
  warning("Cell type selection is performed to reproduce SCoPE2. This should be removed!")
  ## Extract the expression data and metadata as long format
  colname = "colname" ## for compatibility with `group_by_at`
  object %>%
    .assayToLongDF(i = "peptides_filter1", 
                   rowDataCols = c(peptideCol, proteinCol), 
                   colDataCols = c(batchCol, "SampleType")) %>%
    data.frame %>%
    ## Normalize cells with median
    group_by(colname) %>%
    mutate(norm_q1 = value / median(value, na.rm = TRUE)) %>%
    ## Normalized peptides/Set with mean of cell normalized expression
    group_by_at(c(peptideCol, batchCol)) %>%
    mutate(norm_q = value / mean(norm_q1, na.rm = TRUE)) %>%
    ## Filter cell type. 
    ## TODO remove this filtering
    filter(SampleType %in% c("Macrophage", "Monocyte", "Blank")) %>%
    ## Compute the protein CV in every cell
    group_by_at(c(proteinCol, colname)) %>%
    mutate(norm_q_sd = sd(norm_q, na.rm = TRUE),
           norm_q_mean = mean(norm_q, na.rm = TRUE),
           cvq = norm_q_sd / norm_q_mean) %>%
    ## Remove CVs that were computed based on few data points
    group_by_at(c(proteinCol, colname)) %>%
    mutate(cvn = sum(!is.na(norm_q))) %>%
    filter(cvn > 5) %>%
    ## Compute the median CV per cell
    group_by(colname) %>%
    mutate(MedianCV = median(cvq, na.rm = TRUE)) %>%
    ## Store the cell median CV in the coldata
    select(colname, MedianCV) %>%
    unique ->
    CVs
  object@ExperimentList@listData[[i]]$MedianCV <- NA
  colData(object@ExperimentList@listData[[i]])[CVs$colname, "MedianCV"] <- 
    CVs$MedianCV
  return(object)
}

## Internal function to efficiently extract expression data to long format. The 
## efficiency is seen when nNA's are present and `na.rm == TRUE`. Meta
.assayToLongDF <- function(obj, colDataCols, rowDataCols, i, na.rm = TRUE) {
  if (length(i) > 1) stop("Multiple assays are not supported (yet).")
  dat <- assay(obj[[i]])
  if (na.rm) sel <- which(!is.na(dat), arr.ind = TRUE) else sel <- TRUE
  DataFrame(colname = colnames(dat)[sel[, 2]],
            rowname = rownames(dat)[sel[, 1]],
            colData(obj)[colnames(dat)[sel[, 2]], colDataCols, drop = FALSE],
            rowData(obj[[i]])[rownames(dat)[sel[, 1]], rowDataCols, drop = FALSE],
            value = dat[sel])
}


## Internal function to efficiently replace a `SummarizedExperiment` (and 
## subclasses of it) assay in a `QFeatures` object. Watch out this function 
## should only be used when the old assay and the assay to replace have same 
## size. Otherwise, this can lead to inconsistencies in the colData and 
## functionality from `MultiAssayExperiment` should be used. 
## obj: A QFeatures object
## y: A SummarizedExperiment object, or an object that inherits from it
## i: a characther(1) or logical(1) that indicates which assay must be replaced
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
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' @param colDataCol A `character(1)` indicating the variable to take from 
##'     `colData(obj)` that gives the sample annotation. 
##' @param samplePattern A `character(1)` pattern that matches the sample 
##'     encoding in `colDataCol`
##' @param refPattern A `character(1)` pattern that matches the carrier 
##'     encoding in `colDataCol`. Only one match per assay is allowed, otherwise
##'     only the first match is taken
##'
##' @return A `QFeatures` object
##' 
##' @export
##'
##' @examples
##' TODO
divideByReference <- function(obj, 
                              i,
                              colDataCol,
                              samplePattern,
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
##' This function coinsiders any infinite value as missing data. So, any value
##' in assay `i` that return `TRUE` to `is.infinite` will be replaced by `NA`.
##' 
##' @param obj A `QFeatures` object
##' @param i A `numeric()` or `character()` vector indicating from which assays 
##'     the `rowData` should be taken.
##' 
##' @return A `QFeatures` object
##' 
##' @export
##' 
##' @examples 
##' TODO
infIsNA <- function(obj, i) {
  for (ii in i) {
    sel <- is.infinite(assay(obj[[ii]])) 
    assay(obj[[ii]])[sel] <- NA
  }
  obj
}

##' Transfer the `colData` to an Assay
##' 
##' The function transfers the `colData` from a `QFeatures` object to one of the 
##' assays it contains. The transfered data is bound to the existing `colData` 
##' of the target assay.
##'
##' @param obj A `QFeatures` object
##' @param i A `numeric(1)` or `character(1)` indicating which assay to transfer 
##' the `colData` to.
##'
##' @return A `QFeatures` object
##' 
##' @export
##'
##' @examples
##' TODO
transferColDataToAssay <- function (obj, i) {
  cd <- colData(obj)[colnames(obj[[i]]),  ]
  if (all(colnames(cd) %in% colnames(colData(obj@ExperimentList[[i]])))) {
    message("The colData is already present in assay ", i, ".")
    return(obj)
  }
  colData(obj@ExperimentList[[i]]) <- cbind(colData(obj[[i]]), cd)
  return(obj)
}

##' Aggregate features over multiple assays
##' 
##' This function is a wrapper function around [QFeatures::aggregateFeatures]. It
##' allows the user to provide multiple assays for which `aggregateFeatures` 
##' will be applied sequentially. 
##' 
##' @param obj A `QFeatures` object
##' @param i A `numeric(1)` or `character(1)` indicating which assay to transfer 
##'     the `colData` to.
##' @param fcol The feature variables for each assays `i` defining how to 
##'     summarise the QFeatures. If `fcol` has length 1, the variable name is 
##'     assumed to be the same for all assays
##' @param name A `character()` naming the new assay. `name` must have the same
##'     length as `i`. Note that the function will fail if of the names in 
##'     `name` is already present. 
##' @param fun A function used for quantitative feature aggregation. 
##' @param ... Additional parameters passed the `fun`.
##' 
##' @return A `QFeatures` object 
##' 
##' @export
##' 
##' @seealso [QFeatures::aggregateFeatures]
##' 
##' @examples 
##' TODO
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
  el <- obj@ExperimentList[i]
  for (j in seq_along(el)) {
    suppressMessages(
      el[[j]] <- aggregateFeatures(el[[j]], fcol = fcol[j], fun = fun, ...)
    )
    ## Print progress
    cat(paste0("\rAggregated: ", j, "/", length(el)))
    if (j == length(el)) cat ("\n")
    flush.console()
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

