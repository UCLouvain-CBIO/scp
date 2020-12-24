
####---- Internal functions ----####


## Internal function to efficiently extract expression data to long
## format. The efficiency is seen when nNA's are present and `na.rm ==
## TRUE`. Meta
.assayToLongDF <- function(object, 
                           colDataCols, 
                           rowDataCols, 
                           i, 
                           na.rm = TRUE) {
    if (length(i) > 1) stop("Multiple assays are not supported (yet).")
    dat <- assay(object[[i]])
    if (na.rm) sel <- which(!is.na(dat), arr.ind = TRUE) else sel <- TRUE
    DataFrame(colname = colnames(dat)[sel[, 2]],
              rowname = rownames(dat)[sel[, 1]],
              colData(object)[colnames(dat)[sel[, 2]], 
                              colDataCols, 
                              drop = FALSE],
              rowData(object[[i]])[rownames(dat)[sel[, 1]], 
                                   rowDataCols, 
                                   drop = FALSE],
              value = dat[sel])
}

## Internal function: compute coefficient of variation for each column, 
## taking into account the grouping structure of the rows
## This code is inspired from MSnbase::rowsd
##' @importFrom stats median sd
.rowCV <- function(x, 
                   group, 
                   nobs = 2,
                   reorder = FALSE, 
                   na.rm = FALSE) {
    ## Check nobs
    if (nobs < 2)
       stop("'nobs' must be at least 2. No sd can be computed for ",
            "less than 2 observation.")
    if (na.rm) {
        nna <- !is.na(x)
        mode(nna) <- "numeric"
    } else {
        nna <- matrix(1, ncol = ncol(x), nrow = nrow(x))
    }
    ## Get the number of non-NA observations per column per group
    nna <- rowsum(nna, 
                  group = group, 
                  reorder = reorder, 
                  na.rm = na.rm)
    nna[nna < nobs] <- NA ## return NA if nna <= 1 (similar to sd)
    ## Compute mean
    mean <- rowsum(x, 
                   group = group, 
                   reorder = reorder, 
                   na.rm = na.rm) / nna
    ## Compute variance
    var <- rowsum(x * x, 
                  group = group, 
                  reorder = reorder, 
                  na.rm = na.rm) / nna - mean^2L
    ## CV = sd / mean
    sqrt(var * nna / (nna - 1L)) / mean 
}


####---- Exported functions ----####


##' Compute the sample over carrier ratio (SCR)
##' 
##' The function computes the ratio of the intensities of sample
##' channels over the intentisty of the carrier channel for each
##' feature. The ratios are averaged within the assay.
##'
##' @param object A `QFeatures` object.
##' 
##' @param i A `character()` or `integer()` indicating for which
##'     assay(s) the SCR needs to be computed.
##' 
##' @param colDataCol A `character(1)` indicating the variable to take
##'     from `colData(object)` that gives the sample annotation.
##' 
##' @param samplePattern A `character(1)` pattern that matches the
##'     sample encoding in `colDataCol`.
##' 
##' @param carrierPattern A `character(1)` pattern that matches the
##'     carrier encoding in `colDataCol`. Only one match per assay is
##'     allowed, otherwise only the first match is taken
##'
##' @return A `QFeatures` object for which the `rowData` of the given
##'     assay(s) is augmented with the mean SCR (`.meanSCR`
##'     variable). 
##'
##' @export
##'
##' @examples
##' data("scp1")
##' scp1 <- computeSCR(scp1, 
##'                    i = 1,
##'                    colDataCol = "SampleType",
##'                    carrierPattern = "Carrier",
##'                    samplePattern = "Blank|Macrophage|Monocyte")
##' ## Check results
##' rowDataToDF(scp1, 1, ".meanSCR")
##' 
computeSCR <- function(object, 
                       i, 
                       colDataCol, 
                       samplePattern, 
                       carrierPattern) {
    if (!inherits(object, "QFeatures"))
        stop("'object' must be a QFeatures object")
    if (is.numeric(samplePattern)) 
        warning("The pattern is numeric. This is only allowed for replicating ",
                "the SCoPE2 analysis and will later get defunct.")
    
    ## Iterate over the different assay indices
    for (ii in i) {
        annot <- colData(object)[colnames(object[[ii]]), ][, colDataCol]
        ## Get the corresponding indices
        if (is.numeric(samplePattern)) {
            sampIdx <- samplePattern[samplePattern <= length(annot)]
        } else {
            sampIdx <- grep(samplePattern, annot)
        }
        carrIdx <- grep(carrierPattern, annot)
        if (length(carrIdx) > 1) {
            warning("Multiple carriers found in assay '", names(object)[ii], 
                    "'. Only the first match will be used")
            carrIdx <- carrIdx[1]
        } 
        if (any(!c(length(carrIdx), length(sampIdx))))
            stop("Pattern did not match a sample or carrier channel.")
        ## Compute ratios
        carrier <- assay(object[[ii]])[, carrIdx]
        ratio <- assay(object[[ii]])[, sampIdx, drop = FALSE] / carrier
        ## Compute mean sample to carrier ratios
        rowData(object@ExperimentList@listData[[ii]])$.meanSCR <- 
            rowMeans(ratio, na.rm = TRUE)
        ## more efficient than rowData(object[[ii]])$.meanSCR <-  ...
    }
    object
}

##' Compute FDR from posterior error probabilities PEP
##' 
##' The functions takes the posterior error probabilities (PEPs) from 
##' the given assay's `rowData` and adds a new variable to it (called 
##' `.FDR`) that contains the computed false discovery rates (FDRs). 
##'
##' @param object A `QFeatures` object
##' 
##' @param i A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##' 
##' @param groupBy A `character(1)` indicating the variable names in the 
##'     `rowData` that contains the grouping variable. The FDR are usually 
##'     computed for PSMs grouped by peptide ID.
##' 
##' @param PEP A `character(1)` indicating the variable names in the 
##'     `rowData` that contains the PEPs. Since, PEPs are probabilities, the 
##'     variable must be contained in (0, 1).
##'
##' @param colDataName A `character(1)` giving the name of the new 
##'      variable in the `colData` where the computed FDRs will be 
##'      stored. The name cannot already exist in the `colData`.
##'
##' @return A `QFeatures` object.
##' 
##' @export
##'
##' @examples
##' data("scp1")
##' scp1 <- computeFDR(scp1,
##'                    i = 1,
##'                    groupBy = "Sequence",
##'                    PEP = "dart_PEP",
##'                    colDataName = "peptideFDR")
##' ## Check results
##' rowDataToDF(scp1, 1, c("dart_PEP", "peptideFDR"))
##' 
computeFDR <- function(object, 
                       i, 
                       groupBy, 
                       PEP,
                       colDataName = "FDR") {
    if (!inherits(object, "QFeatures"))
        stop("'object' must be a QFeatures object")
    if (is.numeric(i)) i <- names(object)[i]
    ## Check arguments: colDataName is valid
    if (colDataName %in% colnames(colData(object)))
        stop("The colData name '", colDataName, "' already exists. ", 
             "Use another name or remove that column before running ",
             "this function.")
    
    ## Function to compute FDRs from PEPs
    fdrFromPEP <- function(x) ## this is calc_fdr from SCoPE2
        return((cumsum(x[order(x)]) / seq_along(x))[order(order(x))])
    
    ## Get the PEP from all assays
    peps <- rowDataToDF(object, i, vars = c(groupBy, PEP))
    
    ## Check PEP is a probability
    pepRange <- range(peps[, PEP], na.rm = TRUE)
    if (max(pepRange) > 1 | min(pepRange < 0))
        stop(paste0("'", PEP, "' is not a probability in (0, 1)"))
    
    ## Report missing values
    if (anyNA(peps[, PEP]))
        message("The 'PEP' contains missing values. No FDR will ",
                "be computed for missing data.")
    
    ## Compute the FDR for every peptide ID separately
    peps <- group_by(data.frame(peps), .data[[groupBy]])
    peps <- mutate(peps, FDR = fdrFromPEP(.data[[PEP]]))
    
    ## Insert the FDR inside every assay
    pepID <- paste0(peps$.assay, peps$.rowname)
    for (ii in i) {
        rdID <- paste0(ii, rownames(object[[ii]]))
        .FDR <- peps[match(rdID, pepID), ]$FDR
        rowData(object@ExperimentList@listData[[ii]])[, colDataName] <- .FDR
    }
    return(object)
}

##' Compute the median coefficient of variation (CV) per cell
##' 
##' The function computes for each cell the median CV and stores them 
##' accordingly in the `colData` of the `QFeatures` object. The CVs in
##' each cell are computed from a group of features. The grouping is 
##' defined by a variable in the `rowData`. The function can be 
##' applied to one or more assays, as long as the samples (column 
##' names) are not duplicated. Also, the user can supply a minimal 
##' number of observations required to compute a CV to avoid that CVs
##' computed on too few observations influence the distribution within 
##' a cell. The quantification matrix can be optionally normalized 
##' before computing the CVs. Multiple normalizations are possible.
##' 
##' A new column is added to the `colData` of the object. The samples 
##' (columns) that are not present in the selection `i` will get 
##' assigned an NA.
##' 
##' @param object A `QFeatures` object
##'
##' @param i  A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##' 
##' @param groupBy  A `character(1)` indicating the variable name in 
##'     the `rowData` that contains the feature grouping.
##' 
##' @param nobs An `integer(1)` indicating how many observations should
##'     at least be considered for computing the CV. Since no CV can 
##'     be computed for less than 2 observations, `nobs` should at 
##'     least be 2. 
##' 
##' @param na.rm A `logical(1)` indicating whether missing data should
##'      be removed before computation.
##'      
##' @param colDataName A `character(1)` giving the name of the new 
##'      variable in the `colData` where the computed CVs will be 
##'      stored. The name cannot already exist in the `colData`.
##' 
##' @param norm A `character()` of normalization methods that will be 
##'     sequentially applied. Available methods and additional 
##'     information about normalization can be found in 
##'     [MsCoreUtils::normalizeMethods]. `norm = "none"` will not 
##'     normalize the data (default).
##' 
##' @param ... Additional arguments that are passed to the 
##'     normalization method. 
##'     
##' @return A `QFeatures` object. 
##' 
##' @export
##'
##' @importFrom matrixStats colMedians
##'
##' @examples
##' data("scp1")
##' scp1 <- filterFeatures(scp1, ~ !is.na(Proteins))
##' scp1 <- medianCVperCell(scp1, 
##'                         i = 1:3,
##'                         groupBy = "Proteins",
##'                         nobs = 5,
##'                         na.rm = TRUE,
##'                         colDataName = "MedianCV",
##'                         norm = "div.median")
##' ## Check results
##' hist(scp1$MedianCV)
##' 
medianCVperCell <- function(object,
                            i,
                            groupBy,
                            nobs = 5,
                            na.rm = TRUE,
                            colDataName = "MedianCV",
                            norm = "none",
                            ...) {
    ## Check arguments: object
    if (!inherits(object, "QFeatures"))
        stop("'object' must be a QFeatures object")
    ## Check arguments: index
    if (is.numeric(i)) i <- names(object)[i]
    ## Check arguments: colDataName is valid
    if (colDataName %in% colnames(colData(object)))
        stop("The colData name '", colDataName, "' already exists. ", 
             "Use another name or remove that column before running ",
             "this function.")
    ## Check arguments: no redundant columns
    coln <- unlist(colnames(object)[i])
    if (any(duplicated(coln))) 
        stop("Duplicated samples were found in assay(s) 'i'. This would ", 
             "lead to inconsistencies in the 'colData'.")
    
    ## Initiate the vectors with cell median CVs
    medCVs <- rep(NA, length(coln))
    names(medCVs) <- coln
    ## For each assay
    for (ii in i) {
        ## Compute the CV matrix
        cvs <- featureCV(object[[ii]],
                         group = as.factor(rowData(object[[ii]])[, groupBy]),
                         na.rm = na.rm,
                         norm = norm,
                         nobs = nobs,
                         ...)
        ## Compute the median CV per sample
        medCVs[colnames(cvs)] <- colMedians(cvs, na.rm = TRUE)
        ## Warn when nobs is too high
        if (any(is.na(medCVs[colnames(cvs)])))
            warning("The median CV is NA for at least one column in ",
                    "assay ", ii, ". Try a smaller value for 'nobs'.")
    }
    ## Store the medians CVs in the colData
    colData(object)[names(medCVs), colDataName] <- medCVs
    return(object)
}

## TODO discuss whether to export this??
featureCV <- function(x, 
                      group, 
                      na.rm = TRUE,
                      norm = "none",
                      nobs = 2,
                      ...) {
    ## Check object
    if (!inherits(x, "SingleCellExperiment"))
        stop("'x' must inherit from a 'SingleCellExperiment'")
    ## Optional normalization(s)   
    if (!identical(norm, "none")) {
        for(normi in norm)
            x <- .normalizeSCP(x, method = normi, ...)
    }
    ## Compute CVs
    .rowCV(assay(x), 
           group = group, 
           nobs = nobs,
           reorder = TRUE, 
           na.rm = na.rm)
}

##' (Deprecated) Compute the median coefficient of variation (CV) per 
##' cell
##' 
##' This function was implemented using the code provide in the SCoPE2
##' analysis (Specht et al. 2020). The function computes for each cell
##' the median CV. The expression data is normalized twice. First, 
##' cell median expression is used as normalization factor, then, the
##' mean for each batch and peptide. The CV is then computed for each
##' protein in each cell. CV is the standard deviation divided by the 
##' mean expression. The CV is computed only if there are more than 5 
##' observations per protein per cell. 
##' 
##' A new columns, `.medianCV`, is added to the `colData` of the assay 
##' `i` and contains the computed median CVs.
##' 
##' *Watch out* that `peptideCol` and `proteinCol` are feature 
##' variables and hence taken from the `rowData`. `batchCol` is a 
##' sample variable and is taken from the `colData` of the `QFeatures` 
##' object.
##'
##' @param object A `QFeatures` object
##'
##' @param i  A `numeric()` or `character()` vector indicating from which 
##'     assays the `rowData` should be taken.
##' 
##' @param peptideCol  A `character(1)` indicating the variable name in the 
##'     `rowData` that contains the peptide grouping.
##' 
##' @param proteinCol A `character(1)` indicating the variable name in the 
##'     `rowData` that contains the protein grouping.
##' 
##' @param batchCol A `character(1)` indicating the variable name in the 
##'     `colData` of `object` that contains the batch names.
##'     
##' @return A `QFeatures` object. 
##' 
##' @export
##'
##' @importFrom stats median sd
##' @importFrom rlang .data
##'
##' @examples
##' data("scp1")
##' scp1 <- computeMedianCV_SCoPE2(scp1,
##'                                i = "peptides",
##'                                proteinCol = "protein",
##'                                peptideCol = "peptide",
##'                                batchCol = "Set")
##' ## Check results
##' hist(scp1[["peptides"]]$.MedianCV)
##' 
computeMedianCV_SCoPE2 <- function(object, 
                                   i, 
                                   peptideCol, 
                                   proteinCol, 
                                   batchCol) {
    warning("This function is deprecated and is only present in the ",
            "package for replicating the SCoPE2 analysis.")
    ## Extract the expression data and metadata as long format
    object %>%
        .assayToLongDF(i = i, 
                       rowDataCols = c(peptideCol, proteinCol), 
                       colDataCols = c(batchCol)) %>%
        data.frame %>%
        ## Normalize cells with median
        group_by(.data$colname) %>%
        mutate(norm_q1 = .data$value / median(.data$value, na.rm = TRUE)) %>%
        ## Normalize peptides per Set with mean of cell normalized expression
        group_by(.data[[peptideCol]], .data[[batchCol]]) %>%
        mutate(norm_q = .data$value / mean(.data$norm_q1, na.rm = TRUE)) %>%
        ## Compute the protein CV in every cell
        group_by(.data[[proteinCol]], .data$colname) %>%
        mutate(norm_q_sd = sd(.data$norm_q, na.rm = TRUE),
               norm_q_mean = mean(.data$norm_q, na.rm = TRUE),
               cvq = .data$norm_q_sd / .data$norm_q_mean) %>%
        ## Remove CVs that were computed based on few data points
        group_by(.data[[proteinCol]], .data$colname) %>%
        mutate(cvn = sum(!is.na(.data$norm_q))) %>%
        dplyr::filter(.data$cvn > 5) %>%
        ## Compute the median CV per cell
        group_by(.data$colname) %>%
        mutate(.MedianCV = median(.data$cvq, na.rm = TRUE)) %>%
        ## Store the cell median CV in the colData
        select(.data$colname, .data$.MedianCV) %>%
        unique ->
        CVs
    object@ExperimentList@listData[[i]]$.MedianCV <- NA
    colData(object@ExperimentList@listData[[i]])[CVs$colname, ".MedianCV"] <- 
        CVs$.MedianCV
    return(object)
}
