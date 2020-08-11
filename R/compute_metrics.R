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
##' scp1 <- computeSCR(scp1, 
##'                    i = 1,
##'                    colDataCol = "SampleType",
##'                    carrierPattern = "Carrier",
##'                    samplePattern = "Blank|Macrophage|Monocyte")
##' ## Check results
##' rowDataToDF(scp1, 1, "meanSCR")
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
    if (any(!c(length(carrIdx), length(sampIdx)))) stop("Pattern did not match a sample or carrier channel.")
    ## Compute ratios
    ratio <- assay(obj[[ii]])[, sampIdx, drop = FALSE] / assay(obj[[ii]])[, carrIdx]
    ## Compute mean sample to carrier ratios
    rowData(obj@ExperimentList@listData[[ii]])$meanSCR <- 
      rowMeans(ratio, na.rm = TRUE)
    ## more efficient than rowData(obj[[ii]])$meanSCR <- rowMeans(ratio, na.rm = TRUE)
  }
  obj
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
##'     `rowData` that contains the PEPs. Since, PEPs are probabilities, the 
##'     variable must be contained in (0, 1).
##'
##' @return A `QFeatures` object.
##' 
##' @export
##'
##' @examples
##' scp1 <- computeFDR(scp1,
##'                    i = 1,
##'                    groupCol = "Sequence",
##'                    pepCol = "PEP")
##' ## Check results
##' rowDataToDF(scp1, 1, c("PEP", ".FDR"))
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
  if (max(peps$pepCol) > 1 | min(peps$pepCol) < 0) 
    stop(paste0("'", pepCol, "' must link to a probability in (0, 1)"))
  
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
##' scp1 <- computeMedianCV(scp1, 
##'                         i = "peptides", 
##'                         proteinCol = "protein", 
##'                         peptideCol = "peptide", 
##'                         batchCol = "Set")
##' ## Check results
##' hist(scp1[["peptides"]]$MedianCV)
##' 
computeMedianCV <- function(object, 
                            i, 
                            peptideCol, 
                            proteinCol,
                            batchCol) {
  warning("Cell type selection is performed to reproduce SCoPE2. This should be removed!")
  ## Extract the expression data and metadata as long format
  colname = "colname" ## for compatibility with `group_by_at`
  object %>%
    .assayToLongDF(i = i, 
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
