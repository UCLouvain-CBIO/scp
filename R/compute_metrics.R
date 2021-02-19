
####---- Internal functions ----####

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
    ## Return NA if nna < nobs. This behaviour is similar to sd
    nna[nna < nobs] <- NA
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

## Internal function: compute q-values from a vector of PEPs
## this is calc_fdr from SCoPE2
##  x: numeric() containing probabilities 
.pep2qvalue <- function(x) {
    (cumsum(x[order(x)]) / seq_along(x))[order(order(x))]
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
##' @param rowDataName A `character(1)` giving the name of the new 
##'      variable in the `rowData` where the computed SCR will be 
##'      stored. The name cannot already exist in any of the assay
##'      `rowData`.
##'     
##' @return A `QFeatures` object for which the `rowData` of the given
##'     assay(s) is augmented with the mean SCR. 
##'
##' @export
##'
##' @examples
##' data("scp1")
##' scp1 <- computeSCR(scp1, 
##'                    i = 1,
##'                    colDataCol = "SampleType",
##'                    carrierPattern = "Carrier",
##'                    samplePattern = "Blank|Macrophage|Monocyte",
##'                    rowDataName = "MeanSCR")
##' ## Check results
##' rowDataToDF(scp1, 1, "MeanSCR")
##' 
computeSCR <- function(object, 
                       i, 
                       colDataCol, 
                       samplePattern, 
                       carrierPattern,
                       rowDataName = "MeanSCR") {
    if (!inherits(object, "QFeatures"))
        stop("'object' must be a QFeatures object")
    if (is.numeric(samplePattern)) 
        warning("The pattern is numeric. This is only allowed for replicating ",
                "the SCoPE2 analysis and will later get defunct.")
    if (any(unlist(rowDataNames(object)[i]) == rowDataName)) 
        stop("The rowData name '", rowDataName, "' already exists. ", 
             "Use another name or remove that column before running ",
             "this function.")
    
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
        rowData(object@ExperimentList@listData[[ii]])[, rowDataName] <- 
            rowMeans(ratio, na.rm = TRUE)
        ## more efficient than rowData(object[[ii]])[, rowDataName] <-  ...
    }
    object
}

##' Compute q-values
##' 
##' This function computes q-values from the posterior error 
##' probabilities (PEPs). The functions takes the PEPs from the given
##' assay's `rowData` and adds a new variable to it that contains the
##' computed q-values.
##'
##' @param object A `QFeatures` object
##' 
##' @param i A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##' 
##' @param groupBy A `character(1)` indicating the variable name in 
##'     the `rowData` that contains the grouping variable, for 
##'     instance to compute protein FDR. When `groupBy` is not missing,
##'     the best feature approach is used to compute the PEP per group,
##'     meaning that the smallest PEP is taken as the PEP of the group. 
##' 
##' @param PEP A `character(1)` indicating the variable names in the 
##'     `rowData` that contains the PEPs. Since, PEPs are probabilities, the 
##'     variable must be contained in (0, 1).
##'
##' @param rowDataName A `character(1)` giving the name of the new 
##'      variable in the `rowData` where the computed FDRs will be 
##'      stored. The name cannot already exist in any of the assay
##'      `rowData`.
##'
##' @return A `QFeatures` object.
##' 
##' @details 
##' 
##' The q-value of a feature (PSM, peptide, protein) is the minimum 
##' FDR at which that feature will be selected upon filtering 
##' (Savitski et al.). On the other hand, the feature PEP is the 
##' probability that the feature is wrongly matched and hence can be 
##' seen as a local FDR (Kall et al.). While filtering on PEP is 
##' guaranteed to control for FDR, it is usually too conservative. 
##' Therefore, we provide this function to convert PEP to q-values. 
##' 
##' We compute the q-value of a feature as the average of the PEPs 
##' associated to PSMs that have equal or greater identification 
##' confidence (so smaller PEP). See Kall et al. for a visual 
##' interpretation.
##' 
##' We also allow inference of q-values at higher level, for instance 
##' computing the protein q-values from PSM PEP. This can be performed 
##' by supplying the `groupBy` argument. In this case, we adopt the 
##' best feature strategy that will take the best (smallest) PEP for 
##' each group (Savitski et al.).
##' 
##' @references 
##' 
##' Käll, Lukas, John D. Storey, Michael J. MacCoss, and William 
##' Stafford Noble. 2008. “Posterior Error Probabilities and False 
##' Discovery Rates: Two Sides of the Same Coin.” Journal of Proteome
##' Research 7 (1): 40–44.
##' 
##' Savitski, Mikhail M., Mathias Wilhelm, Hannes Hahne, Bernhard 
##' Kuster, and Marcus Bantscheff. 2015. “A Scalable Approach for 
##' Protein False Discovery Rate Estimation in Large Proteomic Data 
##' Sets.” Molecular & Cellular Proteomics: MCP 14 (9): 2394–2404.
##' 
##' @export
##'
##' @examples
##' data("scp1")
##' scp1 <- pep2qvalue(scp1,
##'                    i = 1,
##'                    groupBy = "protein",
##'                    PEP = "dart_PEP",
##'                    rowDataName = "qvalue_protein")
##' ## Check results
##' rowDataToDF(scp1, 1, c("dart_PEP", "qvalue_protein"))
##' 
pep2qvalue <- function(object, 
                       i, 
                       groupBy, 
                       PEP,
                       rowDataName = "qvalue") {
    if (!inherits(object, "QFeatures"))
        stop("'object' must be a QFeatures object")
    if (is.numeric(i)) i <- names(object)[i]
    ## Check arguments: rowDataName is valid
    if (rowDataName %in% unlist(rowDataNames(object)[i]))
        stop("The rowData name '", rowDataName, "' already exists. ", 
             "Use another name or remove that column before running ",
             "this function.")
    
    
    ## Get the PEP from all assays
    vars <- PEP
    if (!missing(groupBy)) vars <- c(vars, groupBy)
    df <- rowDataToDF(object, i, vars = vars)
    
    ## Check PEP is a probability
    pepRange <- range(df[, PEP], na.rm = TRUE)
    if (max(pepRange) > 1 | min(pepRange < 0))
        stop(paste0("'", PEP, "' is not a probability in (0, 1)"))
    
    ## Compute the q-values 
    if (missing(groupBy)) {
        df$qval <- .pep2qvalue(df[, PEP])
    } else { ## Apply grouping if supplied
        p <- split(df[[PEP]], df[[groupBy]])
        p <- sapply(p, min, na.rm = TRUE)
        qval <- .pep2qvalue(p)
        df$qval <- unname(qval[df[[groupBy]]])
    }
    
    ## Insert the q-value inside every assay rowData
    pepID <- paste0(df$.assay, df$.rowname)
    for (ii in i) {
        rdID <- paste0(ii, rownames(object[[ii]]))
        rowData(object@ExperimentList@listData[[ii]])[, rowDataName] <- 
            df[match(rdID, pepID), ]$qval
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
##' @param nobs An `integer(1)` indicating how many observations 
##'     (features) should at least be considered for computing the CV.
##'     Since no CV can be computed for less than 2 observations, 
##'     `nobs` should at least be 2. 
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
##'     [MsCoreUtils::normalizeMethods]. You can also specify
##'     `norm = "SCoPE2"` to reproduce the normalization performed 
##'     before computing the CVs as suggested by Specht et al. 
##'     `norm = "none"` will not normalize the data (default)
##' 
##' @param ... Additional arguments that are passed to the 
##'     normalization method. 
##'     
##' @return A `QFeatures` object. 
##' 
##' @references Specht, Harrison, Edward Emmott, Aleksandra A. Petelski,
##'     R. Gray Huffman, David H. Perlman, Marco Serra, Peter Kharchenko, 
##'     Antonius Koller, and Nikolai Slavov. 2021. “Single-Cell Proteomic
##'      and Transcriptomic Analysis of Macrophage Heterogeneity Using 
##'      SCoPE2.” Genome Biology 22 (1): 50.
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
##' @rdname medianCVperCell
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
    }
    ## Warn the median CV cannot be computed for some cells
    if (any(is.na(medCVs)))
        warning("The median CV could not be computed for one or more ",
                "samples. You may want to try a smaller value for ",
                "'nobs'.\n")
    ## Store the medians CVs in the colData
    colData(object)[names(medCVs), colDataName] <- medCVs
    return(object)
}

## TODO discuss whether to export this??
## 
## @param x A `SingleCellExperiment` object
##     
## @param group A `factor()` that indicates how features (rows) should
##    be grouped. The CVs are computed for each group separately. 
##      
## @rdname medianCVperCell
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
    if (identical(norm, "SCoPE2")) {
        xnorm <- .normalizeSCP(x, method = "div.median")
        assay(x) <- sweep(assay(x), 1, rowMeans(assay(xnorm), na.rm = TRUE), "/")
    } else if (!identical(norm, "none")) {
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
##' This function is deprecated and should no longer be used. To 
##' reproduce the SCoPE2 script, you can now use `medianCVperCell`
##' with the following arguments:
##' 
##' - `norm = "SCoPE2"`
##' - `nobs = 6`
##' 
##' Make sure to provide the peptide data from separate assays so that 
##' the normalization factors are computed per batch.
##' 
##' @param object NULL
##' @param i NULL
##' @param peptideCol NULL
##' @param proteinCol NULL
##' @param batchCol NULL
##'
##' @export
computeMedianCV_SCoPE2 <- function(object, 
                                   i, 
                                   peptideCol, 
                                   proteinCol, 
                                   batchCol) {
    stop("'computeMedianCV_SCoPE2' is deprecated and should no longer be used.\n",
            "To reproduce the SCoPE2 script, you can now use ",
            "'medianCVperCell' with the following arguments:\n",
            " - 'norm' = \"SCoPE2\"\n",
            " - 'nobs' = 6\n",
            "Make sure to provide the peptide data from separate assays ", 
            "so that the normalization factors are computed per batch.\n")
}
