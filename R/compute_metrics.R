
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

## @param x A `SummarizedExperiment` object
##     
## @param group A `factor()` that indicates how features (rows) should
##    be grouped. The CVs are computed for each group separately. 
##      
## @rdname medianCVperCell
featureCV <- function(x, group, na.rm = TRUE, norm = "none", nobs = 2, ...) {
    ## Check object
    if (!inherits(x, "SummarizedExperiment"))
        stop("'x' must inherit from a 'SummarizedExperiment'")
    ## Optional normalization(s)   
    if (identical(norm, "SCoPE2")) {
        xnorm <- .normalizeSCP(x, method = "div.median")
        assay(x) <- sweep(assay(x), 1, 
                          rowMeans(assay(xnorm), na.rm = TRUE), "/")
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

.getDetectionMatrix <- function(object, i) {
    i <- QFeatures:::.normIndex(object, i)
    if (length(i) > 1) stop("Invalid 'i'. You selected multiple assays.")
    x <- zeroIsNA(object[[i]])
    !is.na(assay(x))
}

## Internal function that computes missing value metrics given an
## identification matrix. 
## @param x A matrix of logicals were TRUE indicates the feature i is 
##     found in sample j. 
.computeMissingValueMetrics <- function(x) {
    c(
        LocalSensitivityMean = mean(colSums(x != 0)),
        LocalSensitivitySd = sd(colSums(x != 0)),
        TotalSensitivity = sum(rowSums(x) != 0),
        Completeness = mean(x != 0),
        NumberCells = ncol(x)
    )
}

## Internal function that computes the Jaccard index between each pair
## of columns in a given identification matrix. Note that the Jaccard
## matrix is a symmetric matrix, so we only return the upper triangle.
## WARNING: for large number of cells, eg n > 10,000, this may become
## inefficient. When the time comes this problem is relevant, 
## implement sub-sampling.
## @param x A matrix of logicals were TRUE indicates the feature i is 
##     found in sample j. 
.computeJaccardIndex <- function(x) {
    vectorSizes <- colSums(x)
    pwVectorSizes <- sapply(vectorSizes, function(xx) vectorSizes + xx)
    unionSize <- crossprod(x)
    jacc <- unionSize / (pwVectorSizes - unionSize)
    jacc[upper.tri(jacc)]
}

.getSteps <- function(maxn, nsteps) {
    steps <- seq(1, maxn, length.out = min(maxn, nsteps))
    round(steps)
}

.sampledSensitivity <- function(x, n, i) {
    sel <- sample(1:ncol(x), n)
    s <- sum(rowSums(x[, sel, drop = FALSE]) != 0)
    data.frame(i = i, SampleSize = n, Sensitivity = s)
}

## Internal function that computes the cumulative sensitivity curve
## given an identification matrix. The function will sample an 
## increasing number of cells (depending on nsteps). Each sampling is
## repeated niters time. Optionally, if 'batch' is provided, the 
## function will aggregate columns within each batch. This is usefull
## when dealing with mulitplexed data were detection depends on the 
## batch. 
## @param x A matrix of logicals were TRUE indicates the feature i is 
##     found in sample j. 
.computeCumulativeSensitivityCurve <- function(x,
                                               batch = NULL,
                                               niters = 10, 
                                               nsteps = 30) {
    out <- list()
    if (!is.null(batch)) {
        x <- sapply(unique(batch), function(i) {
            ifelse(rowSums(x[, batch == i, drop = FALSE]) == 0, FALSE, TRUE)
        })
    }
    for (n in .getSteps(ncol(x), nsteps)) {
        for (i in 1:niters) {
            out <- c(out, list(.sampledSensitivity(x, n, i)))
        }
    }
    do.call(rbind, out)
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
##' @param colvar A `character(1)` indicating the variable to take
##'     from `colData(object)` that gives the sample annotation.
##' 
##' @param samplePattern A `character(1)` pattern that matches the
##'     sample encoding in `colvar`.
##' 
##' @param carrierPattern A `character(1)` pattern that matches the
##'     carrier encoding in `colvar`. Only one match per assay is
##'     allowed, otherwise only the first match is taken
##'     
##' @param sampleFUN A `character(1)` or `function` that provides the 
##'     summarization function to use (eg mean, sum, media, max, ...).
##'     Only used when the pattern matches multiple samples. Default 
##'     is `mean`. Note for custom function, `na.rm = TRUE` is passed
##'     to `sampleFUN` to ignore missing values, make sure to provide 
##'     a function that accepts this argument.
##'     
##' @param carrierFUN A `character(1)` or `function` that provides the 
##'     summarization function to use (eg mean, sum, media, max, ...).
##'     Only used when the pattern matches multiple carriers. Default 
##'     is the same function as `sampleFUN`. Note for custom function,
##'     `na.rm = TRUE` is passed to `carrierFUN` to ignore missing 
##'     values, make sure to provide a function that accepts this 
##'     argument.
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
##'                    colvar = "SampleType",
##'                    carrierPattern = "Carrier",
##'                    samplePattern = "Blank|Macrophage|Monocyte",
##'                    sampleFUN = "mean",
##'                    rowDataName = "MeanSCR")
##' ## Check results
##' rowData(scp1)[[1]][, "MeanSCR"]
##' 
computeSCR <- function(object, 
                       i, 
                       colvar, 
                       samplePattern, 
                       sampleFUN = "mean",
                       carrierPattern,
                       carrierFUN = sampleFUN,
                       rowDataName = "SCR") {
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
        annot <- colData(object)[colnames(object[[ii]]), ][, colvar]
        ## Get the corresponding indices
        if (is.numeric(samplePattern)) {
            sampIdx <- samplePattern[samplePattern <= length(annot)]
        } else {
            sampIdx <- grep(samplePattern, annot)
        }
        if (!length(sampIdx)) stop("No match found with 'samplePattern = ", 
                                   samplePattern, "'.")
        carrIdx <- grep(carrierPattern, annot)
        if (!length(carrIdx)) stop("No match found with 'carrierPattern = ", 
                                   carrierPattern, "'.")
        ## Get sample data
        samp <- assay(object[[ii]])[, sampIdx, drop = FALSE]
        if (ncol(samp) > 1)
            samp <- apply(samp, 1, sampleFUN, na.rm = TRUE)
        ## Get carrier data
        carrier <- assay(object[[ii]])[, carrIdx, drop = FALSE]
        if (ncol(carrier) > 1)
            carrier <- apply(carrier, 1, carrierFUN, na.rm = TRUE)
        ## Compute ratio
        ratio <- unname(samp / carrier)
        ## Store ratio in rowData
        rowData(object@ExperimentList@listData[[ii]])[, rowDataName] <- 
            ratio
        ## more efficient than 
        ## rowData(object[[ii]])[, rowDataName] <-  ratio
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
##' rowData(scp1)[[1]][, c("dart_PEP", "qvalue_protein")]
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
    
    ## Get the rowData from all assays
    df <- rbindRowData(object, i)
    
    ## Check groupNy and PEP
    vars <- PEP
    if (!missing(groupBy)) vars <- c(vars, groupBy)
    if (any(!vars %in% colnames(df))) 
        stop("Variable(s) not found in the 'rowData'. Make sure 'PEP'",
             "and 'groupBy' are present in all selected assays.")
    
    ## Check PEP is a probability
    pepRange <- range(df[, PEP], na.rm = TRUE)
    if (max(pepRange) > 1 | min(pepRange < 0))
        stop("'", PEP, "' is not a probability in (0, 1)")
    
    ## Compute the q-values 
    if (missing(groupBy)) {
        df$qval <- .pep2qvalue(df[, PEP])
    } else { ## Apply grouping if supplied
        p <- split(df[[PEP]], df[[groupBy]])
        p <- vapply(p, min, FUN.VALUE = numeric(1), na.rm = TRUE)
        qval <- .pep2qvalue(p)
        df$qval <- unname(qval[df[[groupBy]]])
    }
    
    ## Insert the q-value inside every assay rowData
    pepID <- paste0(df$assay, df$rowname)
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
##'     sequentially applied to each feature (row) in each assay. 
##'     Available methods and additional 
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
medianCVperCell <- function(object, i, groupBy, nobs = 5, na.rm = TRUE,
                            colDataName = "MedianCV", norm = "none", ...) {
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
        cvs <- featureCV(object[[ii]], na.rm = na.rm, norm = norm, nobs = nobs,
                         group = as.factor(rowData(object[[ii]])[, groupBy]),
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

##---- Report missing values ---- 

##' Four metrics to report missing values
##'
##' The function computes four metrics to report missing values in
##' single-cell proteomics. 
##'
##' @param object An object of class [QFeatures].
##' @param i The index of the assay in `object`. The assay must 
##'     contain an identification matrix, that is a matrix where an
##'     entry is `TRUE` if the value is observed and `FALSE` is the
##'     value is missing (see examples). `i` may be numeric, character
##'     or logical, but it must select only one assay.
##' @param by A vector of length equal to the number of columns in 
##'     assay `i` that defines groups for which the metrics should be
##'     computed separately. If missing, the metrics are computed for
##'     the complete assay.
##'
##' @return A `data.frame` with groups as rows and 5 columns: 
##' 
##' - `LocalSensitivityMean`: the average number of features per cell.
##' - `LocalSensitivitySd`: the standard deviation of the local
##'   sensitivity.
##' - `TotalSensitivity`: the total number of features found in the 
##'   dataset. 
##' - `Completeness`: the proportion of values that are not missing in
##'   the data.
##' - `NumberCells`: the number of cells in the dataset.
##' 
##' @export
##'
##' @examples
##' 
##' data("scp1")
##' 
##' ## Define the identification matrix
##' peps <- scp1[["peptides"]]
##' assay(peps) <- !is.na(assay(peps))
##' scp1 <- addAssay(scp1, peps, "id")
##' 
##' ## Report metrics 
##' reportMissingValues(scp1, "id")
##' ## Report metrics by sample type
##' reportMissingValues(scp1, "id", scp1$SampleType)
##' 
##' data
##' 
reportMissingValues <- function(object, i, by = NULL) {
    x <- .getDetectionMatrix(object, i)
    if (is.null(by)) by <- rep("all", ncol(x))
    stopifnot(length(by) == ncol(x))
    out <- sapply(unique(by), function(byi) {
        .computeMissingValueMetrics(x[, by == byi])
    })
    data.frame(t(out))
}

##' Compute the pairwise Jaccard index
##'
##' The function computes the Jaccard index between all pairs of cells.
##' 
##' @param object An object of class [QFeatures].
##' @param i The index of the assay in `object`. The assay must 
##'     contain an identification matrix, that is a matrix where an
##'     entry is `TRUE` if the value is observed and `FALSE` is the
##'     value is missing (see examples).
##' @param by A vector of length equal to the number of columns in 
##'     assay `i` that defines groups for which the Jaccard index
##'     should be computed separately. If missing, the Jaccard indices
##'     are computed for all airs of cells in the dataset.
##'
##' @return  A `data.frame` with as many rows as pairs of cells
##'  and the following column(s): 
##'
##' - `jaccard`: the computed Jaccard index
##' - `by`: if `by` is not `NULL`, the group of the pair of cells
##'   for which the Jaccard index is computed. 
##'
##' @export
##'
##' @examples
##'
##' data("scp1")
##'
##' ## Define the identification matrix
##' peps <- scp1[["peptides"]]
##' assay(peps) <- ifelse(is.na(assay(peps)), FALSE, TRUE)
##' scp1 <- addAssay(scp1, peps, "id")
##' 
##' ## Compute Jaccard indices
##' jaccardIndex(scp1, "id")
##' ## Compute Jaccard indices by sample type
##' jaccardIndex(scp1, "id", scp1$SampleType)
##'
##'
jaccardIndex <- function(object, i, by = NULL) {
    x <- .getDetectionMatrix(object, i)
    if (is.null(by)) by <- rep("all", ncol(x))
    stopifnot(length(by) == ncol(x))
    out <- lapply(unique(by), function(byi) {
        ji <- .computeJaccardIndex(assay(x[, by == byi]))
        data.frame(jaccard = ji, by = byi)
    })
    do.call(rbind, out)
}

##' Cumulative sensitivity curve
##'
##' @description
##'
##' The cumulative sensitivity curve is used to evaluate if the sample
##' size is sufficient to accurately estimate the total sensitivity. 
##' If it is not the case, an asymptotic regression model may provide
##' a prediction of the total sensitivity if more samples would have
##' been acquired.
##'
##' @param object An object of class [QFeatures].
##' @param i The index of the assay in `object`. The assay must 
##'     contain an identification matrix, that is a matrix where an
##'     entry is `TRUE` if the value is observed and `FALSE` is the
##'     value is missing (see examples).
##' @param by A vector of length equal to the number of columns in 
##'     assay `i` that defines groups for a cumulative sensitivity
##'     curve will be computed separately. If missing, the sensitivity
##'     curve is computed for the completd dataset.
##' @param batch  A vector of length equal to the number of columns in 
##'     assay `i` that defines the cell batches. All cells in a batch
##'     will be aggregated to a single sample.
##' @param nsteps The number of equally spaced sample sizes to compute
##'     the sensitivity.
##' @param niters The number of iteration to compute 
##'
##' @return  A `data.frame` with groups as many rows as pairs of cells
##'  and the following column(s): 
##'
##' - `jaccard`: the computed Jaccard index
##' - `by`: if `by` is not `NULL`, the group of the pair of cells
##'   for which the Jaccard index is computed. 
##'
##' @details
##' 
##' As more samples are added to a dataset, the total number of 
##' distinct features increases. When sufficient number of samples are
##' acquired, all peptides that are identifiable by the technology and
##' increasing the sample size no longer increases the set of 
##' identified features. The cumulative sensitivity curve depicts the
##' relationship between sensitivity (number of distinct peptides in
##' the data) and the sample size. More precisely, the curve is built
##' by sampling cells in the data and count the number of distinct
##' features found across the sampled cells. The sampling is repeated
##' multiple times to account for the stochasticity of the approach. 
##' Datasets that have a sample size sufficiently large should have a
##' cumulative sensitivity curve with a plateau.
##'
##' The set of features present in a cell depends on the cell type.
##' Therefore, we suggest to build the cumulative sensitivity curve
##' for each cell type separately. This is possible when providing the
##' `by` argument.
##'
##' For multiplexed experiments, several cells are acquired in a run.
##' In that case, when a features is identified in a cell, it is 
##' frequently also identified in all other cells of that run, and
##' this will distort the cumulative sensitivity curve. Therefore, the
##' function allows to compute the cumulative sensitivity curve at the
##' batches level rather than at the cell level. This is possible when
##' providing the `batch` argument.
##' 
##' Once the cumulative sensitivity curve is computed, the returned 
##' data can be visualized to explore the relationship between the 
##' sensitivity and the sample size. If enough samples are acquired,
##' the curve should plateau at high numbers of samples. If it is not
##' the case, the total sensitivity can be predicted using an
##' asymptotic regression curve. To predict the total sensitivity, the
##' model is extrapolated to infinite sample size. Therefore, the 
##' accuracy of the extrapolation will highly depend on the available
##' data. The closer the curve is to the plateau, the more accurate 
##' the prediction.
##'
##' @export
##'
##' @examples
##'
##' ## Simulate data
##' ## 1000 features in 100 cells
##' library(SingleCellExperiment)
##' id <- matrix(FALSE, 1000, 1000)
##' id[sample(1:length(id), 5000)] <- TRUE
##' dimnames(id) <- list(
##'     paste0("feat", 1:1000),
##'     paste0("cell", 1:1000)
##' )
##' sce <- SingleCellExperiment(assays = List(id))
##' sim <- QFeatures(experiments = List(id = sce))
##' sim$batch <- rep(1:100, each = 10)
##' sim$SampleType <- rep(c("A", "B"), each = 500)
##' sim
##'
##' ## Compute the cumulative sensitivity curve, take batch and sample
##' ## type into account
##' csc <- cumulativeSensitivityCurve(
##'     sim, "id", by = sim$SampleType,
##'     batch = sim$batch
##' )
##' predCSC <- predictSensitivity(csc, nSample = 1:50)
##' 
##' library(ggplot2)
##' ggplot(csc) +
##'     aes(x = SampleSize, y = Sensitivity, colour = by) +
##'     geom_point() +
##'     geom_line(data = predCSC)
##'
##' ## Extrapolate the total sensitivity
##' predictSensitivity(csc, nSamples = Inf)
##' ## (real total sensitivity = 1000)
##'
##' @rdname cumulativeSensitivityCurve
cumulativeSensitivityCurve <- function(object, i,  by = NULL,
                                       batch = NULL, nsteps = 30, 
                                       niters = 10) {
    x <- .getDetectionMatrix(object, i)
    if (is.null(by)) by <- rep("all", ncol(x))
    stopifnot(length(by) == ncol(x))
    csc <- lapply(unique(by), function(byi) {
        out <- .computeCumulativeSensitivityCurve(
            x[, by == byi, drop = FALSE], 
            batch[by == byi], 
            niters, nsteps
        )
        out$by <- byi
        out
    })
    do.call(rbind, csc)
}

##' @param df The output from `cumulativeSensitivityCurve()`.
##' @param nSamples A `numeric()` of samples sizes. If `Inf`, the
##'    prediction provides the extrapolated total sensitivity.
##'
##' @importFrom stats nls predict
##' @export
##' @rdname cumulativeSensitivityCurve
predictSensitivity <- function(df, nSamples) {
    stopifnot(c("SampleSize", "Sensitivity") %in% colnames(df))
    new <- data.frame(SampleSize = nSamples)
    if (!"by" %in% colnames(df)) {
        fit <- .modelSensitivity(df)
        new$Sensitivity <- predict(fit, newdata = new)
        new
    } else {
        out <- lapply(unique(df$by), function(byi) {
            fit <- .modelSensitivity(df[df$by == byi, ])
            new$Sensitivity <- predict(fit, newdata = new)
            new$by <- byi
            new
        })
        do.call(rbind, out)
    }
}

.modelSensitivity <- function(df) {
    nls(
        formula = Sensitivity ~ SSasymp(SampleSize, Asym, R0, lrc),
        data = df,
        weights = df$SampleSize^2 ## data is heterorscedastic
    )
}