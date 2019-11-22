library(dplyr)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(MSnbase)
library(nipals)
library(pRoloc)
library(scpdata)
library(tidyr)
library(vsn)

#' Normalization of single-cell proteomics data
#'
#' The function \code{scp_normalise} (and identically \code{scp_normalize}) 
#' performs normalization of cells or peptides/proteins or both of quantified
#' expression data (\code{MSnSet} objects). Different normalization methods are 
#' implemented (see Details). 
#'
#' @param obj 
#' An MSnSet object
#' @param what 
#' A character indicating whether columns (\code{"col"}), rows (\code{"row"}) or
#' both (\code{"both"}) should be normalized
#' @param method 
#' The type of normalization (see Details)
#'
#' @details 
#' 
#' \code{method == 1}: each row (if \code{what == "row"}) or column (if 
#' \code{what == "col"}) is subtracted by summary of the corresponding values. 
#' The summary function is the mean for rows and the median for columns. 
#'
#' @return
#' 
#' An \code{MSnSet} object similar to the input object \code{obj}, but where the 
#' expression values have been replaced by normalized expression values. 
#' 
#' @export
#'
#' @examples
#' sc <- scpdata("specht2019", type = "peptide")
#' scpNorm <- scp_normalize(sc, what = "both")
#' 
scp_normalize <- scp_normalise <- function(obj, what = "col", method = 1, fun = "-"){
  # Check arguments
  if(!inherits(obj, "MSnSet")) stop("'obj' should be an MSnSet object." )
  what <- match.arg(what, choices = c("col", "row", "both"))
  
  ## Method 1: implemented by Specht et al. 2019
  if(method == 1){
    # Normalize rows
    if(what %in% c("row", "both")) 
      obj <- stat_normalize(obj, margin = 1, stats = mean, fun = fun)
    # Normalize columns
    if(what %in% c("col", "both")) 
      obj <- stat_normalize(obj, margin = 2, stats = median, fun = fun)
  } else {
    stop("Method '", method, "' not implemented.")
  }
  
  return(obj)
}


#' Aggregation of peptide expression data into protein expression data
#'
#' The \code{scp_AggregateByProtein} function takes an \code{MSnSet} object and 
#' merges the peptides (rows) of the expression data into proteins. This merging
#' is performed by taking the median value for each cell (column) for a given
#' protein group. The peptide to protein mapping should be present in the 
#' feature data of the \code{\link{MSnSet}} object.
#'
#' @param obj 
#' An \code{MSnSet} object containing the peptide expression data.
#'
#' @return
#' A new \code{MSnSet} object containing the aggregated protein expression data.
#' Note that the feature data is adapted and will contain only protein 
#' annotations.
#' 
#' @export
#'
#' @examples
#' data(specht2019)
#' scProt <- scp_aggregateByProtein(specht2019)
#' 
scp_aggregateByProtein <- function(obj){
  if(nrow(fData(obj)) == 0) stop("'fData(obj)' cannot be empty")
  return(aggregateByProtein(obj))
}



#' @export
scp_plotStats <- function(obj, what = c("both", "cells", "features"), 
                          xstat = "mean", ystat = "sd"){
  # Check and format arguments 
  what <- match.arg(what, c("both", "cells", "features"))
  if(is.character(xstat)){
    xstat <- switch(xstat, 
                    mean = list(fun = function(x) mean(x, na.rm = TRUE), title = "Mean"),
                    median = list(fun = function(x) round(median(x, na.rm = TRUE), 15), title = "Median"),
                    sd = list(fun = function(x) sd(x, na.rm = TRUE), title = "Standard deviation"),
                    var = list(fun = function(x) var(x, na.rm = TRUE), title = "Variance"))
  }
  if(is.character(ystat)){
    ystat <- switch(ystat, 
                    mean = list(fun = function(x) mean(x, na.rm = TRUE), title = "Mean"),
                    median = list(fun = function(x) round(median(x, na.rm = TRUE), 15), title = "Median"),
                    sd = list(fun = function(x) sd(x, na.rm = TRUE), title = "Standard deviation"),
                    var = list(fun = function(x) var(x, na.rm = TRUE), title = "Variance"))
  }
  M <- exprs(obj)
  
  # Generate plots
  if(what %in% c("both", "cells")){
    p1 <- ggplot(data = data.frame(var1 = apply(M, 2, xstat$fun),
                                   var2 = apply(M, 2, ystat$fun))) + 
      geom_point(aes(x = var1, y = var2), col = rgb(0, 0, 0.5, 0.5)) +
      ggtitle("Distribution statistics for single cells") +
      xlab(xstat$title) + ylab(ystat$title)
  }
  if(what %in% c("both", "features")){
    p2 <- ggplot(data = data.frame(var1 = apply(M, 1, xstat$fun),
                                   var2 = apply(M, 1, ystat$fun))) + 
      geom_point(aes(x = var1, y = var2), col = rgb(0, 0, 0.5, 0.5)) +
      xlab(xstat$title) + ylab(ystat$title) +
      ggtitle("Distribution statistics for features")
  }
  
  # Print plots
  if(what == "both"){
    grid.arrange(p1, p2, nrow = 1)
  } else if (what == "features"){
    print(p1)
  } else if (what == "cells"){
    print(p2)
  }
}

#' @export
scp_plotMissing <- function(obj, what = c("both", "samples", "features")){
  # Check and format arguments 
  what <- match.arg(what[1], c("both", "cells", "features")) 
  
  M <- exprs(obj)
  
  if(what %in% c("both", "samples")){
    mis <- apply(M, 2, function(x) sum(is.na(x))/length(x)) * 100
    p1 <- ggplot(data = data.frame(mis), mapping = aes(x = mis)) + 
      geom_histogram(binwidth = 2.5, fill = "grey60", col = "grey40") +
      geom_vline(aes(xintercept = mean(mis)), col = "red", linetype = "dashed") + 
      xlim(0, 100) + xlab("Missingness (%)") + ggtitle("Missing data among samples")
  }
  if(what %in% c("both", "features")){
    mis <- apply(M, 1, function(x) sum(is.na(x))/length(x)) * 100
    p2 <- ggplot(data = data.frame(mis), mapping = aes(x = mis)) + 
      geom_histogram(binwidth = 2.5, fill = "grey60", col = "grey40") +
      geom_vline(aes(xintercept = mean(mis)), col = "red", linetype = "dashed") + 
      xlim(0, 100) + xlab("Missingness (%)") + ggtitle("Missing data among features")
  }
  
  # Print plots
  if(what == "both"){
    grid.arrange(p1, p2, nrow = 1)
  } else if (what == "features"){
    print(p1)
  } else if (what == "cells"){
    print(p2)
  }
}

#' @export
scp_plotCV <- function(obj, groupBy = "Protein", colorBy = NULL){
  dat <- exprs(obj)
  prots <- as.character(fData(obj)[, "Proteins"])
  CVs <- do.call(rbind, lapply(unique(prots), function(prot){
    .idx <- prots == prot
    if(sum(.idx) < 2) return(NULL)
    xx <- dat[.idx,]
    CV <- apply(xx, 2, sd, na.rm = TRUE)/apply(xx, 2, mean, na.rm = TRUE)
    return(CV)
  }))
  dimnames(CVs) <- list(NULL, NULL)
  CVs <- reshape2::melt(CVs)
  CVs <- CVs[!is.na(CVs$value),]
  # CVs <- cbind(CVs, type = pData(obj)$celltype[CVs$Var2])
  ggplot(data = CVs, aes(x = Var2, y = value)) +
    geom_point(col = rgb(0.3,0.3,0.3,0.2)) + 
    stat_summary(aes(y = value, group = 1, color = "median"), fun.y = median,
                 geom = "point",  size = 0.5, group = 1) +
    scale_color_manual("", values = c("median" = "red")) +
    xlab("Cell index") + ylab("Coefficient of variation") +
    ggtitle("Protein CV distribution per cell")
}

#' @export
show_heatmap <- function(obj, ylab = "Feature index", xlab = "Cell index", 
                         Log2 = FALSE, znorm = FALSE, main, Legend = TRUE,
                         feat_ord = NULL, samp_ord = NULL, trim = NA, ...){
  dat <- t(exprs(obj))
  if(missing(main)){
    main <- paste0("Expression matrix (", ncol(dat), " features x ",
                   nrow(dat), " cells)")
  }
  if(Log2){
    dat <- log2(dat)
    dat[is.infinite(dat)] <- NA
  }
  if(znorm){ # z-normalize features
    dat <- apply(dat, 2, function(x){
      x <- x - mean(x, na.rm = TRUE)
      x <- x / sd(x, na.rm = TRUE)
      return(x)
    })
  }
  if(is.na(trim)) trim <- round(max(abs(dat), na.rm = TRUE), 2)
  dat[dat < -trim] <- -trim
  dat[dat > trim] <- trim
  
  if(!is.null(feat_ord)) dat <- dat[, feat_ord]
  if(!is.null(samp_ord)) dat <- dat[samp_ord, ]
  
  cols <- colorRampPalette(c(rgb(0.55,0.55,1), rgb(1, 0.95, 0.7), rgb(1, 0.45, 0)))(1000)
  init.par <- par()[c("mar", "mai")]
  layout(t(1:2), widths=c(8,1))
  par(mar = c(2, 2, 2, 0))
  image(x = 1:nrow(dat), y = 1:ncol(dat), z = dat, xlab = xlab, ylab = ylab,
        useRaster = TRUE, axes = FALSE, col = cols, mgp = c(1, 0, 0),
        zlim = c(-trim, trim),
        main = main)
  box()
  ds <- dev.size()
  par(mai = c(ds[2]/4, ds[1]/50, ds[2]/4, ds[1]/20))
  image(z = t(seq(-trim, trim, length.out = 1000)), col = cols, axes = FALSE,
        cex = 0.1, mgp = c(0, 1, 0.5), useRaster = TRUE)
  axis(4, at = c(0, 0.5, 1),  labels = c(-trim, 0, trim), cex = 0.1)
  box()
  par(init.par)
  return(invisible(NULL))
}

plotSCoPEset <- function(obj, run, phenotype = NULL){
  if(!all(c("channel", "run", phenotype) %in% colnames(pData(obj))))
    stop(paste0("'pData(obj)' must contain the following fields: ", 
                paste0(c("channel", "run", phenotype), collapse = ", ")))
  if(! run %in% pData(obj)$run) stop(paste0("The run '", run, "' is not present in the MSnSet."))
  # Subset for the desired run
  obj <- obj[, pData(obj)$run == run]
  # Format the data
  df <- data.frame(run = run, channel = pData(obj)$channel, t(exprs(obj)))
  df <- pivot_longer(data = df, cols = -(1:2), values_to = "intensity") 
  # Get counts per channel
  df %>%  group_by(channel) %>% 
    summarise(max = max(intensity, na.rm = TRUE), 
              n = sum(!is.na(intensity)),
              mean = mean(intensity, na.rm = TRUE),
              median = median(intensity, na.rm = TRUE)) -> counts
  labs <- if(!is.null(phenotype))
    array(pData(obj)[,phenotype], dimnames = list(pData(obj)$channel))
  else waiver()
  # Create the plot
  p <- ggplot(data = df, aes(x = channel, y = intensity)) +
    geom_violin(na.rm = TRUE) + scale_y_log10() + ggtitle(run) + 
    geom_point(data = counts, aes(x = channel, y = median, shape = "+"), 
               color = "red", size = 5) + 
    geom_text(data = counts, aes(x = channel, y = max*2, 
                                 label = paste0("n=", n)),
              size = 4, color = "grey50") + 
    scale_shape_manual(values = "+",
                       labels = c(`+` = "Median"),
                       name = "") + 
    scale_x_discrete(limits = as.character(0:10), 
                              labels = labs)
  return(p)
}

plotCorQC <- function(obj, run, na.rm = FALSE){
  X <- exprs(obj)[,pData(obj)$run == run]
  if(na.rm){
    X <- X[-unique(which(is.na(X), arr.ind = T)[,1]),]
    # X[is.na(X)] <- 0
  }
  ct <- pData(obj)$cell_type[pData(obj)$run == run]
  image(cor(X), axes = F, 
        main = paste0("Pearson correlation between channels\n", run))
  mtext(ct, 1, at = seq(0,1,length.out = length(ct)), las = 2)
  mtext(ct, 2, at = seq(0,1,length.out = length(ct)), las = 1)
}

plotSlavovQC <- function(obj)


customPCA <- function(obj, pca, x = "PC1", y = "PC2", color = "cell_type", shape = "batch"){
  PCs <- pca$loadings
  meta <- pData(obj)
  meta$batch <- paste0(meta$lcbatch, "-",
                       sapply(as.character(meta$run), function(x) tail(strsplit(x, "")[[1]], 2)[1]))
  p <- ggplot(data = data.frame(PCs, meta)) +
    geom_point(aes(x = eval(parse(text = x)), y = eval(parse(text = y)), 
                   color = eval(parse(text = color)), 
                   shape = eval(parse(text = shape)))) +
    scale_color_manual(name = color,
                       values = c(carrier_mix = "grey50", 
                                  unused = "grey70",
                                  norm = "darkseagreen",
                                  sc_0 = "bisque3", 
                                  sc_m0 = "cornflowerblue",
                                  sc_u = "coral")) + 
    scale_shape(name = shape) + 
    ggtitle("PCA on the expression data") + xlab(x) + ylab(y)
  return(p)
}



####---- SPECHT ET AL. 2019 FUNCTIONS ----####

# This part of the script contains the code/algorithms used by Specht et al 
# (2019) for processing their SCP data

stat_normalize <- function(obj, margin = 1, stats = mean, fun = "-"){
  stats <- apply(exprs(obj), MARGIN = margin, stats, na.rm = TRUE)
  exprs(obj)  <- sweep(exprs(obj), MARGIN = margin, STATS = stats, FUN = fun)
  return(obj)
}

aggregateByProtein <- function(obj){
  prots <- fData(obj)[,1]
  x <- do.call(rbind, lapply(unique(prots), function(prot){
    xx <- exprs(obj)[prots == prot, , drop = F]
    apply(xx, 2, median, na.rm = TRUE)
  }))
  rownames(x) <- unique(prots)
  obj.new <- MSnSet(exprs = x, 
                    fData = data.frame(protein = unique(prots), row.names =  unique(prots)),
                    pData = pData(obj),
                    experimentData = experimentData(obj))
  return(obj.new)
}


#' Impute missing values using K-nearest neighbours
#' 
#' Internal function 
#' @export
imputeKNN <- function(obj, k = 3){
  dat <- exprs(obj)
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
      }
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ])
      }
    }
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
  }
  exprs(obj) <- dat.imp
  return(obj)
}

#' @export
batchCorrect <- function(obj, batch, target){
  if(is.character(batch)){
    batch <- pData(obj)[, batch]
  } else if (!is.factor(batch)){
    stop("'batch' should be either a column name (character) in pData(obj) or a factor")
  }
  if(is.character(target)){
    target <- model.matrix(~ as.factor(pData(obj)[, target]))
  } else if (!is.matrix(target)){
    stop("'target' should be either a column name (character) in pData(obj) or a design matrix")
  }
  exprs(obj) <- ComBat(dat = exprs(obj), batch = batch, mod = target, par.prior = T)
  return(obj)
}
