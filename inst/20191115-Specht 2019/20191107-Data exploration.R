library(MSnbase)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)
library(sva)
library(pheatmap)
library(nipals)
library(pRoloc)


# laod the peptide data
library(scpdata)
sc <- specht2018_peptide

# Check the data 
show_heatmap(sc, znorm = TRUE, trim = 3)

# Check the distribution of peptide used for protein quantification 
pepPerProt <- table(fData(sc)$Proteins)
tab <- sapply(1:10, function(i){
  if(i == 10){
    return(sum(pepPerProt >= i))
  } else {
    return(sum(pepPerProt == i))
  }
})
names(tab) <- paste0(1:10, c(rep("", 10-1), "+"))
tab

# Check the missingness 
scp_plotMissing(sc)
# about **75%** of the data is missing ! However, the distributions are 
# drastically different depending if we look at peptide level or at cell level. 

# Check the influence of cell type on missingness 
mis <- apply(exprs(sc), 2, function(x) sum(is.na(x))/length(x)) * 100
df <- data.frame(mis = mis, celltype = pData(sc)$celltype)
ggplot(data = df, mapping = aes(x = mis, color = celltype, fill = celltype)) + 
  geom_histogram(binwidth = 1) +
  xlab("Missingness (%)") + ggtitle("Missing data: macrophages vs monocytes") +
  
  scale_color_discrete(guide = FALSE) +
  scale_fill_discrete(name = "Cell type",
                      labels = c("Macrophages", "Monocytes"))
# There is no obvious relation between the rate of missing data and the cell 
# type.
# Compare the amount of missingness bewteen the 2 cell types 
df <- do.call(cbind, lapply(unique(pData(sc)$celltype), function(x){
  .sub <- exprs(sc)[, pData(sc)$celltype == x]
  mis <- rowSums(is.na(.sub))/ncol(.sub)*100
  logFC <- apply(.sub, 1, median, na.rm = TRUE)
  out <- data.frame(mis, logFC)
  colnames(out) <- paste0(c("mis", "logFc"), "_", x)
  return(out)
}))
df$relFC <- df$logFc_sc_m0 - df$logFc_sc_u
df$relFC[df$relFC > 3] <- 3
df$relFC[df$relFC < -3] <- -3
ggplot(data = df, aes(x = mis_sc_m0, y = mis_sc_u, col = relFC), alpha = 0.2) +
  geom_point() + 
  scale_color_gradient2(name = "log2(FC_M0/FC_mono)", low = "darkgreen", 
                        high = "red3", midpoint = 0, mid = "wheat") + 
  ylab("Missingness (%) in monocytes") + xlab("Missingness (%) in macrophages")

# Check relation between missingness and peptide intensity
# TODO 
# we don't have the peptide intensities yet
# df <- data.frame(mis = apply(exprs(sc), 1, function(x) sum(is.na(x))/length(x)) * 100,
#                  int = apply(exprs(sc), 1, median, na.rm = TRUE))
# ggplot(data = df, aes(x = int, y = mis)) + 
#   geom_point(color = "coral", size = 0.75, alpha = 0.35) +
#   xlab("Average intensity (arbitrary units)") + ylab("Missingness (%)") +
#   ggtitle("Relationship between missingness and peptide intensity")


# Basic statistics
scp_plotStats(sc, xstat = "mean", ystat = "var")

# Coefficient of variation
scp_plotCV(sc)

# Plot PCA
pca <- nipals(exprs(sc), ncomp = 5, center = TRUE, scale = TRUE)
PCs <- pca$loadings
meta <- data.frame(celltype = pData(sc)$celltype)
meta$batch <- paste0(pData(sc)$batch_chromatography, "-",
                     sapply(as.character(pData(sc)$raw.file), function(x) tail(strsplit(x, "")[[1]], 2)[1]))
ggplot(data = data.frame(PCs, meta)) +
  geom_point(aes(x = PC2, y = PC3, color = celltype, shape = batch)) +
  ggtitle("PCA on the expression data ") + 
  scale_color_discrete(name = "Cell type",
                       labels = c("Macrophages", "Monocytes"))



