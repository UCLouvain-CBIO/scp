
####---- Description ----####

# This script performs some basic QC on the peptide data from the Specht et al. 
# 2019 paper. 

####---- Setup environment ----####

# Load libraries, functions and data

source("./R/utils.R")
dat <- specht2019_peptide2

# Split data

run_all <- unique(pData(dat)$run)

# QC samples
run_qc <- grep("_QC_", pData(dat)$run, value = TRUE)
dat_qc <- dat[, pData(dat)$run %in% run_qc] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's

# Ladder samples
run_lad <- grep("col24", run_all, value = TRUE)
dat_lad <- dat[, pData(dat)$run %in% run_lad] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's

# 1000 cells samples 
run_1000 <- grep("col23", run_all, value = TRUE)
dat_1000 <- dat[, pData(dat)$run %in% run_1000] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's

# 100 cells samples 
run_100 <- grep("col2[12]", run_all, value = TRUE)
dat_100 <- dat[, pData(dat)$run %in% run_100] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's

# 10 cells samples 
run_10 <- grep("col19|col20", run_all, value = TRUE)
dat_10 <- dat[, pData(dat)$run %in% run_10] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's

# Single cell samples
run_sc <- run_all[!run_all %in% c(run_qc, run_lad, run_1000, run_100, run_10)]
dat_sc <- dat[, pData(dat)$run %in% run_sc] %>%
  filterNA(pNA = 1-1E-6) # Remove rows containing only NA's
# Note specht et al only keep exp FP94, FP97 for sc runs


####---- Check the 1000 cells samples ----####

# Intensity per channel 
pdf("./inst/20191115-Specht 2019/figs/SCoPE2 sets - 1000 cells.pdf")
for(run in unique(pData(dat_1000)$run)){
  p <- plotSCoPEset(obj = dat_1000, 
                    run = run,
                    phenotype = "cell_type")
  print(p)
}
dev.off()


####---- Check the 100 cells samples ----####

pdf("./inst/20191115-Specht 2019/figs/SCoPE2 sets - 100 cells.pdf")
for(run in unique(pData(dat_100)$run)){
  p <- plotSCoPEset(obj = dat_100, 
                    run = run,
                    phenotype = "cell_type")
  print(p)
}
dev.off()


####---- Check the 10 cells samples ----####

pdf("./inst/20191115-Specht 2019/figs/SCoPE2 sets - 10 cells.pdf")
for(run in unique(pData(dat_10)$run)){
  p <- plotSCoPEset(obj = dat_10, 
                    run = run,
                    phenotype = "cell_type")
  print(p)
}
dev.off()


####---- Check the 1 cells samples ----####

# Plots row standard deviations versus row means
meanSdPlot(dat_sc, ranks = TRUE) # from MSnbase

show_heatmap(dat_sc, Log2 = TRUE, znorm = TRUE, trim = 5)
plotNA(dat_sc) # from MSnbase
scp_plotMissing(dat_sc) 

pdf("./inst/20191115-Specht 2019/figs/SCoPE2 sets - single cells.pdf")
for(run in unique(pData(dat_sc)$run)){
  p <- plotSCoPEset(obj = dat_sc, 
                    run = run,
                    phenotype = "cell_type")
  print(p)
}
dev.off()


## Explore the missingness in single cells

mis <- apply(exprs(dat_sc), 2, function(x) sum(is.na(x))/length(x)) * 100
df <- data.frame(mis = mis, pData(dat_sc)) 
# Missingness with respect to sample type
ggplot(data = df, mapping = aes(x = mis, fill = cell_type)) + 
  geom_histogram(binwidth = 1) +
  scale_fill_discrete(name = "Sample type") +
  xlab("Missingness (%)") + ggtitle("Effect of sample type on missingness")
# Missingness with respect to chromatographic batch
ggplot(data = df, mapping = aes(x = mis, fill = lcbatch)) + 
  geom_histogram(binwidth = 1) +
  scale_fill_discrete(name = "Chromatographic batch") +
  xlab("Missingness (%)") + ggtitle("Effect of chromatographic batch on missingness")
# Missingness with respect to channel
ggplot(data = df, mapping = aes(x = mis, fill = channel)) + 
  geom_histogram(binwidth = 1) +
  scale_fill_discrete(name = "TMT Channel index") +
  xlab("Missingness (%)") + ggtitle("Effect of TMT channel on missingness")
# Missingness with respect to digest
ggplot(data = df, mapping = aes(x = mis, fill = digest)) + 
  geom_histogram(binwidth = 1) +
  scale_fill_discrete(name = "Digestion batch") +
  xlab("Missingness (%)") + ggtitle("Effect of digestion on missingness")
# CCL: the digestion batch/technique has a big impact on the amount of missing
# data. Hence we only keep the Q and N batches which contain lower missingness
dat_sc <- dat_sc[, pData(dat_sc)$digest %in% c("Q", "N")]
dat_sc %<>% filterNA(pNA = 1-1E-6)

# Check the missingness in macrophages versus monocytes
df <- do.call(cbind, lapply(c("sc_m0",  "sc_u"), function(x){
  .sub <- exprs(dat_sc)[, pData(dat_sc)$cell_type == x]
  mis <- rowSums(is.na(.sub))/ncol(.sub)*100
  median_intensity <- apply(.sub, 1, median, na.rm = TRUE)
  out <- data.frame(mis, median_intensity)
  colnames(out) <- paste0(c("mis", "median_intensity"), "_", x)
  return(out)
}))
df$logFC <- log2(df$median_intensity_sc_m0/df$median_intensity_sc_u)
df$logFC[df$logFC > 2.5] <- 2.5
df$logFC[df$logFC < -2.5] <- -2.5
ggplot(data = df, aes(x = mis_sc_m0, y = mis_sc_u, col = logFC), alpha = 0.2) +
  geom_point() + geom_abline() +
  scale_color_gradient2(name = "log2(macro/mono)", low = "darkgreen", 
                        high = "red3", midpoint = 0, mid = "wheat") + 
  ylab("Missingness (%) in monocytes") + xlab("Missingness (%) in macrophages")
# Conclusion: peptides that are more missing in monocytes than in macrophage are 
# less expressed in monocytes compared to macrophages

## Explore the peptide missingness

# Only use digestion batch Q and N as the others are biased towards higher 
# missingness (see above)
df <- data.frame(peptide = featureNames(dat_sc),
                 mis = apply(exprs(dat_sc), 1, function(x) sum(is.na(x))/length(x) * 100),
                 n_run = apply(exprs(dat_sc), 1, function(x) length(unique(pData(dat_sc)$run[!is.na(x)]))),
                 int = apply(exprs(dat_sc), 1, mean, na.rm = TRUE))
png("./inst/20191115-Specht 2019/figs/Missingness vs peptide intensity.png",
    res = 300, height = 2000, width = 2000)
p <- ggplot(data = df, aes(y = int, x = mis, color = n_run)) +
  geom_point(size = 0.75, alpha = 0.35, na.rm = TRUE) +
  stat_smooth(method = "glm", formula = y ~ poly(x, 2), se = FALSE, 
              na.rm = TRUE, col = "red3") +
  ylab("Average intensity (arbitrary units)") + xlab("Missingness (%)") +
  ylim(0, 50000) +
  scale_color_gradient2(name = "Found in\n# runs", low = "red3", mid = "bisque2", 
                        high = "darkgreen", midpoint = 10, 
                        breaks = seq(0, max(df$n_run), by = 10)) + 
  ggtitle("Relationship between missingness and peptide intensity")
p
dev.off()
# Strange, the highest peptides show the highest missingness. Maybe due to the 
# fact that with less observations, the mean is more subject to the effect of 
# outliers.
outl.idx <- which.max(df$int * df$mis^3)
p + geom_point(data = df[outl.idx,], aes(y = int, x = mis), color = "red") + 
  ylim(0, 110000)
which(!is.na(exprs(dat_sc)[outl.idx[1],]))
# These are all coming from the same set ! Hence it is probably better to remove
# those highly missing peptides. Specth et al discard peptides and cells with 
# missingness > 99 % .


## Check PCA

pca <- nipals(exprs(dat_sc), ncomp = 10, center = TRUE, scale = TRUE)
PCs <- pca$loadings
meta <- pData(dat_sc)
meta$batch <- paste0(meta$lcbatch, "-",
                     sapply(as.character(meta$run), function(x) tail(strsplit(x, "")[[1]], 2)[1]))
ggplot(data = data.frame(PCs, meta)) +
  geom_point(aes(x = PC2, y = PC1, color = cell_type, shape = batch)) +
  scale_color_manual(values = c(carrier_mix = "grey50", 
                                  unused = "grey70",
                                  norm = "darkseagreen",
                                  sc_0 = "bisque3", 
                                  sc_m0 = "cornflowerblue",
                                  sc_u = "coral")) + 
  ggtitle("PCA on the expression data")
# Conclusion: the PCA plot cannot distinguish between monocytes and macrophages.
# Furthermore, Carrier and normalization mix are contained in two separate 
# clusters at opposite sides of the single cell clusters, whereas the opposite 
# is expected. Performing PCA on the macrophages and monocytes only does not 
# help the separation. 
dat_sc_sub <- dat_sc[,pData(dat_sc)$cell_type %in% c("sc_m0", "sc_u")]
dat_sc_sub <- filterNA(dat_sc_sub, pNA = 0.25)
pca <- nipals(exprs(dat_sc_sub), ncomp = 5, center = TRUE, scale = TRUE)
PCs <- pca$loadings
meta <- pData(dat_sc_sub)
meta$batch <- paste0(meta$lcbatch, "-",
                     sapply(as.character(meta$run), function(x) tail(strsplit(x, "")[[1]], 2)[1]))
ggplot(data = data.frame(PCs, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = cell_type, shape = batch))


####---- Conclusion ----####


# * The specht2019 peptide data set contains different experiments and should be 
#   split accordingly. The lack of experimental design information for 
#   experiments other than single cell makes QC difficult
# * There is a high rate of missing data in the single cell data. Missingness is 
#   very high in the "O" and "P" digestion batches. The associated data is 
#   removed
# * Peptides that are more missing in monocytes than in macrophage are less 
#   expressed in monocytes compared to macrophages
# * There seem to be a inverse relationship betwen peptide intensity and peptide
#   missingness. However, the highest peptides show the highest missingness. 
#   This is because peptides are found at a very high intensities in only a few 
#   runs.
# * The PCA cannot separate macrophages and monocytes, and carrier and 
#   normalization samples cluster appart while they are made of the sames cells. 
#   Subsetting the data for macrophages and monocytes only doesn't help as well 
#   as filtering out highly missing peptides. 