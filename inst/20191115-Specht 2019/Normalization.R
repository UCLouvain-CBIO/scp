
####---- Description ----####

# Different normalization methods are tested. The objective is to retrieve a 
# good separation of the macrophages and the monocytes as seen in the paper and 
# in previous exploration (see scp/vignette/labmeeting/). 


####---- Setup environment ----####

# Load libraries, functions and data
source("./R/utils.R")
dat <- specht2019_peptide2
# Subset the data of interest
run_all <- unique(pData(dat)$run)
run_sc <- run_all[! grepl("_QC_|col19|col2[0-4]", run_all)]
dat_sc <- dat[, pData(dat)$run %in% run_sc] # Keep only single cell runs
dat_sc <- dat_sc[, pData(dat_sc)$digest %in% c("Q", "N")] # keep only good digestion batches

na.thrs <- 0.99 # same threshold as in Specht et al. 

####---- Initial PCA ----#
dat_sub <- dat_sc[, pData(dat_sc)$cell_type %in% c("sc_m0", "sc_u")] %>%
  filterNA(pNA = na.thrs)
ed <- exprs(dat_sub)
ed[ed < 10] <- NA
exprs(dat_sub) <- ed
dat_sub <- filterNA(dat_sub, pNA = na.thrs)
pca_init <- nipals(exprs(dat_sub), ncomp = 5, center = TRUE, scale = TRUE)
png("./inst/20191115-Specht 2019/figs/PCA - peptides processed.png", res = 300, height = 2000, width = 2000)
customPCA(dat_sub, pca_init, x = "PC1", y = "PC5")
dev.off()


# Filtering
# Filtering peptides that are more abundant than expected, ie >10% of carrier 
# channel 
ed <- exprs(dat_sc)
for(run in pData(dat_sc)$run){
  .idx_run <- pData(dat_sc)$run == run
  .idx_carrier <- .idx_run & pData(dat_sc)$cell_type == "carrier_mix"
  .idx_sc <- .idx_run & grepl("sc_", pData(dat_sc)$cell_type)
  ratio <- rowMeans(ed[, .idx_sc]/ed[, .idx_carrier], na.rm = TRUE)
  ed[ratio > 0.1, .idx_sc] <- NA
}
exprs(dat_sc) <- ed
stopifnot(validObject(dat_sc))
# Filtering NA's 
dat_sc <- filterNA(dat_sc, pNA = na.thrs)


####---- Initial PCA plot ----####

# In this PCA plot the macrophages and the monocytes are not separated 
pca <- nipals(exprs(dat_sc), ncomp = 2, center = TRUE, scale = TRUE)
customPCA(dat_sc, pca, x = "PC1", y = "PC2")

####---- Specht et al. method (2019) ----####

# Step 1: normalize every run with the internal reference channel
dat_scn <- dat_sc
ed <- exprs(dat_scn)
for(run in pData(dat_scn)$run){
  .idx_run <- pData(dat_scn)$run == run
  .idx_norm <- .idx_run & pData(dat_scn)$cell_type == "norm"
  ed[, .idx_run] <- ed[, .idx_run]/ed[, .idx_norm]
}
exprs(dat_scn) <- ed
stopifnot(validObject(dat_scn))
# Remove the normalization channels since they all are 1's or NA's
dat_scn <- dat_scn[,pData(dat_scn)$cell_type != "norm"]
dat_scn <- filterNA(dat_scn, pNA = na.thrs)
# Update the PCA plot
pca1 <- nipals(exprs(dat_scn), ncomp = 2, center = TRUE, scale = TRUE)
customPCA(dat_scn, pca1, x = "PC1", y = "PC2")
# PCA results: PCA does not converge... Strange ?

# Step 2: divide rows and columns
dat_scn <- scp_normalise(dat_scn, what = "both", method = 1, fun = "/")
# Update the PCA plot
pca2 <- nipals(exprs(dat_scn), ncomp = 5, center = TRUE, scale = TRUE)
customPCA(dat_scn, pca2, x = "PC1", y = "PC2")
# PCA results: the macrophages and monocytes are still not separated 

# Step 3: log2 transform
dat_scn <- log(dat_scn, base = 2)
# Update the PCA plot
pca3 <- nipals(exprs(dat_scn), ncomp = 5, center = TRUE, scale = TRUE)
customPCA(dat_scn, pca3, x = "PC1", y = "PC2")
# PCA results: the macrophages and monocytes are still not separated, and the 


test <- specht2019_peptide
pcatest <- nipals(exprs(test), ncomp = 2, center = TRUE, scale = TRUE)
customPCA(test, pcatest, x = "PC1", y = "PC2", color = "celltype", shape = "batch_chromatography")

ed_specht <- read.csv("../scpdata/inst/extdata/specht2019/peptides-raw-RI.csv", row.names = 1)
ed_specht <- ed_specht[rowSums(is.na(ed_specht)) != ncol(ed_specht), ]
ed_specht <- ed_specht[, colSums(is.na(ed_specht)) != nrow(ed_specht)]
pcatest <- nipals(ed_specht, ncomp = 5, center = TRUE, scale = TRUE)
PCs <- pcatest$loadings
meta <- t(read.csv("../scpdata/inst/extdata/specht2019/Cells.csv", row.names = 1))
meta <- meta[colnames(ed_specht),]
meta <- as.data.frame(meta)
meta$batch <- paste0(meta$batch_chromatography, "-",
                     sapply(as.character(meta$raw.file), function(x) tail(strsplit(x, "")[[1]], 2)[1]))
png("./inst/20191115-Specht 2019/figs/PCA - peptides-raw-RI csv.png", res = 300, height = 2000, width = 2000)
ggplot(data = data.frame(PCs, meta)) +
  geom_point(aes(x = PC4, y = PC5,
                 color = celltype, 
                 shape = batch)) +
  scale_color_manual(values = c(carrier_mix = "grey50", 
                                unused = "grey70",
                                norm = "darkseagreen",
                                sc_0 = "bisque3", 
                                sc_m0 = "cornflowerblue",
                                sc_u = "coral")) + 
  ggtitle("PCA on the expression data ('peptides-raw-RI.csv')") 
dev.off()






