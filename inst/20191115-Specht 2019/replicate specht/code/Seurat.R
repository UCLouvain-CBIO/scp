# Running seurat tutorial to obtain cluster of monocytes and use their raw copy numbers for computing ion counts
# All code taken from tutorial: https://satijalab.org/seurat/immune_alignment.html


#install.packages("Seurat")
library(Seurat)

# Load raw data: 
ctrl.data <- read.table("C:/Users/HMS-desk21/Desktop/immune_alignment_expression_matrices/immune_control_expression_matrix.txt.gz", 
                        sep = "\t")
stim.data <- read.table("C:/Users/HMS-desk21/Desktop/immune_alignment_expression_matrices/immune_stimulated_expression_matrix.txt.gz", 
                        sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)
# Set up stimulated object
stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim@meta.data$stim <- "STIM"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)

# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))

immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)

immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
                                 dims.align = 1:20)

# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)


new.ident <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", "B", "CD16 Mono", 
               "T activated", "CD8 T", "NK", "DC", "B activated", "Mk", "pDC", "Eryth")
for (i in 0:12) {
  immune.combined <- RenameIdent(object = immune.combined, old.ident.name = i, 
                                 new.ident.name = new.ident[i + 1])
}

# Obtain monocyte cluster subset
cd14.mono <- SubsetData(immune.combined, ident.use = "CD14 Mono", subset.raw = T)

# Save data
saveRDS(cd14.mono, "mono.RData")
mono.rna<-as.matrix(cd14.mono@raw.data)

write.csv(mono.rna, "mono.csv")
