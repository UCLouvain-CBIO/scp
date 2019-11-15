
## Normalization 

The normalization procedure here follows the normalization procedure from the original paper. The first step is to subtract from each row (peptides) the mean intensity of that row. This is performed using our `scp_normalize()` function.

```{r peptide_normalization, echo = TRUE, include=TRUE}
scNorm <- scp_normalize(sc, what = "row")
```

Next, the rows are collapsed by proteins, meaning that measurements for peptides that belong to the same proteins are merged in a single rozw. The rows from different peptides are summarized using the median value of the peptides.

```{r peptide_aggregation}
scNorm <- scp_aggregateByProtein(scNorm)
```

Rows than columns are finally normalized again. For the rows the **mean** value is subtracted; for the columns the **median** value is subtracted. Note that the argument `"both"` implies that rows then colmuns are normalized.

```{r protein_cell_normalization}
scNorm <- scp_normalize(scNorm, what = "both")
```

This can be performed in a single run using pipes:
  ```{r normalization_pipe, eval = FALSE}
sc %>% scp_normalize(what = "row") %>%
  scp_aggregateByProtein() %>%
  scp_normalize(what = "both") -> scNorm
```

In the above section, we looked at the mean and standard deviation distributions of the original data set. We can redraw those plots after data normalization.

```{r, echo = FALSE}
scp_plotStats(scNorm, xstat = "median", ystat = "sd")
```

Since the last step is subracting the median intensity from every cell, it is expected that the median intensity is constant. For individual proteins, the median intensities seem to be normally distributed around 0. It should be noted that the normalization procedure used by Specht and colleagues does not normalize the variance. This should be taken into account when performing downstream analysis such as PCA. 

Other imputation methods could be tested. For instance, the `Linnorm` method is a popular method that has shown good performance on scRNA-Seq data.

```{r combat_bef_imput}
scBc <- batchCorrect(scNorm, batch = "raw.file", target = "celltype")

``` 



## Imputation 

we observed in the data exploration section that the `specht2019` data set is highly sparse and contains a majority of missing entries. This large proportion of missing data hinders downstream analysis and imputation is performed to fill in the gaps. In their paper, Specht et al. implemented the k-nearest neighbour imputation:
  
  ```{r imputation}
scImput <- imputeKNN(scNorm)
cat("Missing data:", 
    round(sum(is.na(exprs(scImput))) / length(exprs(scImput)) * 100, 2), 
    "%")
```

To asses the performance of this imputation, we have a look again at the data using heatmaps. Before imputation (after normalization), the data is very sparse.

```{r before_imput, echo = FALSE}
show_heatmap(scNorm, znorm = TRUE, trim = 3)
```

After imputation:
  
  ```{r after_imput, echo = FALSE}
show_heatmap(scImput, znorm = TRUE, trim = 3)
```

Since the distance measure used for the KNN is Euclidean distance, it is interesting to look at how the similarity between cells is affected by the imputation:
  
  ```{r correl_imput, echo = FALSE}
layout(t(1:3), widths=c(4,4,1))
par(mar =  c(6,2,2,0))

R <- as.matrix( dist(t(exprs(scNorm))))
image(R, xlab = "Cell index 1", ylab = "Cell index 2", main = "Before imputation", 
      useRaster = TRUE, axes = FALSE, mgp = rep(0.5, 3))

par(mar =  c(6,2,2,2))

R <- as.matrix( dist(t(exprs(scImput))))
image(R, xlab = "Cell index 1", ylab = "Cell index 2", main = "After imputation", 
      useRaster = TRUE, axes = FALSE, mgp = rep(0.5, 3))

par(mar = c(6,0,2,0.2))

bdf <- pData(scNorm)[,-3]
bdf$raw.file <- sapply(as.character(bdf$raw.file), function(x) tail(strsplit(x, "")[[1]], 2)[1])
bdf <- apply(bdf, 2, function(x) as.numeric(as.factor(x)))
image(t(bdf), axes = F, col = colorRampPalette(c("wheat", "darkgreen"))(5), 
      useRaster = TRUE, main = "Batch")
axis(1, at = seq(0, 1, length.out = 4), tick = F, las = 2,
     labels = c("Cell type", "Digestion", "Chromatogr.", "\"raw.file\""))
```

We can see the imputation increases the correlation that exists between 3 groups of cells. The groups are exactly correlated with the batches defined by the combination of the **chromatographic batch** and another **undocumented batch** (the batch could be retrieved from the file names...). This strongly suggests that the SCoPE2 pipeline is subject to batch effect that cannot be corrected by the normalization to the reference channel. Therefore, batch correction using `ComBat` is applied (see below). 

Do the imputated data allow to distinguish between macrophages and monocytes ? Let's find this out with a PCA plot 

```{r pca_imput, echo = FALSE}
M <- t(exprs(scImput))
# Center data and scale to unit variance (PCA on correlation matrix)
M <- scale(M, center = TRUE, scale = TRUE)
# Compute weights: weight i is the fraction missing data in all peptides
# belonging to thz protein i across all filtered single cells. 
wi <- apply(exprs(scNorm), 1, function(x)  sum(is.na(x))/length(x))
wi <- wi/sum(wi)
# wi <- rep(1, ncol(M))
# Perform the PCA using SVD
scSVD <- svd(M %*% diag(sqrt(wi)))
eignval <- scSVD$d^2 / (length(scSVD$d) - 1) # convert singular values to eigen values
varexp <- round(eignval / sum(eignval) * 100) # percentage variance explained
PCs <- sweep(scSVD$u[,1:2], 2, scSVD$d[1:2], "*") # aka scores
# Plot PCA
ggplot(data = data.frame(PC = PCs, pData(scBc))) +
  geom_point(aes(x = PC.1, y = PC.2, color = celltype, shape = batch_chromatography)) +
  # geom_segment(data = data.frame(PC = scSVD$v[,1:2]),
  #              aes(x = 0, y = 0, xend = PC.1, yend = PC.2))
  ggtitle("PCA on the expression data ") + 
  scale_color_discrete(name = "Cell type",
                      labels = c("Macrophages", "Monocytes")) + 
  xlab(paste0("PC1 (", varexp[1], " %)")) + 
  ylab(paste0("PC2 (", varexp[2], " %)"))

```

Imputation using KNN has been shown to create artifacts when applied on scRNA-Seq data (see @Tian2019). Better methods could be `Drimpute` or `SAVER`. 

## Batch correction and data integration

Specht et al. performed batch correction of the single cell data using the `ComBat` function from the `sva` package. `ComBat` is a batch correction method using empirical Bayes (EB) procedure developed by Johnson et al. (2007) to remove variation linked to experimental variables (**eg** day of processing) while preserving the other sources of variation, namely the biological variation. The method is based on a previous model called the location and scale adjustment (L/S) model.

```{r batch_correction}
scBc <- batchCorrect(scImput, batch = "raw.file", target = "celltype")
```

Let's look again at the correlation heatmaps after reordering columns by cell type.

```{r correl_batch_correction, echo = FALSE}

.ord <- order(pData(scBc)$celltype)

layout(t(1:2), widths=c(5,1))
par(mar =  c(6,2,2,2))

R <- as.matrix( dist(t(exprs(scBc)[,.ord])))
image(R, xlab = "Cell index 1", ylab = "Cell index 2", main = "Batch corrected data", 
      useRaster = TRUE, axes = FALSE, mgp = rep(0.5, 3))

par(mar = c(6,0,2,0.2))

bdf <- pData(scBc)[.ord,-3]
bdf$raw.file <- sapply(as.character(bdf$raw.file), function(x) tail(strsplit(x, "")[[1]], 2)[1])
bdf <- apply(bdf, 2, function(x) as.numeric(as.factor(x)))
image(t(bdf), axes = F, col = colorRampPalette(c("wheat", "darkgreen"))(5), 
      useRaster = TRUE, main = "Batch")
axis(1, at = seq(0, 1, length.out = 4), tick = F, las = 2,
     labels = c("Cell type", "Digestion", "Chromatogr.", "\"raw.file\""))


```

It seems that the batch effects seen in the previous plots are completely removed. `ComBat` is able to correct for batch effects without affecting the biological variation induced by cell type. Note that 62 batches (that is the 62 SCoPE2 sets) were used for removing the batch effects.

Let's see if correct clusters can be found. 
```{r pheatmap, echo = FALSE}
dat <- exprs(scBc)
bdf <- pData(scBc)[,-(2:3)]
bdf$raw.file <- sapply(as.character(bdf$raw.file), function(x) tail(strsplit(x, "")[[1]], 2)[1])
anncolors <- list(celltype = c(sc_m0 = "#FBB4AE", sc_u = "#B3CDE3"),
                  batch_chromatography = c(LCA9 = "#FBB4AE", LCA10 = "#B3CDE3"),
                  raw.file = c(A = "#FBB4AE", B = "#B3CDE3"))
pheatmap(dat,  show_colnames = FALSE, show_rownames = FALSE, scale = "none",
         annotation_col = bdf, clustering_method = "complete",
         annotation_legend = FALSE, annotation_colors = anncolors)
```

The clustering (hierarchical clustering with complete linkage and Euclidean distance) almost perfectly groups the cell types together. Some batch effect whithin the cell type groups seems still present as batches are more grouped than expected by chance. This could be either due to remaining technical bias or an artifact of the KNN imputation. However, we can see that the main expression patterns are correlated with cell type and indicate that the data contain sufficient protein expression information to discriminate between the undifferentiated monocytes and the macrophages. 


## Conclusion

The data from Specht et al. published in June 2019 is a successful test case for MS-based SCP showing that such a technology can be used for biological purposes. The authors were able with their data to differentiate between macrophages and monocytes. Although the data exhibit high missingness, a clear-cut differences between the cell type could still be observed, and is even more highlighted by PCA. 


```{r svd, echo = FALSE}
M <- t(exprs(scBc))
# Center data and scale to unit variance (PCA on correlation matrix)
M <- scale(M, center = TRUE, scale = TRUE)
# Compute weights: weight i is the fraction missing data in all peptides
# belonging to thz protein i across all filtered single cells. 
wi <- apply(exprs(scNorm), 1, function(x)  sum(is.na(x))/length(x))
wi <- wi/sum(wi)
# wi <- rep(1, ncol(M))
# Perform the PCA using SVD
scSVD <- svd(M %*% diag(sqrt(wi)))
eignval <- scSVD$d^2 / (length(scSVD$d) - 1) # convert singular vectors to eigenvalues
varexp <- round(eignval / sum(eignval) * 100) # percentage variance explained
PCs <- sweep(scSVD$u[,1:2], 2, scSVD$d[1:2], "*") # aka scores
# Plot PCA
ggplot(data = data.frame(PC = PCs, pData(scBc))) +
  geom_point(aes(x = PC.1, y = PC.2, color = celltype)) +
  # geom_segment(data = data.frame(PC = scSVD$v[,1:2]),
  #              aes(x = 0, y = 0, xend = PC.1, yend = PC.2))
  ggtitle("PCA on the expression data ") + 
  scale_color_discrete(name = "Cell type",
                      labels = c("Macrophages", "Monocytes")) + 
  xlab(paste0("PC1 (", varexp[1], " %)")) + 
  ylab(paste0("PC2 (", varexp[2], " %)"))

```


This allows the identification of protein preferentially expressed in one cell type rather than the other. A list of proteins associated to the patterns (*eg* associated with the first PC) can easily be extracted and biological inference and interpretation can be performed. 

```{r list_prots, echo = FALSE}
pc1 <- data.frame(protein = fData(scBc)[,1], svec = scSVD$v[, 1])
pc1 <- pc1[order(abs(pc1$svec), decreasing = TRUE),,drop = FALSE]
top5 <- head(pc1, 5)
top5 <- cbind(top5, 
              annot = c("Peroxisome proliferator-activated receptor gamma coactivator 1-alpha",
                        "Centromere protein R",
                        "Dehydrogenase/reductase SDR family member 13",
                        "Calsyntenin-1",
                        "Guanylate-binding protein 1"))
top5
```

# References

