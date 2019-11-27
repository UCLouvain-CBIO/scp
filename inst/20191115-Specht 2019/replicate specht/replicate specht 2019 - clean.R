
####---- Description ----#### 

# This script will replicate the analysis performed in Specht et al. 2019. The 
# processing of the MaxQuant output file will be reproduced but using clean  
# function implemented in the framework of MSnbase. 


####---- Load the data ----####

library(MSnbase)
setwd("./inst/20191115-Specht 2019/replicate specht/")

# Load meta data
samp <- read.csv("~/tmp/scpdata/inst/extdata/specht2019/annotation_fp60-97.csv", 
                 row.names = 1, check.names = F)
batch <- read.csv("~/tmp/scpdata/inst/extdata/specht2019/batch_fp60-97.csv", 
                  row.names = 1)

# Load the data and create an MSnSet
mq_file <- "~/tmp/scpdata/inst/extdata/specht2019/ev_updated.txt"
coln <- colnames(read.table(mq_file, header = TRUE, sep = "\t", nrow = 1))
sc0 <- readMSnSet2(file = mq_file, fnames = "id", sep = "\t", header = TRUE,
                   ecol = grep("intensity[.]\\d", coln, value = TRUE))
sc <- sc0


####---- Add a peptide sequence-charge field ----####

fData(sc)$sequence_charge <- paste0(fData(sc)$Modified.sequence, fData(sc)$Charge)


####---- Filter based on sample annotation ----####

# Keep only sc runs and runs that have associated sample metadata
sel <- !grepl("blank|col[12][0-9]", fData(sc)$Raw.file)
sc.runs <- unique(fData(sc)$Raw.file[sel])
# Note: the QC runs qre still present ! add "_QC_|" to the grepl to remove them
# Specht et al. select only FP97 and FP94 runs 
sc <- sc[!(sel & !grepl("FP9[47]", fData(sc)$Raw.file)), ]
# Remove sample not described in the metadata
sc <- sc[fData(sc)$Raw.file %in% colnames(samp), ]
# Note: this also removes the QC samples implicitely


####---- Filter based on identification measures ----####

# Remove the reverse hits (from decoy database) and contaminants
sc <- sc[fData(sc)$Reverse != "+", ]
sc <- sc[!grepl("^REV", fData(sc)$Leading.razor.protein), ]
sc <- sc[fData(sc)$Potential.contaminant != "+", ]
sc <- sc[fData(sc)$PIF > 0.8 & !is.na(fData(sc)$PIF), ]

# Remove spectra with poor identification confidence
# The PEP and q-values were updated using DART-ID
# TODO discuss with Laurent ! Why not deleting peptides with high FDR instead of 
# deleting the proteins for which all associated peptides have high FDR 
qprots <- unique(fData(sc)$Leading.razor.protein[fData(sc)$dart_qval < 0.01])
sc <- sc[fData(sc)$Leading.razor.protein %in% qprots, ]
sc <- sc[fData(sc)$dart_PEP < 0.02, ]

# Remove runs with insufficient identified peptides 
pep.t <- table(fData(sc)$Raw.file)
# sc <- sc[fData(sc)$Raw.file %in% names(pep.t[pep.t >= 300]),] 
thres <- sapply(names(pep.t), function(x){
  if(grepl("col19|col2[0-2]", x)){
    return(200)
  } else if(x %in% sc.runs){
    return(300)
  } else {
    return(0)
  }
})
sc <- sc[fData(sc)$Raw.file %in% names(pep.t[pep.t >= thres]),] 
# 25 runs were removed


####--- Filter based on sample to carrier ratio ----####

sc <- scp_filterSCR(sc, samples = 4:11, carrier = 1, thresh = 0.1,
                    sel = fData(sc)$Raw.file %in% sc.runs) # This is only performed single cell runs


####---- Numeric normalization ----####

# Normalize single cell runs to normalization channel
sc_nn <- scp_normalize_rRI(sc, ref_col = 2)

# Normalize single cell runs to carrier
# TODO: discuss with Laurent ! Here the carrier normalization is performed after reference normalization !! 
sc_nc <- scp_normalize_rRI(sc_nn, ref_col = 1)


####---- Format data to a peptide x sample matrix ----####

.keep <- c("Raw.file", "sequence_charge", "Modified.sequence", "Sequence", 
           "Length", "Proteins", "Leading.razor.protein", "Gene.names", 
           "Protein.names", "Mass")

## Reference normalized data

# Convert data to long format 
sc_nn <- scp_exprsToLong(sc_nn, f_sel = .keep)
# Remove duplicates
sc_nn <- sc_nn[!duplicated(fData(sc_nn)[, c("sequence_charge", "Raw.file", "channel")]), ]
# Convert data to wide format
sc_nn0 <- scp_exprsToWide(sc_nn)

## Carrier normalized data

# Convert data to long format 
sc_nc <- scp_exprsToLong(sc_nc, f_sel = .keep)
# Remove duplicates
sc_nc <- sc_nc[!duplicated(fData(sc_nc)[, c("sequence_charge", "Raw.file", "channel")]), ]
# Convert data to wide format
sc_nc0 <- scp_exprsToWide(sc_nc)

####---- Missing data cleaning ----####

# Replace zero's by a low value (needed for later)
# TODO for some strange reason the q30 and median RI filtering occurs on data with this type of trimming
sc_nn_lb <- scp_cleanMissing(sc_nn0, misVal = 0.001)
sc_nc_lb <- scp_cleanMissing(sc_nc0, misVal = 0.001)
# Replace zero's by NA (needed for later)
sc_nn <- scp_cleanMissing(sc_nn0, misVal = NA)
sc_nc <- scp_cleanMissing(sc_nc0, misVal = NA)


####---- Add sample metadata ----####

pd <- data.frame(row.names = sampleNames(sc_nn))
pd$run <- sapply(rownames(pd), function(x) strsplit(x, "-")[[1]][1])
pd$channel <- sapply(rownames(pd), function(x) strsplit(x, "-")[[1]][2])
pd$samp_type <- sapply(1:nrow(pd), function(i) samp[pd$channel[i], pd$run[i]])
pd <- cbind(pd, batch[pd$run, ])
pData(sc_nn) <- pd
pData(sc_nc) <- pd


####---- Keep only single cell data ----####

sel <- pData(sc_nn)$samp_type %in% c("sc_u","sc_m0", "sc_0") &
  pData(sc_nn)$run %in% sc.runs
sc_nn <- sc_nn[, sel]
sc_nc <- sc_nc[, sel]
sc_nc_lb <- sc_nc_lb[, sel]

####--- Filter based on data distributions ----####


sc_nnf <- scp_filterDD(sc_nn, sc_nc, sc_nc_lb, 
                       rowNorm = TRUE, npep = 6, 
                       qprobs = 0.3, q_thres = -2.5,
                       cv_thres = 0.43, 
                       median_thres = -1.3,
                       Plot = TRUE)
# Note1: I was not able to reproduce the CV calculation but had something very 
# similar
# Note2: npep or rowNorm has no big impact on the CV distribution. The 
# distribution changes a lot whether we use the carrier normalized or the 
# reference normalized data 
hist(exprs(sc_nnf), breaks = "FD", xlim = c(0, 2), col = "darkseagreen")


####---- Normalize rows and columns ----####

sc_nnfn <- scp_normalize_stat(sc_nnf, what = "column", fun = "/", stats = median)
sc_nnfn <- scp_normalize_stat(sc_nnfn, what = "row", fun = "/", stats = mean)

hist(exprs(sc_nnfn), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")


####---- Filter based on missingness ----####

sc_nnfn <- scp_filterNA(sc_nnfn, "row", pNA = 0.99)
sc_nnfn <- scp_filterNA(sc_nnfn, "column", pNA = 0.99)

hist(exprs(sc_nnfn), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")


####---- Log2 transform ----####

sc_final <- log(sc_nnfn, base = 2)

hist(exprs(sc_final), breaks = "FD", xlim = c(-2, 2), col = "darkseagreen")

# This sc_final is the peptide data as provided in Specth et al article.


