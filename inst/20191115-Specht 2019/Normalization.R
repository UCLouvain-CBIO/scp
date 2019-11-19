
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

# Filtering
dat_sc <- filterNA(dat_sc, pNA = 0.99) # same threshold as in Specht et al. 


####---- Specht et al. method (2019) ----####

# Step 1: normalize every run with the internal reference channel
lapply(run_sc, function(run){
  
})

# Step 2



