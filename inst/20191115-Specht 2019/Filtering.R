
####---- Description ----####

# This script investigates the impact of filtering (cells and peptides) on the 
# SCP data from Specht et al. 2019. The methodology used is inpspired from the 
# code coming with the article. Note that some filtering was already performed 
# when parsing the data into a data matrix. This includes: 
# - Remove contaminants identified by MaxQuant using a decoy database
# - Remove contaminated spectra where PIF >= 0.8
# - Keep peptides with PEP > 
# The previous filtering was performed on single peptide observations, whereas 
# the filtering hereunder are based on population-based heuristics. 


####---- Setup environment ----####

# Load libraries, functions and data
source("./R/utils.R")
dat <- specht2019_peptide2
# Subset the data of interest
run_all <- unique(pData(dat)$run)
run_sc <- run_all[! grepl("_QC_|col19|col2[0-4]", run_all)]
dat_sc <- dat[, pData(dat)$run %in% run_sc] # Keep only single cell runs


####---- Filter the experiments ----####

# In the article only single-cell experiments with code 94 and 97 are kept. 
# Indeed, in the QC script we could see that the other samples contained more 
# missing data. (see "QC - effect of digestion on missingness.png").
dat_sc <- dat_sc[, pData(dat_sc)$digest %in% c("Q", "N")] 


####---- Filter the 

