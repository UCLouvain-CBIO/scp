
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Clear all previously loaded data
rm(list=ls())


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

### Hard-coded google drive directory to data
path.t<-"G:/My Drive/MS/SCoPE/mPOP/"
#path.t<-"~/GoogleDrv/MS/SCoPE/mPOP/"


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

### Code settigngs: 

# Remove peptides that are more abundant that 10% of the carrier channel from the single cell runs
remove_abundant_peptides<-T

# Use dart-enhanced spectra or spectra
dart_or_spectra<-"dart"

# Normalize the 10 x 10 cell and 10 x 100 cell runs to their median value
norm_10_100_median<-T

# Only look at sets FP94+97 (same experiment, sorted onto two different 384-well plates)
only_fp94_97<-T

# Load the data fresh from search output, else use saved version (.RData, which is quicker to load)
load_fresh<-F

# Correct for isotopic cross-contamination from TMT11 -- this will not work properly because MQ did the correction already 
# for the data I am using
iso_cor<-F

# Minimum number of peptides observed to consider an experiment worth doing further analysis: 
thres_ids_sc <- 300
thres_ids_c10c100 <- 200

# Figure dimensions in inches
w.t<-5
h.t<-5

# Figure dir: 
#save.path<-"G:/My Drive/2018_mPOP/2018_mPOP/Figs/fig3/"

# Filter to less than X missing data per column or row: 
na.col<-0.99
na.row<-0.99

# Imputation
k.t<-3

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Execute all scripts in this order, some later scripts contain variables created in previous scripts

source("code/mPOP_source_functions.R") # 
source("code/import.R") # 
source("code/score.R") # 
source("code/data_transformations.R") # 
source("code/missing_data.R") 

source("code/bulk_markers.R") # 
source("code/dimensional_reduction_colored.R") # 
source("code/spectral_ordering_mac_v_mono.R") # 
source("code/spectral_ordering_mac.R") #

source("code/coverage.R") # 
source("code/density_corr.R") # 
source("code/eigenvector_density.R") # 
source("code/Seurat.R") # 
source("code/ion_counting_v4.R") # 
source("code/figure_2b_v3.R") # 


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

