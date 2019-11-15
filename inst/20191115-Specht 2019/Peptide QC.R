
####---- Description ----####

# This script performs some basic QC on the peptide data from the Specht et al. 
# 2019 paper. 

####---- Setup environment ----####

# Load libraries and data

library(ggplot2)
library(gridExtra)
library(magrittr)
library(MSnbase)
library(scpdata)
library(tidyr)
library(vsn)
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


####---- Check QC samples ----####

# Check intensities per run and per channel


# Convert the original data to wide format
datw <- pivot_longer(data = dat0, cols = intensity.coln, 
                     names_to = "Channel", 
                     values_to = "Reporter.intensity")
# Plot the RI per channel to confirm empty and carrier channels
axis.txt <- array(data = paste0("TMT", 1:11), dim = 11,
                  dimnames = list(intensity.coln))
ggplot(data = datw, aes(x = Channel, y = Reporter.intensity)) +
  geom_violin() + scale_y_log10() +
  scale_x_discrete(labels = axis.txt,
                   limits = paste0("Reporter.intensity.", 0:10))



# Sd vs mean
# Plots row standard deviations versus row means
meanSdPlot(dat_qc, ranks = TRUE) # from MSnbase
# Note: the rank scale distributes the data more evenly along the x-axis and 
# allows a better visual assessment of the standard deviation as a function of 
# the mean.

show_heatmap(dat_qc, Log2 = TRUE, znorm = TRUE, trim = 5)
plotNA(dat_qc) # from MSnbase
scp_plotMissing(dat_qc) 


run <- pData(dat_qc)$run[814]
obj <- dat_qc[, pData(dat_qc)$run == run]
df <- cbind(pData(obj), t(exprs(obj)))
df <- pivot_longer(data = df, cols = -(1:6), values_to = "intensity") 
ggplot(data = df, aes(x = channel, y = intensity)) +
  geom_violin() + scale_y_log10() + ggtitle(run)
  

####---- Check the 1000 cells samples ----####


plotSCoPEset(obj = dat_1000, 
             run = unique(pData(dat_1000)$run)[1],
             phenotype = "cell_type")


####---- Check the 100 cells samples ----####

pdf("plots.pdf", onefile = TRUE)
for (i in seq(length(p))) {
  do.call("grid.arrange", p[[i]])  
}
dev.off()
pdf("./inst/20191115-Specht 2019/figs/SCoPE2 sets - 10 cells.pdf", 
    width = 10, height = 10, onefile = TRUE)
for(run in unique(pData(dat_100)$run)){
  p <- plotSCoPEset(obj = dat_100, 
               run = run,
               phenotype = "cell_type")
  do.call("grid.arrange", p)  
}
dev.off()



####---- Check the 10 cells samples ----####


plotSCoPEset(obj = dat_10, 
             run = unique(pData(dat_10)$run)[1],
             phenotype = "cell_type")


####---- Check the 1 cells samples ----####

plotSCoPEset(obj = dat_sc, 
             run = unique(pData(dat_sc)$run)[12],
             phenotype = "cell_type")

  





