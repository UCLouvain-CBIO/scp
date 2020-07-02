
## This script creates a small toy dataset to showcase the `scp` pipeline

## The dataset is a subset of the SCoPE2 dataset (Specht et al. 2019, BioRXiv)
## https://www.biorxiv.org/content/10.1101/665307v3
## The dataset is contained in the `scpdata` package

library(scp)
library(scpdata)
library(tidyverse)
data("specht2019v2")

## Randomly sample 1 SCoPE2 set from each possible combination of batch variables
set.seed(1000)
specht2019v2 %>%
  colData %>%
  data.frame %>%
  group_by(lcbatch, sortday, digest) %>%
  sample_n(1) %>%
  pull(Set) ->
  sampledRuns
scp1 <- specht2019v2[, , sampledRuns]

## Sample 150 proteins from each run
set.seed(1000)
rowDataToDF(scp1, seq_along(scp1), "protein") %>%
  data.frame %>% 
  rownames_to_column("psm") %>%
  group_by(.assay) %>%
  filter(protein %in% sample(unique(protein), 150)) ->
  l
l <- split(l$psm, f = l$.assay)
scp1 <- scp1[l, ]

## Remove all place holder samples
lapply(experiments(scp1), function(e) colSums(is.na(assay(e))) != nrow(e) ) %>%
  unname %>%
  unlist -> 
  cols
cols <- names(cols)[cols]
scp1 <- scp1[, cols, ]

## Save the data 
format(object.size(scp1), units = "MB", digits = 2)
save(scp1, file = file.path("data/scp1.rda"), 
     compress = "xz", compression_level = 9)

  