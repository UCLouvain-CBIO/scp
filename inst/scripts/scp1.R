
## This script creates a small toy dataset to showcase the `scp` pipeline

## The dataset is a subset of the SCoPE2 dataset (Specht et al. 2019, BioRXiv)
## https://www.biorxiv.org/content/10.1101/665307v3
## The dataset is contained in the `scpdata` package

library(scp)
library(scpdata)
library(tidyverse)
library(QFeatures)

scp1 <- specht2019v3()

## Randomly sample 1 SCoPE2 set from each possible combination of batch variables
set.seed(1000)
scp1[, , 1:177] %>%
    colData %>%
    data.frame %>%
    filter(lcbatch %in% c("LCA9", "LCA10", "LCB3")) %>%
    group_by(lcbatch) %>%
    sample_n(1) %>%
    pull(Set) ->
    sampledRuns
scp1 <- scp1[, , sampledRuns]

## Sample 100 proteins from each run
set.seed(1000)
rbindRowData(scp1, seq_along(scp1)) %>%
    data.frame %>% 
    select(peptide, protein) %>% 
    rownames_to_column("psm") %>%
    group_by(peptide) %>%
    mutate(nprots = length(unique(protein))) %>%
    filter(nprots == 1) %>%
    group_by(assay) %>%
    filter(protein %in% sample(unique(protein), 100)) ->
    l
l <- split(l$psm, f = l$assay)
scp1 <- scp1[l, ]

## Aggregate to peptides
scp1 <- aggregateFeaturesOverAssays(scp1, i = names(scp1), fcol = "peptide",
                                    name = paste0("pep_", names(scp1)),
                                    fun = colMeans, na.rm = TRUE)

## Join peptides 
scp1 <- joinAssays(scp1, i = grep("pep_", names(scp1)), name = "peptides")

## Aggregate to proteins
scp1 <- aggregateFeatures(scp1, i = "peptides", fcol = "protein",
                          name = "proteins", fun = colMedians, na.rm = TRUE)

## Keep only the PSM, joined peptide and protein data
scp1 <- scp1[, , !grepl("pep_", names(scp1))]

## Save the data 
format(object.size(scp1), units = "MB", digits = 2)
save(scp1, file = file.path("data/scp1.rda"), 
     compress = "xz", compression_level = 9)

