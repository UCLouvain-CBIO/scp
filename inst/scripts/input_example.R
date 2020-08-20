library(tidyverse)
library(scp)

## Create an MaxQuant output example file

## The dataset is a subset of the SCoPE2 dataset (Specht et al. 2019,
## BioRXiv) https://www.biorxiv.org/content/10.1101/665307v3 The
## datatable was downloaded from
## https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx
## To run this code, you need to first run scp1.R
read.csv("inst/extdata/evidence_unfiltered.csv", sep = ",", header = TRUE)
ev %>% 
  select(-c("X", "X1", "lcbatch", "sortday",  "digest")) %>%
  ## channel naming should be consistent with metadata
  rename_all(gsub, 
             pattern = "^Reporter[.]intensity[.](\\d*)$", 
             replacement = "RI\\1") %>%
  ## MS set should be consistent with metadata and other data
  dplyr::rename(Set = Raw.file, 
                peptide = modseq,
                protein = "Leading.razor.protein") %>%
  ## Remove "X" at start of batch 
  mutate(Set = gsub("^X", "", Set)) %>%
  ## keep only single cell runs 
  filter(!grepl("col\\d{2}|arrier|Ref|Master|SQC|blank", Set) &
           ## remove experimental sets concurrent with low mass spec performance
           !grepl("FP9[56]|FP103", Set) &
           ## keep only sets selected in scp1
           Set %in% names(scp1) &
           ## keep only a few proteins
           protein %in% rowDataToDF(scp1, 1:3, "protein")$protein) ->
  mqFile
format(object.size(mqFile), units = "MB", digits = 2)
save(mqFile, file = file.path("data/mqFile.rda"), 
     compress = "xz", compression_level = 9)

## Create the associated annotation file
## The datatables are downloaded from https://www.biorxiv.org/content/10.1101/665307v3

cells <- read.csv("inst/extdata/annotation.csv", check.names = FALSE)
batch <- read.csv("inst/extdata/batch.csv", check.names = FALSE)
## Clean the sample metadata so that it meets the requirements for
## `scp::readSCP`. The cell annotation and batch annotation are merge into a 
## table
inner_join(x = cells %>% 
             pivot_longer(-Set, 
                          names_to = "Channel", 
                          values_to = "SampleAnnotation") %>%
             mutate_all(as.character), 
           y = batch %>% 
             dplyr::rename(Set = set) %>%
             mutate_all(as.character),
           by = "Set") -> sampleAnnotation
format(object.size(sampleAnnotation), units = "MB", digits = 2)
save(sampleAnnotation, file = file.path("data/sampleAnnotation.rda"), 
     compress = "xz", compression_level = 9)


