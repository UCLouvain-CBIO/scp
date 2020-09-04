library(tidyverse)
library(scp)

## Create an MaxQuant output example file

## The dataset is a subset of the SCoPE2 dataset (Specht et al. 2019,
## BioRXiv) https://www.biorxiv.org/content/10.1101/665307v3 The
## datatable was downloaded from
## https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx
## To run this code, you need to first run scp1.R
ev <- read.csv("inst/extdata/evidence_unfiltered.csv", sep = ",", header = TRUE) 
  
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
  ## keep sets selected in scp1
  filter(Set %in% names(scp1) &
           ## keep only a few proteins
           protein %in% rowDataToDF(scp1, 1:3, "protein")$protein |
           grepl("FP97_blank_01", Set) ) ->
  mqScpData
format(object.size(mqScpData), units = "MB", digits = 2)
save(mqScpData, file = file.path("data/mqScpData.rda"), 
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
                          values_to = "SampleType") %>%
             mutate(SampleType = recode(SampleType, 
                                        sc_0 = "Blank",
                                        sc_u = "Monocyte",
                                        sc_m0 = "Macrophage",
                                        unused = "Unused",
                                        norm = "Reference",
                                        reference = "Reference",
                                        carrier_mix = "Carrier")) %>%
             mutate_all(as.character), 
           y = batch %>% 
             dplyr::rename(Set = set) %>%
             mutate_all(as.character),
           by = "Set") %>%
  filter(Set %in% mqScpData$Set) %>%
  ## Add the metadata for the blank sample
  add_row(Set = grep("blank", mqScpData$Set, value = TRUE) %>%
            unique %>%
            rep(16),
          Channel = paste0("RI", 1:16),
          SampleType = "Blank",
          lcbatch = "LCA10",
          sortday = "s8") %>%
  data.frame ->
  sampleAnnotation
  
format(object.size(sampleAnnotation), units = "MB", digits = 2)
save(sampleAnnotation, file = file.path("data/sampleAnnotation.rda"), 
     compress = "xz", compression_level = 9)
