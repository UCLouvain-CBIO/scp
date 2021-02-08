library(tidyverse)
library(scp)

## Create an MaxQuant output example file

## To run this code, you need: 
##  1. First run the `scp1.R` script in the `scp/inst/script` folder
data("scp1")
##  2. Downloaded the `raw.RData` file that contains the SCoPE2 data
##     set. The file can be downloaded here:
##     https://drive.google.com/file/d/1sN8CkEnAh3z0CFKegzqmwCQSUZt2mn9M/view?usp=sharing
load("../.localdata/SCP/specht2019/v3/raw.RData")
##  3. Download the `evidence_unfiltered.csv` file. This table 
##     contains additional blank samples that are used as a test case
##     in the vignette. The file can be downloaded here:
##     https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx

## The `ev` table was loaded from the `raw.Rdata` file
ev %>% 
    select(-c("lcbatch", "sortday",  "digest")) %>%
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
    mutate_if(is.logical, as.character) %>%
    ## keep sets selected in scp1
    filter(Set %in% names(scp1) &
               ## keep only a few proteins
               protein %in% rowDataToDF(scp1, 1:3, "protein")$protein |
               grepl("FP97_blank_01", Set) ) ->
    samples

## Add the blank sample
read.csv("../.localdata/SCP/specht2019/v2/evidence_unfiltered.csv",
                header = TRUE) %>%
    select(-c("X", "X1", "lcbatch", "sortday",  "digest")) %>%
    ## extract one of the blank samples
    filter(grepl("FP97_blank_01", Raw.file)) %>%
    ## channel naming should be consistent with metadata
    rename_all(gsub, 
               pattern = "^Reporter[.]intensity[.](\\d*)$", 
               replacement = "RI\\1") %>%
    ## Remove the empty TMT12-16
    select(-paste0("RI", 12:16)) %>%
    ## MS set should be consistent with metadata and other data
    dplyr::rename(Set = Raw.file, 
                  peptide = modseq,
                  protein = Leading.razor.protein) %>%
    ## Remove "X" at start of batch 
    mutate(Set = gsub("^X", "", Set)) ->
    blank
## Adjust data types 
for(col in grep("", colnames(blank), value = TRUE))
    blank[, "remove"] <- as.character(blank[, "remove"])
for(col in grep("Deamidation..NQ.$", colnames(blank), value = TRUE))
    blank[, col] <- as.numeric(blank[, col])
for(col in grep("Reporter.intensity.corrected.", colnames(blank), value = TRUE))
    blank[, col] <- as.double(blank[, col])

mqScpData <- bind_rows(samples, blank)
format(object.size(mqScpData), units = "MB", digits = 2)
save(mqScpData, 
     file = file.path("data/mqScpData.rda"), 
     compress = "xz", 
     compression_level = 9)

## Create the associated annotation file

## The `design` and `batch` tables were loaded from the `raw.RData` 
## file. We clean somewhat the sample metadata so that it meets the 
## requirements for `scp::readSCP`. The cell annotation and batch 
## annotation are merge into a single table
inner_join(x = design %>% 
               mutate(Set = sub("^X", "", Set)) %>%
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
               mutate(Set = sub("^X", "", Set)) %>%
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
