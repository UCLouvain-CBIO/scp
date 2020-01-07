
####---- Description ----####

## The script generates data used for examples


####---- Maxquant example data ----####

## The MaxQuant data is loaded into an MSnSet object and subset for reducing the
## object size. The MaxQuant data used here is the evidence file supplied by the
## Slavov lab, from the Specht et al. 2019 preprint (BioRxiv).

if(F){ # don't run at every build
  library(MSnbase)
  mq_file <- "./inst/20191115-Specht 2019/extdata/ev_updated.txt"
  coln <- colnames(read.table(mq_file, header = TRUE, sep = "\t", nrows = 1))
  mq <- readMSnSet2(file = mq_file, fnames = "id", header = TRUE, sep = "\t",
                    nrows = 1000, ecol = grep("intensity[.]\\d",
                                              coln, value = TRUE))
  save(mq, file = "./data/mq.rda")
}

#' MaxQuant example data
#'
#' Subset of the MaxQuant evidence file from the Specht et al. 2019
#' preprint.
#'
#' @usage
#' data("mq")
#'
#' @format
#' an MSnSet with MS intensities for 1000 peptides in 11 TMT channels
#'
#' @source
#' The original data can be downloaded from the
#' \href{https://scope2.slavovlab.net/docs/data}{Slavov Lab} website, in the
#' \href{https://scope2.slavovlab.net/docs/data}{Data repository 1}.
#'
#' @references
#' Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov.
#' 2019. “High-Throughput Single-Cell Proteomics Quantifies the Emergence of
#' Macrophage Heterogeneity.” bioRxiv (\href{https://doi.org/10.1101/665307}
#' {DOI}).
#'
#' @docType data
#'
#' @keywords datasets
"mq"

