
#### MINIMAL PROCESSING

## This script performs the minimal processing on the leduc2022_pSCoPE
## dataset. 

####---- Loading libraries and preparing data ----####

## libraries
library("scp")
library("scpdata")
library("ggplot2")
library("patchwork")
library("ensembldb")
library("EnsDb.Hsapiens.v86")
library("dplyr")

## data
leduc <- leduc2022_pSCoPE()

## keep only required data
assaysToRemove <- c(
    "peptides", "peptides_log", "proteins_norm2", "proteins_processed"
)
leduc <- removeAssay(leduc, assaysToRemove)
requiredRowData <- c(
    "Sequence", "Leading.razor.protein.symbol", 
    "Leading.razor.protein.id", "Reverse", "Potential.contaminant",
    "Leading.razor.protein", "PIF", "dart_qval"
)
leduc <- selectRowData(leduc, requiredRowData)

## format missing values
leduc <- zeroIsNA(leduc, i = names(leduc))

####---- Feature quality control ----####

leduc <- computeSCR(
    leduc, names(leduc), colvar = "SampleType", 
    samplePattern = "Mel|Macro", carrierPattern = "Carrier",
    sampleFUN = "mean", rowDataName = "MeanSCR"
)
df <- data.frame(rbindRowData(leduc, names(leduc)))
df$ContaminantOrReverse <- df$Reverse != "+" & 
    df$Potential.contaminant != "+" & 
    !grepl("REV|CON", df$Leading.razor.protein)
## Contaminant plot
ggplot(df) + 
    aes(x = ContaminantOrReverse) +
    geom_bar() +
    ## PIF plot
    ggplot(df) + 
    aes(x = PIF) +
    geom_histogram() +
    ## q-value plot
    ggplot(df) + 
    aes(x = log10(dart_qval)) +
    geom_histogram() +
    ## mean SCR plot
    ggplot(df) + 
    aes(x = log10(MeanSCR)) +
    geom_histogram()
leduc <- filterFeatures(
    leduc, ~ Reverse != "+" &
        Potential.contaminant != "+" &
        !grepl("REV|CON", Leading.razor.protein) &
        !is.na(PIF) & PIF > 0.6 &
        dart_qval < 0.01 &
        !is.na(MeanSCR) & MeanSCR < 0.05
)

####---- Sample quality control ----####

## Number of detected peptides per sample
leduc <- countUniqueFeatures(
    leduc, i = names(leduc), groupBy = "Sequence",
    colDataName = "NumberPeptides"
)
## Median intensity per sample
# for (i in names(leduc)) {
#     logAssay <- log2(assay(leduc[[i]]))
#     meds <- colMedians(logAssay, na.rm = TRUE)
#     colData(leduc)[names(med), "MedianIntensity"] <- meds
# }
MedianIntensity <- lapply(experiments(leduc), function(x) {
    out <- colMedians(log(assay(x)), na.rm = TRUE)
    names(out) <- colnames(x)
    out
})
names(MedianIntensity) <- NULL
MedianIntensity <- unlist(MedianIntensity)
colData(leduc)[names(MedianIntensity), "MedianIntensity"] <- MedianIntensity
## Median CV per sample
leduc <- medianCVperCell(
    leduc, i = names(leduc), groupBy = "Leading.razor.protein.symbol",
    nobs = 3, na.rm = TRUE, colDataName = "MedianCV", norm = "SCoPE2"
)

data.frame(colData(leduc)) |>
    ggplot() +
    aes(
        y = MedianIntensity,
        x = NumberPeptides,
        color = MedianCV,
        shape = SampleType
    ) +
    geom_point(size = 2) +
    scale_color_continuous(type = "viridis")

leduc$passQC <- !is.na(leduc$MedianCV) & leduc$MedianCV < 0.6 &
    leduc$MedianIntensity > 6 & leduc$MedianIntensity < 8 &
    leduc$NumberPeptides > 750 &
    grepl("Mono|Mel", leduc$SampleType)
leduc <- subsetByColData(leduc, leduc$passQC)

####---- Building the peptide matrix ----####

## Aggregate PSMs to peptides
peptideAssays <- paste0("peptides_", names(leduc))
leduc <- aggregateFeatures(leduc,
                           i = names(leduc),
                           fcol = "Sequence",
                           name = peptideAssays,
                           fun = colMedians,
                           na.rm = TRUE)
## Apply majority vote for peptide to protein mapping
ppMap <- rbindRowData(leduc, i = grep("^pep", names(leduc))) |> 
    data.frame() |>
    group_by(Sequence) |> 
    ## The majority vote happens here
    mutate(Leading.razor.protein.symbol =
               names(sort(table(Leading.razor.protein.symbol),
                          decreasing = TRUE))[1],
           Leading.razor.protein.id =
               names(sort(table(Leading.razor.protein.id),
                          decreasing = TRUE))[1]) |> 
    dplyr::select(Sequence, Leading.razor.protein.symbol, Leading.razor.protein.id) |>
    dplyr::filter(!duplicated(Sequence, Leading.razor.protein.symbol))
consensus <- lapply(peptideAssays, function(i) {
    ind <- match(rowData(leduc)[[i]]$Sequence, ppMap$Sequence)
    DataFrame(Leading.razor.protein.symbol =
                  ppMap$Leading.razor.protein.symbol[ind],
              Leading.razor.protein.id = 
                  ppMap$Leading.razor.protein.id[ind])
})
names(consensus) <- peptideAssays
rowData(leduc) <- consensus
## Join peptide assays 
leduc <- joinAssays(leduc, i = peptideAssays, 
                    name = "peptides")

proteinIds <- rowData(leduc)[["peptides"]]$Leading.razor.protein.id
## Add gene name information
proteinConversionDf <- transcripts(
    EnsDb.Hsapiens.v86, 
    columns = "gene_name",
    return.type = "data.frame",
    filter = UniprotFilter(proteinIds)
)
matchedIndex <- match(proteinIds, proteinConversionDf$uniprot_id)
geneName <- proteinConversionDf$gene_name[matchedIndex]
rowData(leduc)[["peptides"]]$gene <- geneName

####---- Log-transformation ----####

leduc <- logTransform(leduc, i = "peptides", name = "peptides_log")

####---- Subset cells ----####

## Take the first 3 and the last 3 MS batches
leduc_minimal <- getWithColData(leduc, "peptides_log")
sets <- unique(leduc_minimal$Set)
sel <- c(head(sets, 3), tail(sets, 3))
leduc_minimal <- leduc_minimal[, leduc_minimal$Set %in% sel]

####---- Subset peptides ----####

## Remove all NAs
leduc_minimal <- filterNA(leduc_minimal, pNA = 0.9999)
## Randomly select 200 peptides
set.seed(1234)
sel <- sample(1:nrow(leduc_minimal), 200, replace = FALSE)
leduc_minimal <- leduc_minimal[sel, ]

####---- Example model ----####

f <- ~ 1 + ## intercept
    Channel + Set + ## batch variables
    MedianIntensity +## normalization
    SampleType ## biological variable
leduc_minimal <- scpModelWorkflow(leduc_minimal, formula = f)

####---- Save results ----####

saveDir <- "~/PhD/asca-scp/package/data/"
save(
    leduc_minimal, file = paste0(saveDir, "leduc_minimal.rda"),
    compress = "xz", compression_level = 9
)
