data("scp1")

test_that("rowDataToDF", {
    ## Correct use 
    ## single assay
    test <- rowDataToDF(scp1, i = "190321S_LCA10_X_FP97AG", vars = "peptide")
    expect_identical(nrow(test), nrow(scp1[[1]]))
    expect_identical(dimnames(test), list(rownames(scp1[[1]]), c("peptide", ".assay", ".rowname")))
    ## Multiple assays
    test <- rowDataToDF(scp1, i = 1:3, vars = "peptide")
    expect_identical(nrow(test),
                     as.integer(sum(dims(scp1)[1, 1:3])))
    expect_identical(dimnames(test),
                     list(unlist(lapply(1:3, function(x) rownames(scp1[[x]]))),
                          c("peptide", ".assay", ".rowname")))
    ## Variable not found
    expect_error(rowDataToDF(scp1, i = 1:3, vars = "foo"),
                 regexp = "rowData variable\\(s\\) not found")
})
  

test_that("transferColData", {
    ## Correct Use
    test <- transferColDataToAssay(scp1, i = "peptides")
    expect_identical(colData(test[["peptides"]]), 
                     colData(test)[colnames(test[["peptides"]]), ])
    ## Expect message
    expect_message(transferColDataToAssay(test, i = "peptides"),
                   regexp = "colData is already present in assay")
    
})


test_that("aggregateFeaturesOverAssays", {
    ## Correct use
    test <- aggregateFeaturesOverAssays(scp1, i = 1:3, fcol = "peptide",
                                        name = paste0("peptides", 1:3),
                                        fun = colSums)
    expect_identical(test[["peptides1"]],
                     aggregateFeatures(scp1, i = 1, fcol = "peptide", name = "peptides1",
                                       fun = colSums)[["peptides1"]])
    al <- assayLink(test, "peptides1")
    expect_identical(al@name, "peptides1")
    expect_identical(al@from, "190321S_LCA10_X_FP97AG")
    expect_identical(al@fcol, "peptide")
    expect_identical(al@hits@from, 1:nrow(scp1[[1]]))
    expect_identical(al@hits@to, as.integer(as.factor(rowData(scp1[[1]])$peptide)))
    ## Error: i and names must have same size
    expect_error(aggregateFeaturesOverAssays(scp1, i = 1:3, fcol = "peptide",
                                             name = "foo"),
                 regexp = "'i' and 'name' must have same length")
    ## Error: i and fcoli must have same size
    expect_error(aggregateFeaturesOverAssays(scp1, i = 1:3, name = 1:3,
                                             fcol = rep("peptide", 2)),
                 regexp = "'i' and 'fcol' must have same length")
    ##
    
})
    