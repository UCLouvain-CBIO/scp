data("scp1")

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

