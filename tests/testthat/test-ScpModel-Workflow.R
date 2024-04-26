## ---- ScpModel workflow ----

## Internal function that creates a minimal SE object as expected by
## scplainer for unit testing ScpModel class methods
## @param nr Number of rows
## @param nc Number of columns
.createMinimalData <- function(nr = 10, nc = 5) {
    require("SummarizedExperiment")
    a <- matrix(1, nr, nc)
    rownames(a) <- letters[1:nr]
    colnames(a) <- LETTERS[1:nc]
    SummarizedExperiment(assays = List(assay = a))
}

test_that("scpModelWorkflow", {
    require("SummarizedExperiment")
    se <- SummarizedExperiment()
    ## Empty SE = error
    expect_warning(expect_error(
        scpModelWorkflow(object = se, formula = ~ 1),
        "!.*null.colnames.object.*is not TRUE"
    ))
    ## Test formula with no explanatory variables = error
    se <- .createMinimalData()
    expect_warning(expect_error(
        scpModelWorkflow(se, formula = ~ 1),
        "!.*null.colnames.object.*is not TRUE"
    ))

    ## Test formula with no response variable

    ## Test i as numeric (when assays are unnamed)

    ## Test i as numeric (when assays are named)

    ## Test i as character
    ## Test i as logical

    ## Test no name

    ## Test verbose

    ## Metadata element is already present = warning
    metadata(se)$foo <- "bar"

})