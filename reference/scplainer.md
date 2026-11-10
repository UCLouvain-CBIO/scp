# scplainer: linear models to understand mass spectrometry-based single-cell proteomics data

*scplainer*, standing for SCP-based Linear modelling Approach for
Interpretable aNd Explorable Results, is a principled and standardised
approach for extracting meaningful insights from SCP data. At its core,
the approach performs statistical modelling using linear regression.

The workflow starts from a
[SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object containing SCP data. The data is assumed to be
**log-transformed**. We advise to perform cell and feature quality
control to avoid that failed or outlying cells/feature distort the
results. We also recommend starting at the precursor or the
peptide-level, but the workflow allows protein-level data. Similarly,
the workflow is robust against for missing values, but it also allows
for data where missing values are imputed.

To learn how to import your data, we suggest reading the vignette:
[`vignette("read_scp", package = "scp")`](https://UCLouvain-CBIO.github.io/scp/articles/read_scp.md)

To learn how to process your data, we suggest reading the vignette:
[`vignette("scp", package = "scp")`](https://UCLouvain-CBIO.github.io/scp/articles/scp.md)

**Outline of the workflow**

1.  [ScpModel-Workflow](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-Workflow.md):
    performs the data modelling and filtering using linear regression.

2.  [ScpModel-VarianceAnalysis](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-VarianceAnalysis.md):
    investigate the contribution of each model variable to the data

3.  [ScpModel-DifferentialAnalysis](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-DifferentialAnalysis.md):
    assess the statistical significance of the differences observed
    between group of samples of interest.

4.  [ScpModel-ComponentAnalysis](https://UCLouvain-CBIO.github.io/scp/reference/ScpModel-ComponentAnalysis.md):
    visually explore the data captured by each model variable.

Once the data are modelled and explored, the filtered, normalised and
batch-corrected data can be retrieved for further downstream analysis,
such as clustering or trajectory inference.

You can find a demonstration of the *scplainer* workflow in a dedicated
vignette:
[`vignette("scp_data_modelling", package = "scp")`](https://UCLouvain-CBIO.github.io/scp/articles/scp_data_modelling.md)

## References

scplainer: using linear models to understand mass spectrometry-based
single-cell proteomics data Christophe Vanderaa, Laurent Gatto bioRxiv
2023.12.14.571792; doi: https://doi.org/10.1101/2023.12.14.571792.

## Author

Christophe Vanderaa, Laurent Gatto
