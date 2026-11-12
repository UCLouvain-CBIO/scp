# Minimally processed single-cell proteomics data set

A `SingleCellExperiment` object that has been minimally processed. The
data set is published by Leduc et al. 2022 (see references) and
retrieved using
[`scpdata::leduc2022_pSCoPE()`](https://UCLouvain-CBIO.github.io/scpdata/reference/leduc2022_pSCoPE.html).
The data processing was conducted with `QFeatures` and `scp`. Quality
control was performed, followed by building the peptide data and
log2-transformation. To limit the size of the data, only cells
associated to the 3 first and 3 last MS acquisition runs were kept. For
the same reason, 200 peptides were randomly sampled. Therefore, the data
set consists of 200 peptides and 73 cells. Peptide annotations can be
retrieved from the `rowData` and cell annotations can be retrieved from
the `colData`.

## Usage

``` r
data("leduc_minimal")
```

## Format

An object of class `SingleCellExperiment` with 200 rows and 73 columns.

## Quality control

Any zero value has been replaced by NA.

A peptide was removed from the data set if:

- it matched to a decoy or contaminant peptide

- it had an parental ion fraction below 60 \\

- it had a DART-ID adjusted q-value superior to 1\\

- it had an average sample to carrier ratio above 0.05

A cell was removed from the data set if:

- it had a median coefficient of variation superior to 0.6

- it had a log2 median intensity outside (6, 8)

- it contained less than 750 peptides

## Building the peptide matrix

PSMs belonging to the same peptide were aggregating using the median
value. Some peptides were mapped to a different protein depending on the
MS acquisition run. To solve this issue, a majority vote was applied to
assign a single protein to each peptide. Protein IDs were translated
into gene symbols using the `ensembldb` package.

## References

Leduc, Andrew, R. Gray Huffman, Joshua Cantlon, Saad Khan, and Nikolai
Slavov. 2022. “Exploring Functional Protein Covariation across Single
Cells Using nPOP.” Genome Biology 23 (1): 261.

## Author

Christophe Vanderaa, Laurent Gatto
