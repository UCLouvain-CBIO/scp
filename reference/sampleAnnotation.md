# Single cell sample annotation

A data frame with 48 observations on the following 6 variables.

- Set: a character vector

- Channel: a character vector

- SampleType: a character vector

- lcbatch: a character vector

- sortday: a character vector

- digest: a character vector

## Usage

``` r
data("sampleAnnotation")
```

## Format

An object of class `data.frame` with 64 rows and 6 columns.

## Details

\##' The dataset is a subset of the SCoPE2 dataset (version 2, Specht et
al. 2019, [BioRXiv](https://www.biorxiv.org/content/10.1101/665307v3)).
The input files `batch.csv` and `annotation.csv` were downloaded from a
[Google Drive
repository](https://drive.google.com/drive/folders/1VzBfmNxziRYqayx3SP-cOe2gu129Obgx).
The two files were loaded and the columns names were adapted for
consistency with `mqScpData` table (see
[`?mqScpData`](https://UCLouvain-CBIO.github.io/scp/reference/mqScpData.md)).
The two tables were filtered to contain only sets present in
“mqScpData`. The tables were then merged based on the run ID, hence merging the sample annotation and the batch annotation. Finally, annotation for the blank run was added manually. The data is stored as a `data.frame\`.

## See also

[`readSCP()`](https://UCLouvain-CBIO.github.io/scp/reference/readSCP.md)
to see how this file is used.
