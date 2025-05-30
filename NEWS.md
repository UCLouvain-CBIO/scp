# scp 1.19

## scp 1.19.3

- Nothing yet.

## scp 1.19.2

- feat: improved error message when pattern is not found in
  computeSCR
- feat: use typeMetadata (based on QFeatues 1.19.1)

## scp 1.19.1

- fix: moved from the deprecated `longFormat()` to `longForm()`
- doc: fixed readSingleCellExperiment example (see #94)
- doc: minor documentation recompilation

## scp 1.19.0

- New Bioconductor 3.22 release

# scp 1.18

## scp 1.18.0

- New Bioconductor 3.21 release

# scp 1.17

## scp 1.17.0

- New Bioconductor 3.21 devel

## scp 1.17.1

- Deprecate `aggregateFeaturesOverAssays`, use
  `QFeatures::aggregateFeatures` instead.

## scp 1.17.2

- Recommend the use of `SummarizedExperiment` instead of
`SingleCellExperiment` for performance reasons.
- All functions are now compatible with `SummarizedExperiment`.
- `data("scp1")` now contains `SummarizedExperiment` objects.
- By default, `readSCP` and `readSCPfromDIANN` now create
QFeatures objects that contain `SummarizedExperiment` objects.
- `readSCP` and `readSCPfromDIANN` have a new parameter:
`experimentsAsSce` that indicate if the QFeatures should
contain `SingleCellExperiment`, default is `FALSE`.
- `readSCP` and `readSCPfromDIANN` now return a QFeatures of type
"scp".
- `scpModelWorkflow` now automatically convert a `SummarizedExperiment`
to a `SingleCellExperiment`.

# scp 1.16

## scp 1.16.0

- New Bioconductor 3.20 release

# scp 1.15

## scp 1.15.2

- fix: fixed x-axis direction annotation for volcano plot on contrast
- fix: solved bug in DA when missing contrast level in modelled
  feature (issue #65).
- Add link to Leduc SCP.replication vignette.

## scp 1.15.1

- test: added unit tests for scplainer: ScpModel-Class,
  scpModelFit-Class, ScpModel-Workflow

## scp 1.15.0

- New Bioconductor 3.20 (devel) release

# scp 1.14

## scp 1.14.0

- New Bioconductor 3.19 (stable) release

# scp 1.13

## scp 1.13.5

- fix: fixed small error in degrees of freedom
- fix: break workflow upon infinite values

## scp 1.13.4

- fix: first drop variable before centering numerical variables (see #54)

## scp 1.13.3

- Refactor `readSCP()` and `readSCPFromDIANN()` to use new QFeatures
  implementations.

## scp 1.13.2

- New *scplainer* workflow and citation
- Add `addReducedDims()` (see #52)
- fix: logFC and associated SE are now correctly computed

## scp 1.13.1

- docs: fixed vignette about reporting missing values
- docs: add a QFeatures figure to the nutshell vignette

## scp 1.13.0

- New Bioconductor 3.18 (devel) release

# scp 1.12

## scp 1.12.0

- New Bioconductor 3.18 (stable) release

## scp 1.11.3

- Major feat: added SCP data modelling workflow + documentation
- feat: added readSCPfromDIANN() that creates a QFeatures object from
  DIANN output tables.

## scp 1.11.2

- Nothing yet.

## scp 1.11.2

- feat: added reportMissingValues(), jaccardIndex(),
  cumulativeSensitivityCurve() and predictSensitivity() to facilitate
  reporting missing values. The vignette is also adapted with the new
  functionality.
- docs: created vignette about reporting missing values in SCP
- fix failing unit test.

## scp 1.11.1

- Updated citation

## scp 1.11.0

- New Bioconductor 3.18 (devel) release

# scp 1.10

## scp 1.10.0

- New Bioconductor 3.17 (stable) release

# scp 1.9

## scp 1.9.2

- Updated citation

## scp 1.9.1

- Fix minor typo in `readSCP()` man page

## scp 1.9.0

- New Bioconductor 3.17 (devel) release

# scp 1.8

## scp 1.8.0

- New Bioconductor 3.16 (stable) release

# scp 1.7

## scp 1.7.5

- Updated CITATION
- Added sticker

## scp 1.7.4

- Updated CITATION

## scp 1.7.3

- refactor: package complies with BiocCheck
- docs: fixed bug in vignette

## scp 1.7.2

- Add CC-BY-SA license for vignettes.

## scp 1.7.1

- refactor: removed deprecated function `rowDataToDF()`
- tests: fixed some tests failing because of SCE version differences.
- feat: users can now specify `sep` when sample names are automatically
  generated.

## scp 1.7.0

- New devel (Bioc 3.16)

# scp 1.6

- New stable release (Bioc 3.15)

# scp 1.5

## scp 1.5.1

- Added CITATION <2021-10-29>

## scp 1.5.0

- New devel (Bioc 3.15)

# scp 1.4

- New stable release (Bioc 3.14)

# scp 1.3

## scp 1.3.3

- docs: included `QFeatures` plot in the vignette
- docs: created a vignette about advanced usage of `scp`

## scp 1.3.2

- feat: `computeSCR` now allows for user supplied function that will
  summarize the values from multiple samples and multiple carrier.
- docs: used more standard variable names in scp vignette.
- docs: created a `QFeatures` recap vignette

## scp 1.3.1

- refactor: deprecated `rowDataToDF`. This function is now replaced by
  `QFeatures::rbindRowData`.

## scp 1.3.0

- New devel (Bioc 3.14)

# scp 1.2

- New stable release (Bioc 3.13)

# scp 1.1

## scp 1.1.6

- feature: `readSCP` now allows for a `suffix` argument to better
  customize the sample names. <2021-03-17>

## scp 1.1.5

- deprecation: thanks to the new normalization method in `medianCVperCell`,
  'computeMedianCV_SCoPE2' is now deprecated and should no longer be
  used. <2021-02-19>
- feat: added a new normalization method to `medianCVperCell`. The
  `SCoPE2` normalization method can now reproduce the results from
  SCoPE2. <2021-02-19>
- docs: improved vignette <2021-02-16>
- feat: added a `rowDataName` argument to `computeSCR`
  <2021-02-08>

## scp 1.1.4

- fix: removed bug in vignette header
  <2021-02-06>
- data: update the example data with the latest release of SCoPE2
  <2021-02-06>
- feat: added `removeEmptyCol` argument in `readSCP` to automatically
  remove columns that contain only NA's
  <2021-02-06>

## scp 1.1.3

- docs: improved the manual page for `pep2qvalue` and the
  corresponding section in the vignette.
  <2021-01-23>
- refactor: reimplemented the `computeFDR` to catch up with the new
  release of SCoPE2. `computeFDR` was renamed to `pep2qvalue`. This
  is more in line with the theory. Also adapted the unit tests.
  <2021-01-23>

## scp 1.1.2

- docs: improved the description of the `scp` data structure in the
  vignette
  <2021-01-05>
- refactor: renamed the `groupCol` to `groupBy` and `pepCol` to `PEP`
  in `computeFDR`.
  <2020-12-08>
- refactor: renamed `computeMedianCV` to `computeMedianCV_SCoPE2`
  and deprecated the function. The function will be preserved for
  backward compatibility with the replication of the SCoPE2 analysis
  (Specht et al. 2020). Instead, a new function is implemented and
  called `medianCVperCell`. See
  [issue#7](https://github.com/UCLouvain-CBIO/scp/issues/7) for more
  information
  <2020-12-07>

## scp 1.1.1

- Fix news file

## scp 1.1.0

- New devel (Bioc 3.13)

# scp 1.0

## scp 1.0.0

- New stable release (Bioc 3.12)

# scp 0.99

## scp 0.99.4

- Update installation instructions <2020-10-14 Wed>

## scp 0.99.3

- fix: solved 'invalid subsetting' issue
  <2020-10-14 Wed>
- Adapted the vignette to remove warnings and fix missing PCA plot.
  <2020-10-14 Wed>
- `README.md`: extended the installation guide, providing both a
  stable and a devel installation. <2020-10-13 Tue>
- Removed the `LazyLoad` from the `DESCRIPTION` file and adapted the
  data loading (eg `data(scp1)` to `data("scp1")`)
  <2020-10-13 Tue>
- Documentation: added data collection description for the 3 example
  datasets
  <2020-10-13 Tue>

## scp 0.99.2

- fix: `computeFDR` can handle missing values (see issue #12)
  <2020-10-02 Fri>

## scp 0.99.1

- Maintainer subscribed to bioc-devel mailing list
- Removed `infIsNA`, the implementation was moved to the `QFeatures`
packages

## scp 0.99.0

- Bioconductor submission