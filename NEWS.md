# scp 0.1

## scp 0.1.3

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

## scp 0.1.2

- deprecated: removed the `transferColDataToAssay`. You should better 
  use the `getWithColData` function from `MultiAssayExperiment`.
  <2020-12-01>

## scp 0.1.1

- refactor: renamed all `obj` arguments by `object`
  <2020-12-01>
- feat: new function `normalizeSCP` that allow normalizing an assay in a 
  `QFeatures` object that contains `SingleCellExperiment` objects
  <2020-11-30>

## scp 0.1.0

- `scp` package acceptance on Bioconductor!!
  <2020-10-15 Thu>

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
