# scp 1.1

## scp 1.1.5

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
