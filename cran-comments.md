## Submission summary

Feature and robustness update of PeerPerformance (2.3.2 -> 2.4.0). It adds
cross-group and rolling-window screening, bootstrap confidence intervals for
the peer performance ratios, a factor-exposure heterogeneity measure, a set of
S3 methods (print / summary / plot / as.data.frame / confint), a vignette, and
a reproducible Monte-Carlo validation script. It also fixes a bootstrap
resampling bug on unbalanced panels and several robustness issues; see NEWS.md
for the full list.

## Test environments

* Local: macOS 26.5 (aarch64-apple-darwin20), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are specific to the local machine and are not package issues:

```
* checking top-level files ... NOTE
Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
```
pandoc is not on the PATH of the local R session.

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```
The local `tidy` is the version shipped with macOS.

Neither note is expected on the CRAN check machines.

## Reverse dependencies

None on CRAN.
