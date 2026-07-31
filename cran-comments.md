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

NOT YET RE-RUN FOR 2.4.0 -- do before submitting:
* win-builder, R-release and R-devel
* R-hub: Linux and Windows

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 1 note

The note is specific to this machine, where pandoc is not on the PATH of the
R session; it disappears where pandoc is available:

```
* checking top-level files ... NOTE
Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
```

## Reverse dependencies

None on CRAN.
