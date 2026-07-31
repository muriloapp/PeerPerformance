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

0 errors | 0 warnings | 0 notes

## Reverse dependencies

None on CRAN.
