## Submission summary

This is a bug-fix release, submitted shortly after 2.4.0. Apologies for the
quick turnaround: 2.4.0 introduced no new defect, but a pre-existing one was
found immediately after its acceptance and it silently biases a statistical
result, so it seemed better to correct it than to leave it in place.

`bootIndices()` generated the circular block bootstrap by drawing
`floor(T / bBoot)` blocks. When the block length did not divide the number of
concordant observations, the last `T %% bBoot` indices of every bootstrap
sample kept their initial value of zero, and because the callers remap indices
with `1 + ids %% T` each of those zeros selected the first observation. At
`T = 50` and `bBoot = 6` this over-sampled observation 1 by a factor of about
three and biased the bootstrap null distribution, with no error, warning or
`NA` to signal it. The routine now draws `ceiling(T / bBoot)` blocks and keeps
the first `T` indices.

The default is `bBoot = 1` (i.i.d. resampling), which was never affected, and
results are unchanged whenever `bBoot` divides `T`.

A second fix is included: a fund with no usable observation no longer aborts a
screening (`infoFund()` called `lm()` on every column, so an all-`NA` fund
stopped the run with "0 (non-NA) cases"); its summary statistics are now
returned as `NA`.

Both fixes have regression tests.

## Test environments

* Local: macOS 26.5 (aarch64-apple-darwin20), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 2 notes

```
* checking CRAN incoming feasibility ... NOTE
Days since last update: 1
```
Explained above.

```
* checking for future file timestamps ... NOTE
unable to verify current time
```
Local clock/network issue, not a package issue.

## Reverse dependencies

None on CRAN.
