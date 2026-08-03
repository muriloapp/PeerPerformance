## ---------------------------------------------------------------------------
## Extended Monte-Carlo validation of the PeerPerformance estimators
##
## The companion script validation.R checks the estimators under i.i.d. Gaussian
## returns with independent funds. That is the setting in which the luck
## correction is easiest; it is not the setting the package exists for. This
## script covers the two questions that regime cannot answer:
##
##   (D) What does the estimator return under an exact null at the design of a
##       realistic application (N = 100 funds, T = 60 months)? Without this the
##       empirical mean pi0 of about 0.78 on hfdata has no benchmark: a reader
##       cannot tell whether 22% non-ties is more than the estimator reports
##       when nothing is there.
##
##   (E) How does it behave when the assumptions are violated in the ways fund
##       returns actually violate them -- fat tails, serial correlation, and
##       cross-sectional correlation (hfdata has a mean pairwise correlation of
##       about 0.25)?
##
##   (F) Does the modified Sharpe test keep its size under non-normality and
##       autocorrelation, and does the studentized circular bootstrap
##       (type = 2, bBoot = 0) repair it where the asymptotic test fails? The
##       article recommends the bootstrap for exactly these cases, so the
##       recommendation should rest on evidence.
##
## This script is slower than validation.R (tens of minutes). Run with:
##   source(system.file("scripts", "validation-extended.R",
##                      package = "PeerPerformance"))
## ---------------------------------------------------------------------------

library("PeerPerformance")
set.seed(1234)

## fastAdjust only changes the pi0 bias correction by a few 1e-5 (see NEWS) and
## makes the data-driven lambda affordable at these sample sizes.
ctr <- list(nCore = 1, fastAdjust = TRUE)

RS <- 50L    # replications for the screening scenarios
RT <- 400L   # replications for the test-size scenarios

## helper: mean of the triple over funds, for one screening
ratios <- function(X) {
  sc <- alphaScreening(X, control = ctr)
  c(pi0 = mean(sc$pizero, na.rm = TRUE),
    pip = mean(sc$pipos,  na.rm = TRUE),
    pim = mean(sc$pineg,  na.rm = TRUE))
}

## ===========================================================================
## (D) Null calibration across the design grid
## ===========================================================================
cat("\n(D) Exact null: what the estimator returns when nothing is there\n")
cat("    (true pi0 = 1; the non-tie mass pi+ + pi- is the false-discovery floor)\n\n")
cat(sprintf("    %5s %5s   %14s %14s\n", "N", "T", "pi0 (s.e.)", "pi+ + pi- (s.e.)"))
gridD <- expand.grid(N = c(30L, 100L), T = c(60L, 120L))
outD <- vector("list", nrow(gridD))
for (g in seq_len(nrow(gridD))) {
  N <- gridD$N[g]; TT <- gridD$T[g]
  m <- matrix(NA_real_, RS, 2)
  for (r in seq_len(RS)) {
    X <- matrix(stats::rnorm(TT * N, 0.005, 0.04), TT, N)
    z <- ratios(X)
    m[r, ] <- c(z["pi0"], z["pip"] + z["pim"])
  }
  outD[[g]] <- c(N = N, T = TT, colMeans(m),
                 apply(m, 2, stats::sd) / sqrt(RS))
  cat(sprintf("    %5d %5d   %.3f (%.3f)   %.3f (%.3f)\n", N, TT,
              mean(m[, 1]), stats::sd(m[, 1])/sqrt(RS),
              mean(m[, 2]), stats::sd(m[, 2])/sqrt(RS)))
}
cat("\n    For comparison, the empirical mean pi0 on hfdata (N = 100, T = 60)\n")
data("hfdata")
set.seed(1234)
emp <- mean(alphaScreening(hfdata, control = ctr)$pizero, na.rm = TRUE)
cat(sprintf("    is %.3f, i.e. a non-tie mass of %.3f.\n", emp, 1 - emp))

## ===========================================================================
## (E) Departures from the i.i.d. Gaussian null, at the empirical design
## ===========================================================================
cat("\n(E) Exact null at N = 100, T = 60 under violated assumptions\n\n")
N <- 100L; TT <- 60L

genGauss <- function() matrix(stats::rnorm(TT * N, 0.005, 0.04), TT, N)
genT5    <- function() matrix(0.005 + 0.04 * stats::rt(TT * N, df = 5)/sqrt(5/3),
                              TT, N)
genAR1   <- function(rho = 0.3) {
  e <- matrix(stats::rnorm(TT * N, 0, 0.04 * sqrt(1 - rho^2)), TT, N)
  X <- matrix(NA_real_, TT, N)
  X[1, ] <- stats::rnorm(N, 0, 0.04)
  for (t in 2:TT) X[t, ] <- rho * X[t - 1, ] + e[t, ]
  X + 0.005
}
genFac   <- function(rho = 0.25) {           # common factor => cross-sectional rho
  f <- stats::rnorm(TT, 0, 0.04 * sqrt(rho))
  X <- matrix(stats::rnorm(TT * N, 0, 0.04 * sqrt(1 - rho)), TT, N)
  0.005 + X + f
}
## AR(1) idiosyncratic returns plus an AR(1) common factor. Calibrated to
## reproduce BOTH the serial and the cross-sectional dependence of hfdata, so
## that the resulting floor is the right benchmark for the empirical section.
genBoth  <- function(rho = 0.20, rc = 0.25) {
  ar1 <- function(n, s) {
    e <- stats::rnorm(n, 0, s * sqrt(1 - rho^2))
    x <- numeric(n); x[1] <- stats::rnorm(1, 0, s)
    for (t in 2:n) x[t] <- rho * x[t - 1] + e[t]
    x
  }
  f <- ar1(TT, 0.04 * sqrt(rc))
  X <- matrix(0, TT, N)
  for (j in seq_len(N)) X[, j] <- ar1(TT, 0.04 * sqrt(1 - rc))
  0.005 + X + f
}
gens <- list("i.i.d. Gaussian (reference)" = genGauss,
             "t(5) innovations"            = genT5,
             "AR(1), rho = 0.2"            = function() genAR1(0.20),
             "AR(1), rho = 0.3"            = genAR1,
             "common factor, rho = 0.25"   = genFac,
             "AR(1) 0.2 + factor 0.25"     = genBoth)
cat(sprintf("    %-28s %14s %14s\n", "return process", "pi0 (s.e.)", "pi+ + pi- (s.e.)"))
for (nm in names(gens)) {
  m <- matrix(NA_real_, RS, 2)
  for (r in seq_len(RS)) {
    z <- ratios(gens[[nm]]())
    m[r, ] <- c(z["pi0"], z["pip"] + z["pim"])
  }
  cat(sprintf("    %-28s %.3f (%.3f)   %.3f (%.3f)\n", nm,
              mean(m[, 1]), stats::sd(m[, 1])/sqrt(RS),
              mean(m[, 2]), stats::sd(m[, 2])/sqrt(RS)))
}

## ===========================================================================
## (E2) The maximum pi+ over the cross-section, under the calibrated null
##
## Picking the best fund out of N is an extreme-value operation: max_i pi+_i
## must be judged against the distribution of that maximum under the null,
## not against zero. Without this, the single most striking number in any
## applied screening has no benchmark at all.
## ===========================================================================
cat("\n(E2) Largest pi+ across the cross-section, calibrated null (AR(1) 0.2 +",
    "factor 0.25)\n\n")
mx <- numeric(RS)
for (r in seq_len(RS)) {
  sc <- alphaScreening(genBoth(), control = ctr)
  mx[r] <- max(sc$pipos, na.rm = TRUE)
}
data("hfdata")
set.seed(1234)
empmax <- max(alphaScreening(hfdata, control = ctr)$pipos, na.rm = TRUE)
cat(sprintf("    null max pi+ : mean %.3f (s.d. %.3f), 95th pct %.3f\n",
            mean(mx), stats::sd(mx), stats::quantile(mx, 0.95)))
cat(sprintf("    hfdata max pi+ : %.3f -> %.0f%% of null replications exceed it\n",
            empmax, 100 * mean(mx >= empmax)))

## Dependence actually present in hfdata, for reference
ac1 <- apply(hfdata, 2, function(z) {
  z <- z[is.finite(z)]; stats::cor(z[-1], z[-length(z)])
})
cc <- stats::cor(hfdata, use = "pairwise.complete.obs")
cat(sprintf("    hfdata dependence: mean lag-1 AC %.3f (%.0f%% above 0.2),",
            mean(ac1, na.rm = TRUE), 100 * mean(ac1 > 0.2, na.rm = TRUE)))
cat(sprintf(" mean pairwise cor %.3f\n", mean(cc[upper.tri(cc)], na.rm = TRUE)))

## ===========================================================================
## (F) Size of the modified Sharpe test: asymptotic vs studentized bootstrap
## ===========================================================================
cat("\n(F) Empirical size of the modified Sharpe equality test (nominal 5%)\n")
cat("    asymptotic (type = 1, the default) vs studentized circular bootstrap\n")
cat("    (type = 2, bBoot = 0, i.e. data-driven block length)\n\n")
binse <- function(p, R) sqrt(p * (1 - p)/R)
TT <- 120L
pairGauss <- function() cbind(stats::rnorm(TT, 0.006, 0.04),
                              stats::rnorm(TT, 0.006, 0.04))
pairT5    <- function() cbind(0.006 + 0.04 * stats::rt(TT, 5)/sqrt(5/3),
                              0.006 + 0.04 * stats::rt(TT, 5)/sqrt(5/3))
pairAR1   <- function(rho = 0.3) {
  sim <- function() as.numeric(stats::arima.sim(list(ar = rho), TT,
                                                sd = 0.04 * sqrt(1 - rho^2))) + 0.006
  cbind(sim(), sim())
}
pairs <- list("i.i.d. Gaussian" = pairGauss,
              "t(5) innovations" = pairT5,
              "AR(1), rho = 0.3" = pairAR1)
cat(sprintf("    %-20s %18s %18s\n", "return process",
            "asymptotic (s.e.)", "bootstrap (s.e.)"))
for (nm in names(pairs)) {
  rej1 <- rej2 <- 0L
  for (r in seq_len(RT)) {
    xy <- pairs[[nm]]()
    p1 <- msharpeTesting(xy[, 1], xy[, 2], level = 0.90)$pval
    rej1 <- rej1 + isTRUE(p1 < 0.05)
    p2 <- msharpeTesting(xy[, 1], xy[, 2], level = 0.90,
                         control = list(type = 2, bBoot = 0, nBoot = 199))$pval
    rej2 <- rej2 + isTRUE(p2 < 0.05)
  }
  cat(sprintf("    %-20s %.3f (%.3f)      %.3f (%.3f)\n", nm,
              rej1/RT, binse(rej1/RT, RT), rej2/RT, binse(rej2/RT, RT)))
}

cat("\nsettings: replications", RS, "(screening) and", RT,
    "(test size); seed 1234; fastAdjust = TRUE\n")
