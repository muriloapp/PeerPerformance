# Screening plotting helpers + S3 plot method for class "SCREENING"


# helper: pick the sorted measure from screening output 
#' @keywords internal
#' @noRd
.pp_get_measure <- function(x, measure = NULL, coef = 1L) {
  stopifnot(!is.null(measure), measure %in% names(x))
  
  m <- x[[measure]]
  
  # If screen_beta=TRUE in alphaScreening, alpha can become an array-like object
  # where coef=1 corresponds to alpha
  if (!is.null(dim(m)) && length(dim(m)) >= 2) {
    m <- m[coef, ]
  }
  as.numeric(m)
}

# helper: bucketization + bucket averages 
#' @keywords internal
#' @noRd
bucketize_screening <- function(x, nbucket = 50L, measure = NULL, annualize = 1,
                                decreasing = TRUE, coef = 1L) {
  
  m <- .pp_get_measure(x, measure = measure, coef = coef)
  
  # required peer ratios
  stopifnot(all(c("pipos","pizero","pineg") %in% names(x)))
  
  pi_plus <- x$pipos
  if (!is.null(dim(pi_plus)) && length(dim(pi_plus)) >= 2) pi_plus <- pi_plus[coef, ]
  pi_plus <- as.numeric(pi_plus)
  
  pi_zero <- x$pizero
  if (!is.null(dim(pi_zero)) && length(dim(pi_zero)) >= 2) pi_zero <- pi_zero[coef, ]
  pi_zero <- as.numeric(pi_zero)
  
  pi_minus <- x$pineg
  if (!is.null(dim(pi_minus)) && length(dim(pi_minus)) >= 2) pi_minus <- pi_minus[coef, ]
  pi_minus <- as.numeric(pi_minus)
  
  ok <- is.finite(m) & is.finite(pi_plus) & is.finite(pi_zero) & is.finite(pi_minus)
  m <- m[ok]; pi_plus <- pi_plus[ok]; pi_zero <- pi_zero[ok]; pi_minus <- pi_minus[ok]
  
  N <- length(m)
  nbucket <- as.integer(min(nbucket, N))
  ord <- order(m, decreasing = decreasing)
  
  # nearly equal-sized buckets
  grp <- cut(seq_len(N), breaks = nbucket, labels = FALSE, include.lowest = TRUE)
  
  out <- lapply(seq_len(nbucket), function(b) {
    idx <- ord[grp == b]
    a   <- mean(m[idx], na.rm = TRUE) * annualize
    p1  <- mean(pi_plus[idx],  na.rm = TRUE)
    p0  <- mean(pi_zero[idx],  na.rm = TRUE)
    p_1 <- mean(pi_minus[idx], na.rm = TRUE)
    
    # numerical safety
    v <- pmax(0, pmin(1, c(p1, p0, p_1)))
    s <- sum(v)
    if (s > 0) v <- v / s
    
    data.frame(bucket = b, metric = a, pipos = v[1], pizero = v[2], pineg = v[3])
  })
  
  do.call(rbind, out)
}


#' Null-coalescing operator
#' Returns `a` if it is not NULL, otherwise returns `b`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


#' @export
#' @noRd
# plot method 
plot.pp_screening <- function(x, nbucket = 50L,
                              annualize = attr(x, "annualize") %||% 1,
                              coef = 1L,
                              show_diag = TRUE
                              ) {
  
  # infer the measure directly from x
  if ("alpha" %in% names(x)) {
    measure <- "alpha"
  } else if ("sharpe" %in% names(x)) {
    measure <- "sharpe"
  } else if ("msharpe" %in% names(x)) {
    measure <- "msharpe"
  } else {
    stop("No recognized measure found in x (expected one of: alpha, sharpe, msharpe).")
  }
  
  B <- bucketize_screening(x, nbucket = nbucket, measure = measure,
                           annualize = annualize, coef = coef)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # Left: bucket mean metric vs bucket rank 
  y <- B$bucket  
  plot(B$metric, y,
       pch = 16, cex = 0.7,
       xlab = "", ylab = "",
       yaxt = "n", main = measure,
       )
  axis(2, at = pretty(y))
  grid()
  
  lines(B$metric, y)
  
  # Right: stacked bars per bucket 
  nb <- nrow(B)
  
  plot(NA,
       xlim = c(0, 1), ylim = c(0, 1),
       xaxs = "i", yaxs = "i",    
       xlab = "", ylab = "",
       axes = FALSE, main = "")
  
  # fill each bucket as a horizontal strip in [0,1]
  for (b in seq_len(nb)) {
    y0 <- (b - 1) / nb
    y1 <- b / nb
    
    x1 <- sort(B$pipos, decreasing = TRUE)[b]
    x2 <- x1 + sort(B$pizero, decreasing = TRUE)[b]
    
    rect(0,  y0, x1, y1, col = "black",  border = NA)
    rect(x1, y0, x2, y1, col = "grey80", border = NA)
    rect(x2, y0, 1,  y1, col = "grey40", border = NA)
  }
  
  axis(1, at = seq(0, 1, by = 0.25), labels = paste0(seq(0, 100, 25), "%"))
  mtext(expression(hat(pi)^"+" ~ "/" ~ hat(pi)^0 ~ "/" ~ hat(pi)^"-"),
        side = 3, line = 0.5)
  
  box(lwd = 2)
  
  # diagonals: parallel 45-degree guides in the unit square
  if (show_diag) {
    for (cst in c(2/3, 1, 4/3)) { 
      abline(a = cst, b = -1, lty = 2)
    }
  }
  
  invisible(B)
}

# S3 method: makes plot(name) work when class(name) includes "SCREENING" 
#' Plot SCREENING objects
#'
#' S3 method for \code{plot()} for objects of class \code{"SCREENING"}.
#'
#' @param x A SCREENING object (returned by alphaScreening/sharpeScreening/msharpeScreening).
#' @param ... Passed to \code{plot.pp_screening()}.
#' @export
plot.SCREENING <- function(x, ...) {
  plot.pp_screening(x, ...)
}
