#' Summarize SCREENING objects
#'
#' S3 method for \code{summary()} for objects of class \code{"SCREENING"}.
#'
#' @param object A SCREENING object.
#' @param coef If the selected measure is matrix/array-like, which row/coef to use (default 1).
#' @param top How many best/worst funds to show.
#' @param p_level Significance level for win/loss counts (requires pval+tstat in object).
#' @return An object of class \code{"summary.SCREENING"}.
#' @export
#' @method summary SCREENING
summary.SCREENING <- function(x,
                              coef = 1L,
                              top = 5L,
                              p_level = 0.05
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
  
  est <- .pp_get_measure(x, measure = measure, coef = coef)
  
  funds <- names(x$n)
  if (is.null(funds) || length(funds) != length(est)) {
    funds <- paste0("Fund ", seq_along(est))
  }
  
  # helper to align vectors by fund names when available
  pick_vec <- function(v) {
    if (is.null(v)) return(rep(NA_real_, length(funds)))
    if (!is.null(names(v))) return(as.numeric(v[funds]))
    as.numeric(v)
  }
  
  tab <- data.frame(
    fund     = funds,
    n        = pick_vec(x$n),
    npeer    = pick_vec(x$npeer),
    estimate = as.numeric(est),
    lambda   = pick_vec(x$lambda),
    pizero   = pick_vec(x$pizero),  # π^0
    pipos    = pick_vec(x$pipos),   # π^+
    pineg    = pick_vec(x$pineg),   # π^-
    stringsAsFactors = FALSE
  )
  
  # win/loss counts from pairwise tests (needs pval+tstat)
  if (!is.null(x$pval) && !is.null(x$tstat) &&
      is.matrix(x$pval) && is.matrix(x$tstat) &&
      nrow(x$pval) == length(funds) && ncol(x$pval) == length(funds)) {
    
    pmat <- x$pval
    tmat <- x$tstat
    
    tab$wins_p   <- rowSums((pmat < p_level) & (tmat > 0), na.rm = TRUE)
    tab$losses_p <- rowSums((pmat < p_level) & (tmat < 0), na.rm = TRUE)
    tab$net_p    <- tab$wins_p - tab$losses_p
  }
  
  # ranks 
  tab$rank <- rank(-tab$estimate, ties.method = "first")
  tab <- tab[order(tab$rank), , drop = FALSE]
  
  # quick stats
  est_stats <- c(
    min  = min(tab$estimate, na.rm = TRUE),
    q25  = unname(stats::quantile(tab$estimate, 0.25, na.rm = TRUE)),
    med  = stats::median(tab$estimate, na.rm = TRUE),
    mean = mean(tab$estimate, na.rm = TRUE),
    q75  = unname(stats::quantile(tab$estimate, 0.75, na.rm = TRUE)),
    max  = max(tab$estimate, na.rm = TRUE)
  )
  
  top <- as.integer(top)
  top <- max(1L, min(top, nrow(tab)))
  
  # Top by π+ and π- 
  top_pi_plus  <- NULL
  top_pi_minus <- NULL
  if (all(c("pipos", "pineg") %in% names(tab))) {
    ord_plus  <- order(-tab$pipos,  -tab$estimate)
    ord_minus <- order(-tab$pineg,  -tab$estimate)
    top_pi_plus  <- utils::head(tab[ord_plus,  , drop = FALSE], top)
    top_pi_minus <- utils::head(tab[ord_minus, , drop = FALSE], top)
  }
  
  res <- list(
    call         = match.call(),
    measure      = measure,
    coef         = coef,
    p_level      = p_level,
    n_funds      = nrow(tab),
    stats        = est_stats,
    table        = tab,
    top          = utils::head(tab, top),
    bottom       = utils::tail(tab, top),
    top_pi_plus  = top_pi_plus,
    top_pi_minus = top_pi_minus
  )
  class(res) <- "summary.SCREENING"
  res
}

#' @export
#' @method print summary.SCREENING
print.summary.SCREENING <- function(x, ...) {
  cat("SCREENING summary\n")
  cat("  Measure:", x$measure, "\n")
  cat("  Funds:", x$n_funds, "\n\n")
  
  # formatting helpers 
  fmt5 <- function(z) formatC(as.numeric(z), format = "f", digits = 5)
  round_df5 <- function(df) {
    if (is.null(df)) return(df)
    cols <- intersect(c("estimate", "pizero", "pipos", "pineg"), names(df))
    for (cc in cols) df[[cc]] <- fmt5(df[[cc]])
    df
  }
  
  s <- x$stats
  
  cat(sprintf("%s summary\n", x$measure))
  cat(sprintf("  Max   : %s\n\n", fmt5(s["max"])))
  cat(sprintf("  75%%   : %s\n", fmt5(s["q75"])))
  cat(sprintf("  Mean  : %s\n", fmt5(s["mean"])))
  cat(sprintf("  Median: %s\n", fmt5(s["med"])))
  cat(sprintf("  25%%   : %s\n", fmt5(s["q25"])))
  cat(sprintf("  Min   : %s\n", fmt5(s["min"])))
  cat("\n")
  
  
  if (!is.null(x$top_pi_plus)) {
    cat("Top by outperformance ratio (π+):\n")
    print(round_df5(x$top_pi_plus), row.names = FALSE, right = TRUE)
    cat("\n")
  }
  
  if (!is.null(x$top_pi_minus)) {
    cat("Top by underperformance ratio (π-):\n")
    print(round_df5(x$top_pi_minus), row.names = FALSE, right = TRUE)
    cat("\n")
  }
  
  cat(sprintf("Top funds (by %s estimate)\n", x$measure))
  print(round_df5(x$top), row.names = FALSE, right = TRUE)
  
  cat(sprintf("\nBottom funds (by %s estimate)\n", x$measure))
  print(round_df5(x$bottom), row.names = FALSE, right = TRUE)
  
  invisible(x)
}


