## Set of R function for targetPeerPerformance

# @name .targetPeerPerformance
# @description See targetPeerPerformance
.targetPeerPerformance <- function(X,
                                   funds,
                                   method = c("alpha", "sharpe", "msharpe"),
                                   factors = NULL,
                                   level = 0.9,
                                   na.neg = TRUE,
                                   control = list()) {
  
  method <- match.arg(method)
  ctr <- processControl(control)
  
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  
  if (is.null(colnames(X))) colnames(X) <- paste0("Fund ", seq_len(N))
  fund_names <- colnames(X)
  
  # normalize funds -> indices (preserve user order)
  idx <- if (is.numeric(funds)) {
    as.integer(funds)
  } else {
    m <- match(funds, fund_names)
    if (anyNA(m)) stop("Unknown fund(s): ", paste(funds[is.na(m)], collapse = ", "))
    m
  }
  
  if (any(idx < 1L | idx > N)) stop("Some 'funds' indices are out of range.")
  if (anyDuplicated(idx)) stop("Duplicate entries in 'funds' are not allowed.")
  m <- length(idx)
  
  # basic alpha-specific checks (ONLY if factors provided, as in your original code)
  if (method == "alpha" && !is.null(factors)) {
    factors <- as.matrix(factors)
    if (nrow(factors) != T) stop("'factors' must have the same number of rows as 'X'.")
  }
  
  # permute so focal funds are first (so Screeningi compares them vs ALL peers to the right)
  perm <- c(idx, setdiff(seq_len(N), idx))
  Xp <- X[, perm, drop = FALSE]
  
  # allocate pairwise objects (mostly NA; we compute only what we need)
  pval  <- matrix(NA_real_, N, N)
  ddiff <- matrix(NA_real_, N, N)
  tstat <- matrix(NA_real_, N, N)
  
  if (N >= 2L) {
    
    # precompute bootstrap ids if needed
    bsids <- NULL
    if (method %in% c("sharpe", "msharpe")) {
      bsids <- bootIndices(T, ctr$nBoot, ctr$bBoot)
    }
    
    # only need to run i = 1..m (focals) as long as i < N
    ids <- seq_len(m)
    ids <- ids[ids < N]
    
    if (length(ids) > 0L) {
      for (i in ids) {
        
        if (method == "alpha") {
          
          outi <- alphaScreeningi(
            i = i, rdata = Xp, factors = factors, T = T, N = N,
            hac = ctr$hac, screen_beta = FALSE
          )
          
          pvi <- outi$pvali[1, ]
          dvi <- outi$dalphai[1, ]
          tvi <- outi$tstati[1, ]
          
          pval[i,  i:N] <- pval[i:N, i] <- pvi[i:N]
          ddiff[i, i:N] <- dvi[i:N]
          ddiff[i:N, i] <- -dvi[i:N]
          tstat[i, i:N] <- tvi[i:N]
          tstat[i:N, i] <- -tvi[i:N]
          
        } else if (method == "sharpe") {
          
          outi <- sharpeScreeningi(
            i = i, rdata = Xp, T = T, N = N,
            nBoot = ctr$nBoot, bsids = bsids, minObs = ctr$minObs,
            type = ctr$type, hac = ctr$hac, b = ctr$bBoot,
            ttype = ctr$ttype, pBoot = ctr$pBoot
          )
          
          pvi <- outi$pvali
          dvi <- outi$dsharpei
          tvi <- outi$tstati
          
          pval[i,  i:N] <- pval[i:N, i] <- pvi[i:N]
          ddiff[i, i:N] <- dvi[i:N]
          ddiff[i:N, i] <- -dvi[i:N]
          tstat[i, i:N] <- tvi[i:N]
          tstat[i:N, i] <- -tvi[i:N]
          
        } else { # msharpe
          
          outi <- msharpeScreeningi(
            i = i, rdata = Xp, level = level, T = T, N = N,
            nBoot = ctr$nBoot, bsids = bsids, minObs = ctr$minObs,
            na.neg = na.neg, type = ctr$type, hac = ctr$hac, b = ctr$bBoot,
            ttype = ctr$ttype, pBoot = ctr$pBoot
          )
          
          pvi <- outi$pvali
          dvi <- outi$dmsharpei
          tvi <- outi$tstati
          
          pval[i,  i:N] <- pval[i:N, i] <- pvi[i:N]
          ddiff[i, i:N] <- dvi[i:N]
          ddiff[i:N, i] <- -dvi[i:N]
          tstat[i, i:N] <- tvi[i:N]
          tstat[i:N, i] <- -tvi[i:N]
        }
      }
    }
  }
  
  # focal positions in permuted matrix are 1..m
  sel <- seq_len(m)
  nm  <- colnames(Xp)[sel]
  
  # compute pi ONLY for focals
  pval_sel  <- pval[sel, , drop = FALSE]
  ddiff_sel <- ddiff[sel, , drop = FALSE]
  tstat_sel <- tstat[sel, , drop = FALSE]
  
  # lambda handling: allow scalar/NULL, or length-N vectors (subset to focals)
  lambda_sel <- ctr$lambda
  if (!is.null(lambda_sel) && length(lambda_sel) > 1L) {
    if (length(lambda_sel) != N) stop("'lambda' must be length 1, NULL, or length ncol(X).")
    lambda_sel <- lambda_sel[perm][sel]
  }
  
  pi <- if (N >= 2L) {
    computePi(
      pval   = pval_sel,
      dalpha = ddiff_sel,
      tstat  = tstat_sel,
      lambda = lambda_sel,
      nBoot  = ctr$nBoot
    )
  } else {
    list(
      pizero = rep(NA_real_, m),
      pipos  = rep(NA_real_, m),
      pineg  = rep(NA_real_, m)
    )
  }
  
  # info ONLY for focals
  Xf <- Xp[, sel, drop = FALSE]
  info <- switch(
    method,
    alpha   = infoFund(Xf, factors = factors, screen_beta = FALSE),
    sharpe  = infoFund(Xf),
    msharpe = infoFund(Xf, level = level, na.neg = na.neg)
  )
  
  metric_vec <- switch(
    method,
    alpha   = info$alpha,
    sharpe  = info$sharpe,
    msharpe = info$msharpe
  )
  
  # ---- build focal-vs-all matrices in ORIGINAL (input) column order ----
  inv_perm <- match(seq_len(N), perm)  # position of original cols inside permuted cols
  
  pval_out  <- pval_sel[,  inv_perm, drop = FALSE]
  ddiff_out <- ddiff_sel[, inv_perm, drop = FALSE]
  tstat_out <- tstat_sel[, inv_perm, drop = FALSE]
  
  rownames(pval_out)  <- nm
  rownames(ddiff_out) <- nm
  rownames(tstat_out) <- nm
  
  colnames(pval_out)  <- fund_names
  colnames(ddiff_out) <- fund_names
  colnames(tstat_out) <- fund_names
  
  # number of available peers for each focal (exclude self)
  self_non_na <- !is.na(pval_out[cbind(seq_len(m), idx)])  # idx are ORIGINAL fund indices
  npeer <- rowSums(!is.na(pval_out)) - as.integer(self_non_na)
  
  # difference object name depends on method
  dname <- switch(
    method,
    alpha   = "dalpha",
    sharpe  = "dsharpe",
    msharpe = "dmsharpe"
  )
  
  # lambda to return (what was actually used for the focals)
  lambda_out <- lambda_sel
  if (!is.null(lambda_out) && length(lambda_out) > 1L) names(lambda_out) <- nm
  
  # base output (focal-only vectors + focal-vs-all matrices)
  res <- list(
    n      = setNames(info$nObs, nm),
    npeer  = setNames(npeer, nm),
    
    pval   = pval_out,
    tstat  = tstat_out,
    lambda = lambda_out,
    
    pizero = setNames(pi$pizero, nm),
    pipos  = setNames(pi$pipos, nm),
    pineg  = setNames(pi$pineg, nm)
  )
  
  # add metric under its true name (alpha/sharpe/msharpe)
  res[[method]] <- setNames(metric_vec, nm)
  
  # add corresponding difference matrix (dalpha/dsharpe/dmsharpe)
  res[[dname]] <- ddiff_out
  
  return(res)
}


#' @name targetPeerPerformance
#' @title Targeted peer-performance screening for selected funds
#' @description
#' Performs peer-performance screening for a user-selected subset of funds (the
#' "focal" funds) against all available peers, returning detailed outputs for
#' those focal funds only.
#'
#' @details
#' Compared to \code{\link{alphaScreening}}, \code{\link{sharpeScreening}}, and
#' \code{\link{msharpeScreening}} (which compute all pairwise tests), this function
#' computes pairwise tests only for each focal fund versus all peers.
#'
#' The matrices \code{pval}, \code{tstat}, and \code{dalpha}/\code{dsharpe}/\code{dmsharpe}
#' are returned with rows corresponding to focal funds (in the order provided by
#' \code{funds}) and columns corresponding to all funds in \code{X} (original input order).
#'
#' @param X Matrix \eqn{(T \times N)}{(TxN)} of \eqn{T} returns for \eqn{N} funds.
#' \code{NA} values are allowed.
#' @param funds Integer indices or character column names identifying the focal funds.
#' The output preserves the order provided in \code{funds}.
#' @param method Screening method: \code{"alpha"}, \code{"sharpe"}, or \code{"msharpe"}.
#' @param factors Optional matrix \eqn{(T \times K)}{(TxK)} of factor returns (used only when
#' \code{method="alpha"}). If supplied, it must have the same number of rows as \code{X}.
#' It is passed through to the underlying \code{alphaScreeningi()} / \code{infoFund()} calls.
#' @param level Modified Value-at-Risk level (only used when \code{method="msharpe"}).
#' @param na.neg Logical; if \code{TRUE}, returns \code{NA} when a negative modified
#' Value-at-Risk is obtained (only used when \code{method="msharpe"}).
#' @param control Control parameters passed to \code{processControl()} (see
#' \code{\link{alphaScreening}}, \code{\link{sharpeScreening}}, \code{\link{msharpeScreening}}
#' for typical fields, including \code{nCore}, \code{type}, \code{hac}, \code{nBoot}, \code{bBoot},
#' \code{ttype}, \code{pBoot}, and \code{lambda}).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{n}: named numeric vector (length \code{m = length(funds)}), number of non-\code{NA} observations for each focal fund.
#'   \item \code{npeer}: named integer vector (length \code{m}), number of available peers for each focal fund (excluding self; based on non-\code{NA} p-values).
#'   \item \code{pval}: numeric matrix \eqn{(m \times N)}{(m x N)} of p-values for focal-vs-peer tests.
#'   \item \code{tstat}: numeric matrix \eqn{(m \times N)}{(m x N)} of test statistics for focal-vs-peer tests.
#'   \item \code{lambda}: \code{NULL}, a scalar, or a named numeric vector (length \code{m}) containing the \code{lambda} values used for the focal funds in \code{computePi()}.
#'   \item \code{pizero}: named numeric vector (length \code{m}), probability of equal performance.
#'   \item \code{pipos}: named numeric vector (length \code{m}), probability of outperformance.
#'   \item \code{pineg}: named numeric vector (length \code{m}), probability of underperformance.
#'   \item one metric vector: \code{alpha} or \code{sharpe} or \code{msharpe} (named numeric vector, length \code{m}).
#'   \item one difference matrix: \code{dalpha} or \code{dsharpe} or \code{dmsharpe} (numeric matrix \eqn{(m \times N)}{(m x N)}).
#' }
#'
#' @seealso \code{\link{alphaScreening}}, \code{\link{sharpeScreening}}, \code{\link{msharpeScreening}}.
#'
#' @examples
#' data("hfdata")
#' rets <- hfdata[, 1:10]
#'
#' ## focal Sharpe screening (funds by index)
#' out <- targetPeerPerformance(rets, funds = c(2, 5, 7), method = "sharpe",
#'                              control = list(nCore = 1, type = 2))
#' out$sharpe
#' out$dsharpe
#' out$pval
#' out$tstat
#'
#' @export
#' @importFrom compiler cmpfun
targetPeerPerformance <- compiler::cmpfun(.targetPeerPerformance)
