# Internal helpers -- a faithful R port of the shipped reread.js engine.
# None of these are exported. Semantics match the browser tool so that the
# deterministic components reproduce its scores exactly (the permutation
# baseline matches to Monte Carlo error).

# Coerce a data frame / matrix of responses to a numeric matrix.
.as_matrix <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  storage.mode(data) <- "double"
  if (is.null(dim(data)) || ncol(data) < 2L)
    stop("`data` must be a matrix or data frame of item responses (>= 2 columns).")
  data
}

# Proportion / other rescaling + median imputation of missing cells.
# Returns the rescaled matrix P and the per-item rescaled minimum minP.
.prepare <- function(mat, rescale = "proportion") {
  J <- ncol(mat)
  mx  <- apply(mat, 2, max, na.rm = TRUE)
  mn  <- apply(mat, 2, min, na.rm = TRUE)
  med <- apply(mat, 2, stats::median, na.rm = TRUE)
  mx[!is.finite(mx) | mx == 0] <- 1
  mn[!is.finite(mn)] <- 0
  med[!is.finite(med)] <- 0
  P <- mat
  for (j in seq_len(J)) {
    v <- P[, j]; v[!is.finite(v)] <- med[j]; P[, j] <- v
  }
  P <- switch(rescale,
    proportion = sweep(P, 2, mx, "/"),
    none       = P,
    minmax     = sweep(sweep(P, 2, mn, "-"), 2, pmax(mx - mn, 1e-9), "/"),
    zscore     = scale(P),
    stop("Unknown `rescale`: ", rescale))
  P <- as.matrix(P)
  # rescaled per-item min/max for reverse-coding: reverse = (maxP + minP) - v.
  # For "proportion" this is exactly (1 + min/max) - v/max, matching the tool.
  list(P = P, minP = apply(P, 2, min), maxP = apply(P, 2, max))
}

# Sample-level correlations of every item pair (upper triangle).
.correlation_pairs <- function(P) {
  R <- suppressWarnings(stats::cor(P))
  R[!is.finite(R)] <- 0
  ut <- which(upper.tri(R), arr.ind = TRUE)
  list(pa = ut[, 1], pb = ut[, 2], pr = R[upper.tri(R)], n_pairs = nrow(ut))
}

# |cor(u, v)| down each respondent's row over a set of coupled pairs, with the
# partner item reverse-coded when the sample correlation of the pair is negative.
# Swap-invariant per-respondent coherence over exchangeable coupled pairs: the
# intraclass correlation = Pearson on the symmetrized (a,b)+(b,a) data. Unlike a plain
# cor(A,B) it does not depend on the arbitrary left/right labelling of each pair.
.ind_cors <- function(P, minP, maxP, pa, pb, pr, idx) {
  A <- P[, pa[idx], drop = FALSE]
  B <- P[, pb[idx], drop = FALSE]
  neg <- pr[idx] < 0
  if (any(neg)) {
    bcols <- pb[idx][neg]
    B[, neg] <- sweep(-B[, neg, drop = FALSE], 2, (maxP[bcols] + minP[bcols]), "+")
  }
  k <- ncol(A)
  m <- (rowSums(A) + rowSums(B)) / (2 * k)
  num <- rowSums(A * B) / k - m * m
  den <- (rowSums(A * A) + rowSums(B * B)) / (2 * k) - m * m
  ok <- den > 1e-12
  out <- numeric(nrow(P))
  out[ok] <- abs(num[ok] / den[ok])
  out
}

# Robust within-dataset z (median / MAD), matching the JS upper-median and the
# mean/sd fallback when the MAD is ~0.
.robust_z <- function(x) {
  n <- length(x)
  med <- sort(x)[floor(n / 2) + 1L]
  mad <- sort(abs(x - med))[floor(n / 2) + 1L] * 1.4826
  if (!(mad > 1e-9)) {
    m <- mean(x); sdv <- stats::sd(x); if (!(sdv > 0)) sdv <- 1
    return((x - m) / sdv)
  }
  (x - med) / mad
}

# Profile informativeness = Var(item means) / mean(item variance) on P.
.profile_info <- function(P) {
  cm <- colMeans(P)
  mv <- mean(apply(P, 2, stats::var))
  vm <- mean((cm - mean(cm))^2)
  if (mv > 0) vm / mv else 0
}

# Partner reliability gate g in [0, 1] from profile informativeness.
.gate <- function(profile_info) {
  max(0, min(1, (profile_info - .RR$gate_offset) / .RR$gate_scale))
}
