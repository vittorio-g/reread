# The rc engine: prepare -> pair correlations -> coupled coherence ->
# permutation baseline -> per-respondent z. Returns everything the ensemble
# needs so nothing is recomputed.
.rc_engine <- function(mat, cor_prop = .RR$cor_prop, iterations = .RR$iterations,
                       min_pairs = .RR$min_pairs, seed = 1L,
                       rescale = "proportion") {
  n <- nrow(mat); J <- ncol(mat)
  if (J < 6L)  stop("Need at least 6 item columns (found ", J, ").")
  if (n < 10L) stop("Need at least 10 respondents (found ", n, ").")

  prep <- .prepare(mat, rescale); P <- prep$P; minP <- prep$minP; maxP <- prep$maxP
  cp <- .correlation_pairs(P); pa <- cp$pa; pb <- cp$pb; pr <- cp$pr; nP <- cp$n_pairs

  k <- min(max(min_pairs, round(cor_prop * nP)), nP)
  ord <- order(abs(pr), decreasing = TRUE)
  top <- ord[seq_len(k)]
  coupled <- .ind_cors(P, minP, maxP, pa, pb, pr, top)

  # structure diagnostic (noise-calibrated signal strength)
  top_mean <- mean(abs(pr[top])); noise_top <- 2.4 / sqrt(max(n, 2))
  strength <- top_mean / noise_top

  # permutation baseline
  if (!is.null(seed)) set.seed(seed)
  s <- numeric(n); ss <- numeric(n)
  for (b in seq_len(iterations)) {
    rc <- .ind_cors(P, minP, maxP, pa, pb, pr, sample.int(nP, k))
    s <- s + rc; ss <- ss + rc * rc
  }
  mu <- s / iterations
  va <- pmax(ss / iterations - mu * mu, 0)
  sdv <- sqrt(va)
  z <- ifelse(sdv > 1e-12, (coupled - mu) / sdv, NA_real_)
  # degenerate respondents (no coherence baseline) -> most-careless observed z
  z[!is.finite(z)] <- if (any(is.finite(z))) min(z[is.finite(z)]) else 0
  # per-respondent Monte Carlo standard error of the permutation z-score
  z_se <- sqrt((1 + z^2 / 2) / iterations)

  list(z = z, z_se = z_se, coupled = coupled, k = k, n_pairs = nP, P = P, minP = minP,
       signal_strength = strength, top_mean_r = top_mean)
}

#' The rc index of individual response coherence
#'
#' Computes \emph{rc} (resampled reliability), a permutation-based per-respondent
#' z-score of how coherently each respondent answered the sample's most strongly
#' correlated item pairs, relative to a within-person null built from random item
#' pairs. \strong{Low values indicate potential carelessness}; a high value means
#' the respondent is more coherent than their own chance baseline.
#'
#' The index needs no declared factor structure: the coupled pairs are read off
#' the data's own correlation matrix. Reverse-keyed items are handled
#' automatically (a pair with negative sample correlation has its partner
#' reverse-coded), and mixed response scales are handled by the \code{rescale}
#' step. Respondents with no response variance over the coupled items
#' (straight-liners) receive the most-careless score.
#'
#' @param data A matrix or data frame of numeric item responses (respondents in
#'   rows, items in columns; at least 6 items and 10 respondents). Missing values
#'   are median-imputed per item.
#' @param cor_prop Proportion of top-\eqn{|r|} item pairs retained as coupled
#'   pairs (default 0.03).
#' @param iterations Number of permutation rounds for the null (default 200).
#' @param min_pairs Floor on the number of coupled pairs (default 15).
#' @param rescale Per-item rescaling before correlations: \code{"proportion"}
#'   (default, divide by the item maximum; leaves homogeneous scales unchanged),
#'   \code{"none"}, \code{"minmax"}, or \code{"zscore"}.
#' @param seed Integer seed for the permutation baseline (default 1); set to
#'   \code{NULL} to leave the RNG untouched.
#'
#' @return A numeric vector of \emph{rc} z-scores, one per respondent (low =
#'   careless), with attributes \code{n_pairs}, \code{k} (coupled pairs used),
#'   and \code{signal_strength} (a noise-calibrated measure of usable structure;
#'   values below about 1.3 mean the questionnaire has no multi-construct
#'   structure and rc is not suitable).
#'
#' @examples
#' set.seed(1)
#' d <- simulate_clean(n_factors = 8, items_per_factor = 6, n = 200)
#' rc <- rc_index(d$data)
#' summary(rc)
#' @seealso \code{\link{reread}} for the full ensemble.
#' @export
rc_index <- function(data, cor_prop = 0.03, iterations = 400L, min_pairs = 15L,
                     rescale = c("proportion", "none", "minmax", "zscore"),
                     seed = 1L) {
  rescale <- match.arg(rescale)
  mat <- .as_matrix(data)
  eng <- .rc_engine(mat, cor_prop, iterations, min_pairs, seed, rescale)
  z <- eng$z
  attr(z, "n_pairs") <- eng$n_pairs
  attr(z, "k") <- eng$k
  attr(z, "signal_strength") <- eng$signal_strength
  attr(z, "se") <- eng$z_se
  z
}
