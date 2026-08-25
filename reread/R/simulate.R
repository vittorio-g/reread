#' Simulate clean questionnaire data
#'
#' Generates attentive (clean) multi-factor Likert responses from a common-factor
#' model: random loadings (about 40\% reverse-keyed), a weakly correlated factor
#' structure, per-item intercepts so the item-mean profile resembles a real
#' questionnaire, and equiprobable Likert discretization. Useful for examples,
#' teaching, and testing (combine with \code{\link{inject_careless}}).
#'
#' @param n_factors Number of latent factors.
#' @param items_per_factor Items per factor (total items = \code{n_factors *
#'   items_per_factor}).
#' @param n Number of respondents.
#' @param rho Common-factor correlation (default 0.3).
#' @param n_cat Number of Likert categories (default 5).
#' @param seed Optional integer seed.
#' @return A list with \code{data} (an \code{n} by \code{J} integer matrix of
#'   responses), the item \code{loadings}, the \code{factor} index of each item,
#'   and a logical \code{reverse} (reverse-keyed items).
#' @examples
#' d <- simulate_clean(n_factors = 10, items_per_factor = 6, n = 200, seed = 1)
#' dim(d$data)
#' @export
simulate_clean <- function(n_factors, items_per_factor, n, rho = 0.3,
                           n_cat = 5L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  F <- n_factors; ipf <- items_per_factor; J <- F * ipf; K <- n_cat
  thr <- stats::qnorm(seq_len(K - 1L) / K)
  load <- 0.4 + 0.4 * stats::runif(J)
  flip <- stats::runif(J) < 0.4; load[flip] <- -load[flip]
  fac <- rep(seq_len(F), each = ipf)
  tau <- 0.45 * stats::rnorm(J)
  sq <- sqrt(rho); sq1 <- sqrt(1 - rho)
  data <- matrix(0L, n, J)
  for (i in seq_len(n)) {
    common <- stats::rnorm(1); f <- sq * common + sq1 * stats::rnorm(F)
    y <- tau + load * f[fac] + sqrt(1 - load^2) * stats::rnorm(J)
    data[i, ] <- 1L + rowSums(outer(y, thr, ">"))
  }
  colnames(data) <- sprintf("i%02d", seq_len(J))
  list(data = data, loadings = load, factor = fac, reverse = flip)
}

# make blocks of identical values over the k items in `idx` (longstring/fatigue).
.inject_chunks <- function(row, idx, levels, sequential) {
  k <- length(idx)
  nc <- if (k >= 4L) sample(2:min(5L, k), 1L) else if (k >= 2L) 2L else 1L
  asg <- (seq_len(k) - 1L) %% nc
  if (!sequential) asg <- sample(asg)
  for (ch in seq_len(nc) - 1L) {
    v <- sample(levels, 1L); row[idx[asg == ch]] <- v
  }
  row
}

#' Inject careless responding into clean data
#'
#' Corrupts a share of respondents with one of six prototypical careless patterns
#' -- \code{"random"}, \code{"longstring"}, \code{"pure_straight"},
#' \code{"acquiescent"}, \code{"mixed"}, and \code{"fatigue"} -- each at a random
#' corruption fraction. Two of the six (\code{pure_straight}, \code{acquiescent})
#' are \emph{consistent} patterns that rc is deliberately blind to, so the set is
#' not stacked in the ensemble's favour.
#'
#' @param data A clean numeric matrix or data frame (e.g. from
#'   \code{\link{simulate_clean}}).
#' @param prevalence Fraction of respondents to make careless.
#' @param patterns Character vector of patterns to draw from (default: all six).
#' @param corruption Numeric vector of corruption fractions to draw from (default
#'   0.5 to 1.0).
#' @param seed Optional integer seed.
#' @return A list with the corrupted \code{data}, a 0/1 \code{labels} vector, a
#'   \code{severity} vector (corruption fraction, 0 for clean), and the
#'   \code{pattern} applied to each careless respondent (\code{NA} for clean).
#' @examples
#' clean <- simulate_clean(10, 6, 200, seed = 1)
#' inj <- inject_careless(clean$data, prevalence = 0.2, seed = 2)
#' table(inj$labels)
#' @export
inject_careless <- function(data, prevalence = 0.2,
                            patterns = c("random", "longstring", "pure_straight",
                                         "acquiescent", "mixed", "fatigue"),
                            corruption = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  patterns <- match.arg(patterns, several.ok = TRUE)
  mat <- .as_matrix(data); N <- nrow(mat); J <- ncol(mat)
  lo <- min(mat, na.rm = TRUE); hi <- max(mat, na.rm = TRUE)
  levels <- lo:hi
  high_levels <- levels[levels >= lo + 0.65 * (hi - lo)]
  clamp <- function(x) pmax(lo, pmin(hi, x))
  all_idx <- seq_len(J)

  M <- round(N * prevalence)
  ids <- sample.int(N, M)
  labels <- integer(N); sev <- numeric(N); pat_out <- rep(NA_character_, N)
  out <- mat
  for (m in seq_len(M)) {
    i <- ids[m]; pat <- patterns[((m - 1L) %% length(patterns)) + 1L]
    lvl <- sample(corruption, 1L); k <- min(J, max(1L, round(J * lvl)))
    row <- out[i, ]
    if (pat == "random") {
      idx <- sample(all_idx, k); row[idx] <- sample(levels, k, replace = TRUE)
    } else if (pat == "longstring") {
      row <- .inject_chunks(row, sample(all_idx, k), levels, FALSE)
    } else if (pat == "pure_straight") {
      v <- sample(levels, 1L); row[sample(all_idx, k)] <- v
    } else if (pat == "acquiescent") {
      b <- sample(high_levels, 1L)
      idx <- sample(all_idx, k)
      row[idx] <- clamp(b + sample(c(-1L, 0L, 0L, 0L, 1L), k, replace = TRUE))
    } else if (pat == "mixed") {
      kl <- max(1L, round(k / 2)); idx <- sample(all_idx, k)
      row <- .inject_chunks(row, idx[seq_len(kl)], levels, FALSE)
      rest <- idx[(kl + 1L):k]; row[rest] <- sample(levels, length(rest), replace = TRUE)
    } else if (pat == "fatigue") {
      idx <- (J - k + 1L):J; row <- .inject_chunks(row, idx, levels, TRUE)
    }
    out[i, ] <- row; labels[i] <- 1L; sev[i] <- lvl; pat_out[i] <- pat
  }
  list(data = out, labels = labels, severity = sev, pattern = pat_out)
}
