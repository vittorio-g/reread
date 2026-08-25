# Auxiliary careless-detection indices, computed on the proportion-rescaled
# matrix so they match the ensemble's internal features and the browser tool.

#' Longstring index
#'
#' The longest run of identical consecutive responses for each respondent, in the
#' column order of \code{data}. High values indicate straight-lining.
#'
#' @param data A matrix or data frame of numeric item responses.
#' @return A numeric vector, one value per respondent (high = careless).
#' @examples
#' d <- simulate_clean(6, 6, 100)
#' longstring(inject_careless(d$data, prevalence = 0.2)$data)
#' @export
longstring <- function(data) {
  mat <- .as_matrix(data)
  vapply(seq_len(nrow(mat)), function(i) {
    v <- mat[i, ]; v <- v[!is.na(v)]
    if (length(v) < 2L) return(NA_real_)
    max(rle(as.numeric(v))$lengths)
  }, numeric(1))
}

#' Person-total correlation
#'
#' The correlation between each respondent's profile and the sample's mean item
#' profile. \strong{Low} values indicate answering against the grain (careless).
#'
#' @param data A matrix or data frame of numeric item responses.
#' @return A numeric vector, one value per respondent (low = careless).
#' @export
person_total <- function(data) {
  P <- .prepare(.as_matrix(data))$P
  cm <- colMeans(P)
  apply(P, 1, function(r) suppressWarnings(stats::cor(r, cm)))
}

#' Intra-individual response variability (IRV)
#'
#' The within-person standard deviation of the (rescaled) responses. Careless
#' respondents are at the extremes: very low for straight-lining, high for
#' erratic responding.
#'
#' @param data A matrix or data frame of numeric item responses.
#' @return A numeric vector, one value per respondent.
#' @export
irv <- function(data) {
  P <- .prepare(.as_matrix(data))$P
  apply(P, 1, stats::sd)
}

#' Mahalanobis squared distance
#'
#' The Mahalanobis \eqn{D^2} of each respondent from the sample centroid, using a
#' Moore-Penrose pseudo-inverse of the covariance (so it is defined when the
#' number of items approaches or exceeds the number of respondents). High values
#' indicate multivariate outliers. Requires the \pkg{MASS} package.
#'
#' @param data A matrix or data frame of numeric item responses.
#' @return A numeric vector, one value per respondent (high = careless).
#' @export
mahalanobis_d2 <- function(data) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("mahalanobis_d2() requires the 'MASS' package.")
  X <- .prepare(.as_matrix(data))$P
  X <- X[, apply(X, 2, stats::sd) > 0, drop = FALSE]
  cm <- colMeans(X); Si <- MASS::ginv(stats::cov(X))
  apply(X, 1, function(r) { d <- r - cm; as.numeric(t(d) %*% Si %*% d) })
}
