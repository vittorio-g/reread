#' Benchmark the reread procedure against other careless indices
#'
#' Computes a set of careless-detection indices on the same data, orients them so
#' that higher means more careless, and compares them by the area under the ROC
#' curve against a supplied ground-truth label, with DeLong confidence intervals
#' and paired DeLong tests against a reference method. Reproduces the style of the
#' head-to-head comparison in the paper. Requires the \pkg{pROC} package.
#'
#' @param data A matrix or data frame of item responses.
#' @param truth A 0/1 (or logical) vector of true careless labels, one per
#'   respondent.
#' @param methods Indices to include; any of \code{"reread"}, \code{"rc"},
#'   \code{"longstring"}, \code{"person_total"}, \code{"irv"},
#'   \code{"mahalanobis_d2"}. Defaults to all.
#' @param reference Method to test the others against with the DeLong test
#'   (default \code{"reread"}).
#' @param seed Seed passed to the rc engine (default 1).
#' @return A data frame with, per method, the \code{auc} and its 95\% DeLong
#'   confidence interval (\code{lo}, \code{hi}), and the DeLong \code{p} value
#'   versus the reference method, sorted by AUC.
#' @examples
#' \donttest{
#' clean <- simulate_clean(12, 6, 300, seed = 1)
#' inj <- inject_careless(clean$data, prevalence = 0.25, seed = 2)
#' if (requireNamespace("pROC", quietly = TRUE))
#'   benchmark_indices(inj$data, inj$labels)
#' }
#' @export
benchmark_indices <- function(data, truth,
                              methods = c("reread", "rc", "longstring",
                                          "person_total", "irv", "mahalanobis_d2"),
                              reference = "reread", seed = 1L) {
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("benchmark_indices() requires the 'pROC' package.")
  methods <- match.arg(methods, several.ok = TRUE)
  mat <- .as_matrix(data)
  truth <- as.integer(truth)
  if (length(truth) != nrow(mat)) stop("`truth` must have one label per respondent.")

  # oriented scores (higher = more careless)
  scores <- list()
  if (any(c("reread", "rc") %in% methods)) {
    f <- .compute_features(mat, seed = seed)
    if ("reread" %in% methods) scores[["reread"]] <- f$eta
    if ("rc" %in% methods)       scores[["rc"]] <- -f$rc
  }
  if ("longstring" %in% methods)     scores[["longstring"]]   <- longstring(mat)
  if ("person_total" %in% methods)   scores[["person_total"]] <- -person_total(mat)
  if ("irv" %in% methods)            scores[["irv"]]          <- -irv(mat)
  if ("mahalanobis_d2" %in% methods) scores[["mahalanobis_d2"]] <- mahalanobis_d2(mat)

  roc_of <- function(s) {
    ok <- is.finite(s)
    if (length(unique(truth[ok])) < 2) return(NULL)
    pROC::roc(truth[ok], s[ok], direction = "<", quiet = TRUE)
  }
  rocs <- lapply(scores, roc_of)
  ref_roc <- if (reference %in% names(rocs)) rocs[[reference]] else NULL

  rows <- lapply(names(rocs), function(m) {
    r <- rocs[[m]]
    if (is.null(r)) return(data.frame(method = m, auc = NA, lo = NA, hi = NA, p = NA))
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    p <- NA_real_
    if (!is.null(ref_roc) && m != reference) {
      ok <- is.finite(scores[[m]]) & is.finite(scores[[reference]])
      p <- tryCatch(pROC::roc.test(
        pROC::roc(truth[ok], scores[[reference]][ok], direction = "<", quiet = TRUE),
        pROC::roc(truth[ok], scores[[m]][ok], direction = "<", quiet = TRUE),
        method = "delong", paired = TRUE)$p.value, error = function(e) NA_real_)
    }
    data.frame(method = m, auc = round(as.numeric(pROC::auc(r)), 3),
               lo = round(ci[1], 3), hi = round(ci[3], 3), p = signif(p, 3))
  })
  res <- do.call(rbind, rows)
  res <- res[order(-res$auc), ]
  rownames(res) <- NULL
  attr(res, "reference") <- reference
  res
}
