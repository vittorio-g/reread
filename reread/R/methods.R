#' Methods for the reread procedure results
#'
#' \code{print} gives a one-screen summary; \code{summary} adds the component-score
#' quantiles, the flagged rows, and any diagnostic notes; \code{as.data.frame}
#' returns the per-respondent score table.
#'
#' @param x,object A \code{"reread"} object (from \code{\link{reread}}).
#' @param ... Ignored.
#' @return \code{print} and \code{summary} return their argument invisibly;
#'   \code{as.data.frame} returns the per-respondent score data frame.
#' @name reread-methods
#' @export
print.reread <- function(x, ...) {
  cat("the reread procedure careless-response detection\n")
  cat(sprintf("  Respondents: %d   Items: %d\n", x$n, x$J))
  if (!is.null(x$score) && !identical(x$score, "ensemble"))
    cat(sprintf("  Score cut: %s (not the ensemble)\n", x$score))
  sens <- if (is.numeric(x$sensitivity)) sprintf("z = %.2f", x$sensitivity) else x$sensitivity
  cat(sprintf("  Flagged: %d (%.1f%%)   Sensitivity: %s\n",
              x$n_flagged, 100 * x$flagged_fraction, sens))
  cat(sprintf("  Estimated careless prevalence: %.1f%%%s\n", 100 * x$prevalence,
              if (isTRUE(x$prevalence > 0.25)) " (lower bound)" else ""))
  cat(sprintf("  Reliability gate: %.2f   Signal strength: %.2f\n",
              x$gate, x$signal_strength))
  if (!is.null(x$n_borderline) && x$n_borderline > 0)
    cat(sprintf("  Borderline flags (within ~2 SE of the threshold): %d\n", x$n_borderline))
  if (length(x$warnings))
    cat(sprintf("  %d note%s (see summary())\n", length(x$warnings),
                if (length(x$warnings) > 1) "s" else ""))
  invisible(x)
}

#' @rdname reread-methods
#' @export
summary.reread <- function(object, ...) {
  print(object)
  cat("\nComponent scores (higher eta = more careless):\n")
  s <- object$scores
  tab <- rbind(
    rc           = stats::quantile(s$rc, c(0, .5, 1), na.rm = TRUE),
    longstring   = stats::quantile(s$longstring, c(0, .5, 1), na.rm = TRUE),
    person_total = stats::quantile(s$person_total, c(0, .5, 1), na.rm = TRUE),
    eta          = stats::quantile(s$eta, c(0, .5, 1), na.rm = TRUE))
  colnames(tab) <- c("min", "median", "max")
  print(round(tab, 3))
  if (object$n_flagged > 0) {
    cat(sprintf("\nFlagged respondents (%d): rows %s\n", object$n_flagged,
        paste(utils::head(which(object$flagged), 10), collapse = ", ")))
    if (object$n_flagged > 10) cat("  ... (first 10 shown)\n")
  }
  if (length(object$warnings)) {
    cat("\nNotes:\n")
    for (w in object$warnings) cat(strwrap(paste0("- ", w), exdent = 2), sep = "\n")
    cat("\n")
  }
  invisible(object)
}

#' @rdname reread-methods
#' @export
as.data.frame.reread <- function(x, ...) x$scores

#' Plot a the reread procedure result
#'
#' Draws the distribution of the ensemble log-odds with the flagging threshold
#' and the flagged region highlighted. Requires the \pkg{ggplot2} package.
#'
#' @param x A \code{"reread"} object.
#' @param bins Number of histogram bins (default 40).
#' @param ... Ignored.
#' @return A \pkg{ggplot2} object (invisibly printed).
#' @examples
#' \donttest{
#' d <- inject_careless(simulate_clean(12, 6, 300)$data, prevalence = 0.2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(reread(d$data))
#' }
#' @export
plot.reread <- function(x, bins = 40, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot.reread() requires the 'ggplot2' package.")
  df <- data.frame(eta = x$scores$eta, flagged = x$scores$flagged)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, fill = .data$flagged)) +
    ggplot2::geom_histogram(bins = bins, colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#8a8a8a", `TRUE` = "#8f3535"),
                               labels = c("attentive", "flagged"), name = NULL) +
    ggplot2::labs(x = "Ensemble score (log-odds, higher = more careless)",
                  y = "Respondents",
                  title = "the reread procedure: score distribution and flags") +
    ggplot2::theme_minimal(base_size = 12)
  if (is.finite(x$threshold))
    p <- p + ggplot2::geom_vline(xintercept = x$threshold, linetype = "dashed",
                                 colour = "#33475f")
  print(p)
  invisible(p)
}
