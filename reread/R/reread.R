#' Detect careless responding (the reread procedure)
#'
#' The main entry point. Scores every respondent with the the reread procedure ensemble
#' (\emph{rc} + longstring + person-total, reliability-gated logistic combiner with
#' frozen weights), flags careless respondents with an automatic, label-free
#' threshold, and estimates the careless prevalence. One call from a matrix or
#' data frame of item responses.
#'
#' @inheritParams rc_index
#' @param score Which score the threshold is applied to. \code{"ensemble"}
#'   (default) cuts the reread ensemble log-odds. \code{"person_total"} cuts the
#'   oriented person-total correlation instead, using the identical calibration.
#'   On external data the two are close to interchangeable in both ranking and
#'   decision quality, so the alternative is offered rather than hidden; the
#'   ensemble is the default because it proved the more robust of the two --- it
#'   separates far better on collected data with clean ground truth, and it
#'   removes a steadier share of the sample --- and because the probability scale
#'   and the borderline band are defined only for it (both are \code{NA} under
#'   \code{"person_total"}).
#' @param sensitivity The shipped cut, \code{"standard"} (a conservative
#'   \eqn{2.5\sigma} tail on the attentive mode). It is the only preset and is the
#'   better choice up to a careless rate of about 35\%; for samples beyond that,
#'   pass a numeric \code{z} to \code{\link{auto_flag}} or supply
#'   \code{prevalence}.
#' @param prevalence Optional known/expected careless rate (a fraction in (0, 1)).
#'   If supplied, the automatic threshold is bypassed and the top
#'   \code{prevalence} share of respondents by score is flagged instead.
#' @param verbose If \code{TRUE}, print the diagnostic warnings.
#'
#' @return An object of class \code{"reread"}: a list with \code{flagged} (a
#'   logical vector), \code{score} (which score was cut),
#'   \code{scores} (a data frame of per-respondent \code{rc},
#'   \code{longstring}, \code{person_total}, the log-odds \code{eta} and its
#'   Monte Carlo standard error \code{eta_se}, the careless \code{prob} and its
#'   error \code{prob_se}, \code{flagged}, and \code{borderline} -- whether the
#'   score lies within about two standard errors of the flagging threshold), the
#'   calibrated \code{prevalence}, \code{n_flagged}, \code{n_borderline},
#'   \code{flagged_fraction}, \code{threshold}, \code{sensitivity}, the reliability
#'   \code{gate}, \code{profile_info}, \code{signal_strength},
#'   \code{two_component}, and any \code{warnings}. See \code{\link{print.reread}},
#'   \code{\link{summary.reread}}, and \code{\link{plot.reread}}.
#'
#' @examples
#' set.seed(1)
#' clean <- simulate_clean(n_factors = 12, items_per_factor = 6, n = 300)
#' dat   <- inject_careless(clean$data, prevalence = 0.2)
#' fit   <- reread(dat$data)
#' fit
#' head(fit$scores)
#' # the same label-free cut applied to person-total alone
#' fit_pt <- reread(dat$data, score = "person_total")
#' table(ensemble = fit$flagged, person_total = fit_pt$flagged)
#' @seealso \code{\link{rc_index}}, \code{\link{ensemble_score}}, \code{\link{auto_flag}}
#' @export
reread <- function(data, score = c("ensemble", "person_total"),
                     sensitivity = "standard",
                     prevalence = NULL, rescale = c("proportion", "none", "minmax", "zscore"),
                     iterations = 200L, cor_prop = 0.03, min_pairs = 15L,
                     seed = 1L, verbose = FALSE) {
  score <- match.arg(score)
  sensitivity <- .resolve_sensitivity(sensitivity)
  rescale <- match.arg(rescale)
  mat <- .as_matrix(data)
  n <- nrow(mat); J <- ncol(mat)

  f <- .compute_features(mat, cor_prop, iterations, min_pairs, seed, rescale,
                         .RR$weights)
  ## The threshold is a property of the calibration, not of the ensemble: it is
  ## equivariant to affine transformations of whatever score it is given, so the
  ## same rule can be applied to a single index. `score = "person_total"` cuts the
  ## oriented person-total correlation instead of the ensemble log-odds -- the two
  ## are close to interchangeable on external data (see the package vignette and
  ## the paper's head-to-head), and this makes that choice available rather than
  ## implicit. The ensemble remains the default: it was the more robust of the two
  ## on collected data with clean ground truth, and it removes a steadier share of
  ## the sample.
  eta <- if (score == "ensemble") f$eta else f$o_person_total

  if (!is.null(prevalence)) {
    if (!is.numeric(prevalence) || prevalence <= 0 || prevalence >= 1)
      stop("`prevalence` must be a fraction in (0, 1).")
    k <- max(0L, min(n, round(prevalence * n)))
    thr <- if (k == 0L) Inf else sort(eta, decreasing = TRUE)[k]
    flagged <- eta >= thr & k > 0L
    fl <- list(flagged = flagged, n_flagged = sum(flagged),
               flagged_fraction = mean(flagged), prevalence = prevalence,
               threshold = if (k == 0L) NA_real_ else thr, two_component = NA,
               separation = NA_real_)
  } else {
    fl <- auto_flag(eta, sensitivity)
  }

  warns <- character(0)
  if (score == "ensemble") {
    if (f$signal_strength < 1.3)
      warns <- c(warns, "No usable multi-construct structure (signal strength < 1.3): rc is not suitable for this dataset; rely on longstring / IRV / response times.")
    else if (f$signal_strength < 1.8)
      warns <- c(warns, "Weak multi-construct structure: treat scores as a ranking, not a hard flag, and corroborate with other indices.")
    if (J < 60L)
      warns <- c(warns, sprintf("Short questionnaire (%d items): below the method's recommended envelope (~60+ items); rc is weaker here and the partners carry more of the signal.", J))
    if (f$gate <= 0)
      warns <- c(warns, "Flat item-mean profile: the partner detectors were gated off and the score fell back to rc alone.")
  } else {
    warns <- c(warns, "Scoring on person-total alone: the ensemble is the default because it was the more robust of the two, and it removes a steadier share of the sample. Person-total is offered because on external data the two are close to interchangeable; prefer it only with a reason.")
    if (f$profile_info < 0.02)
      warns <- c(warns, "Flat item-mean profile: person-total is correlated against a nearly flat reference here, so it carries little information on this battery. The ensemble would gate it off.")
  }
  if (is.null(prevalence) && fl$prevalence > 0.35)
    warns <- c(warns, "Estimated prevalence above ~35%: beyond the range where the shipped cut is the best choice; it is deliberately conservative and may under-flag here. Consider auto_flag(eta, z = 1.5) or supplying `prevalence`.")
  if (is.null(prevalence) && identical(fl$status, "none-past-cut"))
    warns <- c(warns, "No respondent reaches the cut, so nothing was flagged. That is the expected result on clean data, but carelessness that shifts a profile slightly without producing extreme scores looks the same: use the continuous ranking, or supply `prevalence`, before concluding the data are clean.")

  ## prob and the borderline band are properties of the fitted ensemble: the
  ## logistic scale is calibrated only for eta, and the only stochastic component
  ## is rc's permutation baseline. Scoring on person-total alone, both are
  ## undefined and are reported as NA rather than as a number that looks
  ## comparable but is not.
  is_ens <- score == "ensemble"
  prob <- if (is_ens) stats::plogis(eta) else rep(NA_real_, n)
  eta_se <- if (is_ens) f$eta_se else rep(NA_real_, n)
  thr <- fl$threshold
  borderline <- if (is_ens && is.finite(thr)) abs(eta - thr) < 2 * f$eta_se
                else rep(FALSE, n)
  scores <- data.frame(
    rc = f$rc, longstring = f$longstring, person_total = f$person_total,
    eta = eta, eta_se = eta_se, prob = prob,
    prob_se = if (is_ens) f$eta_se * prob * (1 - prob) else rep(NA_real_, n),
    flagged = fl$flagged, borderline = borderline)

  out <- list(
    flagged = fl$flagged, score = score, scores = scores, prevalence = fl$prevalence,
    n_flagged = fl$n_flagged, n_borderline = sum(borderline),
    flagged_fraction = fl$flagged_fraction,
    threshold = fl$threshold, two_component = fl$two_component,
    sensitivity = if (is.null(prevalence)) sensitivity else "manual",
    gate = f$gate, profile_info = f$profile_info,
    signal_strength = f$signal_strength, separation = fl$separation,
    n = n, J = J, warnings = warns, call = match.call())
  class(out) <- "reread"
  if (verbose && length(warns)) for (w in warns) warning(w, call. = FALSE)
  out
}
