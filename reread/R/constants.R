# Frozen constants of the shipped the reread procedure model (n = 157, 109 items).
# These are the exact values used by the browser tool at ca.re-re.re and reported
# in the paper; the package reproduces its scores. Do not edit without refitting
# and re-validating (see `reread_weights`).

.RR <- list(
  # Logistic weights on the oriented robust-z features (higher = more careless).
  # irv and d2 carry zero weight (see the package documentation).
  weights = c(b0 = -0.1057, rc = 2.7645, longstring = 1.2689,
              person_total = 1.4452, irv = 0, d2 = 0),
  # Reliability gate: g = clamp((profile_info - offset) / scale, 0, 1).
  gate_offset = 0.01,
  gate_scale  = 0.05,
  # The single shipped flagging cut: attentive-tail multiplier z. Anchoring the
  # cut on the one-sided attentive fit extends its useful range to a ~35% careless
  # rate, so the second ("high", 1.5) preset of the earlier two-setting scheme was
  # removed; `z` in auto_flag() remains for the out-of-envelope case.
  sens_z = c(standard = 2.5),
  # rc engine defaults. The operational permutation count is 400 (linear cost);
  # the per-respondent Monte Carlo standard error of the rc z-score is about
  # sqrt((1 + z^2 / 2) / iterations).
  cor_prop = 0.03, min_pairs = 15L, iterations = 400L
)

#' Frozen model constants
#'
#' The logistic weights, reliability-gate parameters, and sensitivity presets of
#' the shipped the reread procedure model, exposed for inspection and citation. These values
#' are fixed (fitted and leave-one-out validated once on the validation study of
#' Guerrieri, Gallucci, and Passarelli, 2026) and are applied unchanged to every
#' questionnaire; the package does not refit them.
#'
#' @format A named list with elements \code{weights} (the logistic coefficients
#'   on the oriented robust-z features), \code{gate_offset} and \code{gate_scale}
#'   (the reliability gate), \code{sens_z} (the shipped attentive-tail
#'   multiplier), and the \code{rc} engine defaults.
#' @examples
#' reread_weights$weights
#' @export
reread_weights <- .RR
