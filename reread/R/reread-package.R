#' reread: Careless Respondent Recognition via Resampled Ensemble
#'
#' Detects careless or insufficient-effort responding in questionnaire data. The
#' method is built on \emph{rc} (resampled reliability), a permutation-based index
#' of individual response coherence that needs no declared factor structure,
#' combined with two complementary detectors (longstring and person-total
#' correlation) in a compact, reliability-gated logistic ensemble. A BIC-gated
#' two-component mixture sets the flagging threshold without labels and returns a
#' calibrated prevalence estimate.
#'
#' The main entry point is \code{\link{reread}}. The components are exposed as
#' \code{\link{rc_index}}, \code{\link{ensemble_score}}, and \code{\link{auto_flag}};
#' \code{\link{simulate_clean}} / \code{\link{inject_careless}} generate test data;
#' and \code{\link{benchmark_indices}} compares against other indices.
#'
#' @keywords internal
"_PACKAGE"

# `.data` is the ggplot2/rlang pronoun used in plot.reread(); declaring it here
# silences the "undefined global" NOTE without importing rlang.
utils::globalVariables(".data")
