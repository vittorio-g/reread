# Column-order longstring on a (rescaled) matrix: longest run of equal values.
.longstring_mat <- function(P) {
  vapply(seq_len(nrow(P)), function(i) {
    v <- P[i, ]; v <- v[!is.na(v)]
    if (length(v) < 2L) return(NA_real_)
    max(rle(v)$lengths)
  }, numeric(1))
}

# Compute the ensemble features, orient them (higher = more careless), apply the
# reliability gate, and return the log-odds. The single place features are built.
.compute_features <- function(mat, cor_prop = .RR$cor_prop, iterations = .RR$iterations,
                              min_pairs = .RR$min_pairs, seed = 1L,
                              rescale = "proportion", weights = .RR$weights) {
  eng <- .rc_engine(mat, cor_prop, iterations, min_pairs, seed, rescale)
  P <- eng$P
  ls_raw <- .longstring_mat(P)
  cm <- colMeans(P)
  pt_raw <- apply(P, 1, function(r) suppressWarnings(stats::cor(r, cm)))
  irv_raw <- apply(P, 1, stats::sd)
  # A zero-variance (straight-lining) respondent has an undefined person-total
  # correlation; person-total is oriented low = careless, so impute the most
  # careless finite value. Guard the others for finiteness too.
  .impute_low <- function(v) { fin <- is.finite(v)
    if (!all(fin)) v[!fin] <- if (any(fin)) min(v[fin]) else 0; v }
  .impute_high <- function(v) { fin <- is.finite(v)
    if (!all(fin)) v[!fin] <- if (any(fin)) max(v[fin]) else 0; v }
  pt_raw  <- .impute_low(pt_raw)       # low  person-total  = careless
  ls_raw  <- .impute_high(ls_raw)      # high longstring     = careless
  irv_raw <- .impute_low(irv_raw)      # low  variability    = straight-lining

  o_rc <- -.robust_z(eng$z)             # low rc  = careless -> high feature
  o_ls <- .robust_z(ls_raw)             # high    = careless
  o_pt <- -.robust_z(pt_raw)            # low     = careless -> high feature

  pinfo <- .profile_info(P)
  g <- .gate(pinfo)
  W <- weights
  eta <- as.numeric(W[["b0"]] + W[["rc"]] * o_rc +
         g * (W[["longstring"]] * o_ls + W[["person_total"]] * o_pt))
  # Guarantee a finite score for every respondent: an unscoreable one (all
  # features degenerate) is treated as most careless.
  if (any(!is.finite(eta))) {
    fin <- is.finite(eta)
    eta[!fin] <- if (any(fin)) max(eta[fin]) else 0
  }
  # Propagate the rc Monte Carlo error to the log-odds. rc is the only stochastic
  # component; o_rc = -(z - med)/mad, so SE(eta) ~= |w_rc| * SE(z) / mad.
  nn <- length(eng$z)
  medz <- sort(eng$z)[floor(nn / 2) + 1L]
  madz <- sort(abs(eng$z - medz))[floor(nn / 2) + 1L] * 1.4826
  if (!(madz > 1e-9)) { madz <- stats::sd(eng$z); if (!(madz > 0)) madz <- 1 }
  eta_se <- abs(W[["rc"]]) * eng$z_se / madz

  list(eta = eta, eta_se = eta_se, rc = eng$z, rc_se = eng$z_se, longstring = ls_raw,
       person_total = pt_raw, irv = irv_raw,
       o_rc = o_rc, o_longstring = o_ls, o_person_total = o_pt,
       gate = g, profile_info = pinfo, signal_strength = eng$signal_strength,
       k = eng$k, n_pairs = eng$n_pairs, n = nrow(mat), J = ncol(mat))
}

#' Ensemble carelessness score (log-odds)
#'
#' Computes the the reread procedure ensemble score for each respondent: the log-odds
#' \eqn{\eta} of a reliability-gated logistic combination of \emph{rc}, longstring,
#' and person-total, using the frozen shipped weights (\code{\link{reread_weights}}).
#' Higher \eqn{\eta} means more careless. When the item-mean profile is flat
#' (uninformative partners), the reliability gate shrinks the two auxiliaries
#' toward zero and the score falls back to \emph{rc} alone.
#'
#' This function returns the continuous score only; use \code{\link{reread}}
#' for automatic flagging and a prevalence estimate.
#'
#' @inheritParams rc_index
#' @param weights Optional named vector of logistic weights (\code{b0}, \code{rc},
#'   \code{longstring}, \code{person_total}); defaults to the frozen shipped
#'   weights.
#' @return A list with \code{eta} (the per-respondent log-odds), a data frame
#'   \code{features} of the raw and oriented component scores, the reliability
#'   \code{gate} \eqn{g}, \code{profile_info}, and \code{signal_strength}.
#' @examples
#' d <- simulate_clean(10, 6, 200)
#' es <- ensemble_score(inject_careless(d$data, prevalence = 0.2)$data)
#' str(es$eta)
#' @seealso \code{\link{reread}}, \code{\link{auto_flag}}
#' @export
ensemble_score <- function(data, cor_prop = 0.03, iterations = 200L, min_pairs = 15L,
                           rescale = c("proportion", "none", "minmax", "zscore"),
                           seed = 1L, weights = reread_weights$weights) {
  rescale <- match.arg(rescale)
  f <- .compute_features(.as_matrix(data), cor_prop, iterations, min_pairs, seed,
                         rescale, weights)
  features <- data.frame(
    rc = f$rc, longstring = f$longstring, person_total = f$person_total, irv = f$irv,
    o_rc = f$o_rc, o_longstring = f$o_longstring, o_person_total = f$o_person_total)
  list(eta = f$eta, eta_se = f$eta_se, features = features, gate = f$gate,
       profile_info = f$profile_info, signal_strength = f$signal_strength)
}
