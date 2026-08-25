#' Two-component Gaussian mixture (EM)
#'
#' Fits a one-dimensional two-component Gaussian mixture to a numeric vector by
#' the EM algorithm and returns the fitted parameters together with the Bayesian
#' information criterion of the one- and two-component fits. This is a
#' \emph{diagnostic}, not a decision: \code{\link{auto_flag}} reports how far apart
#' a two-component solution puts the groups, but neither the flagging cut nor the
#' prevalence estimate comes from it (see \code{\link{left_fit}}). An earlier
#' version used this test as a gate that suppressed all flags when it found no
#' distinct careless component; that gate was removed after it was found to leave
#' the worst case on clean data unchanged while withholding detection at low
#' prevalence.
#'
#' @param x A numeric vector (typically the ensemble log-odds).
#' @return A list with the mixing proportion \code{pi} of the high (careless)
#'   component, component means \code{m1} < \code{m2} and variances \code{v1},
#'   \code{v2}, the per-element posterior \code{post} of the careless component,
#'   and \code{bic1}, \code{bic2} (one- vs two-component BIC; a two-component fit
#'   is preferred when \code{bic1 - bic2 > 2}).
#' @examples
#' x <- c(rnorm(180, 0), rnorm(20, 4))
#' gauss_mixture(x)$pi
#' @export
gauss_mixture <- function(x) {
  x <- as.numeric(x); n <- length(x)
  mu <- mean(x); varr <- mean((x - mu)^2)
  sd <- sqrt(varr); if (!(sd > 0)) sd <- 1
  vfloor <- max(1e-4, varr * 0.02)
  m1 <- mu - 0.7 * sd; m2 <- mu + 0.7 * sd; v1 <- varr; v2 <- varr; pihat <- 0.3
  dn <- function(v, m, s2) exp(-((v - m)^2) / (2 * s2)) / sqrt(2 * pi * s2)
  post <- numeric(n)
  for (it in seq_len(300)) {
    a <- pihat * dn(x, m2, v2); b <- (1 - pihat) * dn(x, m1, v1)
    post <- a / (a + b + 1e-300); sp <- sum(post)
    pihat <- min(0.999, max(1e-3, sp / n))
    m2 <- sum(post * x) / (sp + 1e-9); m1 <- sum((1 - post) * x) / (n - sp + 1e-9)
    v2 <- max(vfloor, sum(post * (x - m2)^2) / (sp + 1e-9))
    v1 <- max(vfloor, sum((1 - post) * (x - m1)^2) / (n - sp + 1e-9))
  }
  if (m2 < m1) { t <- m1; m1 <- m2; m2 <- t; t <- v1; v1 <- v2; v2 <- t
                 pihat <- 1 - pihat; post <- 1 - post }
  ll2 <- sum(log(pihat * dn(x, m2, v2) + (1 - pihat) * dn(x, m1, v1) + 1e-300))
  ll1 <- sum(log(dn(x, mu, varr) + 1e-300))
  list(pi = pihat, post = post, m1 = m1, m2 = m2, v1 = v1, v2 = v2,
       bic1 = -2 * ll1 + 2 * log(n), bic2 = -2 * ll2 + 5 * log(n))
}

#' One-sided fit of the attentive mode
#'
#' Locates the dominant mode \eqn{m} of the score distribution by a
#' Gaussian kernel density estimate, and estimates the attentive scale
#' \eqn{\sigma} from the \emph{left side only} --- the half-width at half maximum
#' of that peak's left flank, deconvolved of the kernel bandwidth, with the median
#' absolute deviation of the points below \eqn{m} as a fallback. The left side of the ensemble score is
#' careless-free by construction at any realistic prevalence, so neither
#' quantity is contaminated when careless respondents are numerous. The careless
#' themselves form a heavy, multi-type right tail with no mode of its own, which
#' is why the attentive mode is the only anchor the score distribution offers
#' (a full two-component fit collapses at high prevalence: its low component
#' absorbs careless observations and both the cut and the prevalence estimate
#' shrink).
#'
#' \eqn{m} is the global maximum of the density. An earlier version took the
#' first local maximum clearing 15\% of the peak, which let a spurious bump in
#' the low tail win over the dominant mode; \eqn{\sigma} was then estimated from
#' the few points below it and collapsed, taking the cut with it. The change is
#' free on real data (identical flags on the validation sample and on six
#' labelled external datasets) and strictly safer on small clean samples.
#'
#' @param eta A numeric vector of ensemble log-odds.
#' @return A list with \code{m} (the attentive mode), \code{sigma} (the
#'   left-side scale), \code{sigma_mad} (the MAD-based fallback, for comparison),
#'   and \code{n_below} (points below the mode).
#' @examples
#' x <- c(rnorm(180, 0), rnorm(120, 4, 2))
#' left_fit(x)
#' @seealso \code{\link{auto_flag}}
#' @export
left_fit <- function(eta) {
  x <- as.numeric(eta); n <- length(x)
  s <- sort(x)
  # 0-based indexing and upper-median conventions mirror the reference JS engine
  q1 <- s[floor(0.25 * n) + 1L]; q3 <- s[floor(0.75 * n) + 1L]
  iqr <- q3 - q1
  sd0 <- sqrt(mean((x - mean(x))^2))
  h <- 0.9 * min(if (sd0 > 0) sd0 else 1, if (iqr > 0) iqr / 1.34 else if (sd0 > 0) sd0 else 1) * n^(-0.2)
  if (!(h > 0)) h <- (if (sd0 > 0) sd0 else 1) * 0.3
  G <- 512L
  lo <- s[1] - h; hi <- s[n] + h
  grid <- seq(lo, hi, length.out = G)
  dens <- vapply(grid, function(g) sum(exp(-0.5 * ((g - x) / h)^2)), numeric(1))
  gm <- which.max(dens)
  dmax <- dens[gm]
  m <- grid[gm]
  ## Scale from the LEFT half-width at half maximum of the attentive peak. The
  ## median absolute deviation of the whole left side is inflated when eta has a
  ## long coherent tail (|z| grows with battery length), which pushed the cut past
  ## everyone and drove the one-sided prevalence negative. The HWHM measures the
  ## shoulder of the peak itself: the far tail lies below half maximum and cannot
  ## influence it, and the careless mass sits on the right, which is never used.
  ## The kernel widens the apparent shoulder, so its bandwidth is deconvolved
  ## (floored at half the raw width, to keep the estimate from collapsing when the
  ## half-width is of the order of the bandwidth). MAD is the fallback for a
  ## degenerate density.
  sig_hw <- 0
  if (gm > 1L) {
    gl <- gm
    while (gl > 1L && dens[gl] > 0.5 * dmax) gl <- gl - 1L
    if (dens[gl] <= 0.5 * dmax) {
      hw <- (m - grid[gl]) / 1.17741
      sig_hw <- sqrt(max(hw * hw - h * h, 0.25 * hw * hw))
    }
  }
  below <- sort(m - x[x < m])
  sig_mad <- if (length(below)) 1.4826 * below[floor(length(below) / 2) + 1L] else 0
  if (!(sig_mad > 0)) sig_mad <- if (sd0 > 0) sd0 else 1e-9
  sig <- if (sig_hw > 0) sig_hw else sig_mad
  list(m = m, sigma = sig, sigma_mad = sig_mad, n_below = length(below))
}

# Resolve a sensitivity preset. Only "standard" ships; the names of the earlier
# two-setting scheme still resolve to it so old scripts keep running, but "high"
# warns, because it used to select a genuinely different (more liberal) cut and a
# silent change of behaviour would be worse than a visible one.
.resolve_sensitivity <- function(s) {
  s <- if (is.character(s) && length(s)) tolower(s[[1]]) else "standard"
  if (identical(s, "high"))
    warning("sensitivity = \"high\" no longer exists: the shipped 2.5-sigma cut is used instead. ",
            "It is the better choice up to a ~35% careless rate; beyond that pass z = 1.5 explicitly.",
            call. = FALSE)
  if (s %in% c("high", "low", "medium")) s <- "standard"
  if (!s %in% names(.RR$sens_z))
    stop("`sensitivity` must be \"standard\"; pass a numeric `z` for a custom cut.", call. = FALSE)
  s
}

#' Automatic, label-free flagging from ensemble scores
#'
#' Turns a vector of scores into per-respondent careless flags without requiring
#' the (unknown) prevalence. Both the cut and the prevalence estimate come from
#' the one-sided fit of the attentive mode (\code{\link{left_fit}}): the cut is
#' that mode's \eqn{z\sigma} upper tail, \eqn{m + z\,\sigma}, and the prevalence
#' is estimated by counting the observations below \eqn{m + \sigma} (a region the
#' careless cannot reach), dividing by \eqn{\Phi(1)} to recover the full attentive
#' mass, and taking the complement. Neither quantity is distorted when careless
#' respondents are numerous, which is where a full two-component fit collapses.
#'
#' The cut is applied \emph{unconditionally}, so a careless-free sample is not
#' guaranteed to yield no flags: it yields its own \eqn{z\sigma} tail, on the
#' order of one to two per cent of respondents. That is a cost we report rather
#' than a promise we make --- see the package vignette for the measured envelope
#' by battery length and sample size.
#'
#' The rule is equivariant to affine transformations of the score, so it can be
#' applied to any index oriented so that high means careless, not only to the
#' ensemble log-odds.
#'
#' @param eta A numeric vector of ensemble log-odds (from \code{\link{ensemble_score}}).
#' @param sensitivity The shipped cut, \code{"standard"} (a \eqn{2.5\sigma} tail on
#'   the attentive mode). It is the only preset: anchoring the cut on the one-sided
#'   fit keeps it the best choice up to a careless rate of about 35\%, so the
#'   second, more liberal level of the earlier two-setting scheme was removed.
#' @param z Optional numeric override for the attentive-tail multiplier (advanced;
#'   ignores \code{sensitivity}). Intended for samples beyond about a third
#'   careless, where a smaller \code{z} (around 1.5) recovers more of them.
#'
#' @return A list with \code{flagged} (logical vector), \code{n_flagged},
#'   \code{flagged_fraction} (the operating-point count as a fraction),
#'   \code{prevalence} (the one-sided estimate of how many careless are present;
#'   reported only when a careless component is detected, otherwise 0),
#'   \code{threshold} on the log-odds scale, \code{two_component} (whether a
#'   careless component was detected), \code{status} --- one of
#'   \code{"flagged"}, \code{"no-subpopulation"} (no careless group found: the
#'   expected result on clean data) or \code{"no-one-past-cut"} (a group was found
#'   but none of its members reaches the cut, the signature of carelessness that
#'   shifts the profile without producing extreme scores) --- the gate \code{separation},
#'   \code{attentive_mode} and \code{attentive_sigma} (the one-sided fit), and
#'   the fitted gate \code{mixture}.
#' @examples
#' d <- simulate_clean(10, 6, 300)
#' eta <- ensemble_score(inject_careless(d$data, prevalence = 0.2)$data)$eta
#' fl <- auto_flag(eta)
#' fl$n_flagged; fl$prevalence
#' @seealso \code{\link{reread}}, \code{\link{left_fit}}, \code{\link{gauss_mixture}}
#' @export
auto_flag <- function(eta, sensitivity = "standard", z = NULL) {
  eta <- as.numeric(eta)
  zz <- if (!is.null(z) && is.numeric(z) && is.finite(z)) z
        else {
          sensitivity <- .resolve_sensitivity(sensitivity)
          .RR$sens_z[[sensitivity]]
        }
  ## The cut is applied unconditionally. A structure gate used to precede it and
  ## flag nobody when a two-component fit found no distinct careless component; it
  ## was removed after measurement showed it left the worst case on clean data
  ## untouched (identical 95th percentile in 42 of 48 design cells) while refusing
  ## to act on a third of low-prevalence samples that did contain careless
  ## respondents. The mixture is still fitted and returned as a diagnostic.
  g <- gauss_mixture(eta)
  pooled <- sqrt((g$v1 + g$v2) / 2); if (!(pooled > 0)) pooled <- 1
  sep <- (g$m2 - g$m1) / pooled
  two <- (g$bic1 - g$bic2 > 2) && g$pi > 0.01 && g$pi < 0.7 && sep > 1.0
  lf <- left_fit(eta)
  cut <- lf$m + zz * lf$sigma
  flagged <- eta > cut
  n_below <- sum(eta < lf$m + lf$sigma)
  prev <- max(0, min(1, 1 - (n_below / stats::pnorm(1)) / length(eta)))
  ## An empty flag list has two different meanings and they must not be conflated:
  ## either no careless subpopulation was detected, or one was detected but none of
  ## its members reaches the cut. The second is possible because the gate and the
  ## cut are estimated from different objects (an EM mixture; the one-sided fit), so
  ## a minority component can sit entirely inside the attentive mode's shoulder.
  status <- if (sum(flagged) == 0) "none-past-cut" else "flagged"
  list(flagged = flagged, n_flagged = sum(flagged),
       flagged_fraction = mean(flagged), prevalence = prev,
       threshold = cut, two_component = two,
       status = status,
       separation = sep, attentive_mode = lf$m, attentive_sigma = lf$sigma,
       mixture = g,
       sensitivity = if (!is.null(z)) z else sensitivity)
}
