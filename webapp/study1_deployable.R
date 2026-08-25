## study1_deployable.R — Study 1 ablations scored the way the method is used.
##
## Both Study-1 ablation tables reported an MCC obtained with the labels: tab:ablation
## at a fixed 5% false-positive rate (which needs to know who the true negatives are)
## and tab:ablrpr at the best threshold. Neither is available to an analyst. Here every
## model's leave-one-out score is cut by the SHIPPED label-free calibration
## (reread::auto_flag, m + 2.5*sigma on the one-sided attentive fit),
## the same rule the tool applies to a fresh dataset.
##
## The cut runs on the log-odds, which is the scale the calibration is defined on; the
## LOO probability is a monotone transform of it, so ranking (AUC) is untouched.
##
##   Rscript study1_deployable.R   ->  study1_deployable.csv
.libPaths(c("C:/Users/vitto/AppData/Local/Temp/claude/C--Users-vitto-Desktop-ReReReRe/8f426989-5c28-4009-a838-215f9a1e2ca0/scratchpad/rlib", .libPaths()))
suppressMessages({library(pROC); library(reread); library(careless)})
set.seed(1)

mat0 <- as.matrix(read.csv("_study_matrix.csv", check.names = FALSE))
lab  <- read.csv("_study_labels.csv")[[1]]
it   <- read.csv("_study_items.csv", stringsAsFactors = FALSE)
sc   <- read.csv("_study_scores.csv")
n    <- length(lab)
stopifnot(nrow(mat0) == n, nrow(sc) == n)

impute <- function(m) { for (j in seq_len(ncol(m))) { v <- m[, j]; mm <- median(v, na.rm = TRUE)
  if (is.na(mm)) mm <- 0; m[is.na(v), j] <- mm }; m }
mat <- impute(mat0); matR <- mat
for (j in which(it$reverse == 1)) matR[, j] <- (it$rmax[j] + it$rmin[j]) - mat[, j]
domains <- unique(it$domain[it$domain != "Random"])

## resampled personal reliability, identical to ablation_rpr_study1.R
B <- 200; prr_acc <- matrix(0, n, B)
for (b in seq_len(B)) {
  A <- matrix(NA, n, length(domains)); Bm <- A
  for (di in seq_along(domains)) {
    cols <- which(it$domain == domains[di]); cols <- sample(cols)
    h <- floor(length(cols) / 2); if (h < 1) next
    a <- cols[seq_len(h)]; bb <- cols[(h + 1):length(cols)]
    ma <- rowMeans(matR[, a, drop = FALSE]); mb <- rowMeans(matR[, bb, drop = FALSE])
    A[, di] <- (ma - mean(ma)) / (sd(ma) + 1e-9); Bm[, di] <- (mb - mean(mb)) / (sd(mb) + 1e-9)
  }
  prr_acc[, b] <- vapply(seq_len(n), function(i) {
    x <- A[i, ]; y <- Bm[i, ]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok])) }, numeric(1))
}
prr_v <- rowMeans(prr_acc, na.rm = TRUE)

D <- data.frame(y = as.integer(lab), rr = sc$rr, ls = sc$longstring, pt = sc$person_total,
                irv = sc$irv, d2 = sc$d2, rpr = -prr_v)
D$rpr[is.na(D$rpr)] <- median(D$rpr, na.rm = TRUE)

## Ridge-penalised logistic (lambda = 1), the SAME estimator that produced the shipped
## weights. This matters here in a way it does not for AUC: an unpenalised glm on these
## near-separable data drives the coefficients to the hundreds, so the log-odds saturate
## at +/-30 and the one-sided fit -- which is defined on the scale of a regularised linear
## combiner -- is handed a distribution it was never meant to see. Ranking is unaffected
## either way; the cut is not.
ridge_fit <- function(X, y, lambda = 1) {
  P <- ncol(X); b <- rep(0, P)
  for (it in seq_len(60)) {
    e <- as.vector(X %*% b); pr <- 1 / (1 + exp(-e)); w <- pmax(pr * (1 - pr), 1e-6)
    g <- as.vector(t(X) %*% (y - pr)); H <- t(X * w) %*% X
    pen <- rep(lambda, P); pen[1] <- 0                      # intercept unpenalised
    g <- g - pen * b; diag(H) <- diag(H) + pen
    s <- tryCatch(solve(H, g), error = function(e) rep(0, P))
    b <- b + s; if (max(abs(s)) < 1e-9) break
  }
  b
}
loo <- function(form) {
  X <- model.matrix(form, D); y <- D$y
  eta <- numeric(n)
  for (i in seq_len(n)) {
    b <- ridge_fit(X[-i, , drop = FALSE], y[-i])
    eta[i] <- sum(X[i, ] * b)
  }
  eta                                                        # LOO linear predictor (log-odds)
}
mccf <- function(f, y) { tp <- sum(f & y == 1); tn <- sum(!f & y == 0); fp <- sum(f & y == 0); fn <- sum(!f & y == 1)
  d <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)); if (d == 0) 0 else (tp * tn - fp * fn) / d }
aucd <- function(y, s) as.numeric(pROC::auc(pROC::roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE)))
ciauc <- function(y, s) { r <- roc(y, s, direction = "<", quiet = TRUE); a <- ci.auc(r, method = "delong")
  c(as.numeric(auc(r)), a[1], a[3]) }

## the deployable cut: the shipped calibration on the LOO log-odds
autocut <- function(eta) {
  af <- reread::auto_flag(eta)
  list(flag = af$flagged, gate = af$two_component, nflag = af$n_flagged)
}
## bootstrap CI on the deployable MCC: the cut is re-derived inside each resample,
## because in use it would be re-derived on whatever sample the analyst has
boot_ci <- function(eta, y, Bn = 1000) {
  v <- numeric(Bn)
  for (b in seq_len(Bn)) {
    ix <- sample(n, n, replace = TRUE)
    if (length(unique(y[ix])) < 2) { v[b] <- NA; next }
    ac <- autocut(eta[ix]); v[b] <- mccf(ac$flag, y[ix])
  }
  quantile(v, c(.025, .975), na.rm = TRUE)
}

MODELS <- list(
  "Shipped ensemble (rc + LongString + PT)"      = y ~ rr + ls + pt,
  "refit adding IRV + D2 (5 features)"           = y ~ rr + ls + pt + irv + d2,
  "Auxiliaries only, no rc (IRV + LS + D2 + PT)" = y ~ irv + ls + d2 + pt,
  "rc alone"                                     = y ~ rr,
  "RPR alone"                                    = y ~ rpr,
  "LongString + PT"                              = y ~ ls + pt,
  "RPR + LongString + PT"                        = y ~ rpr + ls + pt,
  "rc + RPR + LongString + PT"                   = y ~ rr + rpr + ls + pt
)

## The induced sample is 55% careless by design, far above any deployment condition and
## above the range the shipped cut is calibrated for. Reporting the automatic cut on it as
## it stands would say more about the base rate than about the models. So the deployable
## MCC is read where the method is meant to run: all 70 attentive respondents are kept and
## real careless respondents are sampled in to a target rate, the cut is re-derived on each
## mixture, and the result is averaged over the realistic 10-25% band. The leave-one-out
## scores are reused unchanged, so no respondent is scored by a model fitted on them.
RATES <- c(.10, .15, .20, .25); REPS <- 200
careful_ix <- which(D$y == 0); careless_ix <- which(D$y == 1)
band <- function(eta) {
  res <- list()
  for (rt in RATES) {
    kk <- round(length(careful_ix) * rt / (1 - rt))
    mm <- numeric(REPS); ff <- numeric(REPS)
    for (b in seq_len(REPS)) {
      ix <- c(careful_ix, sample(careless_ix, kk))
      ac <- autocut(eta[ix]); mm[b] <- mccf(ac$flag, D$y[ix]); ff[b] <- ac$nflag / length(ix)
    }
    res[[as.character(rt)]] <- c(mcc = mean(mm), flag = mean(ff))
  }
  c(mcc = mean(sapply(res, function(v) v["mcc"])),
    flag = mean(sapply(res, function(v) v["flag"])),
    lo = min(sapply(res, function(v) v["mcc"])), hi = max(sapply(res, function(v) v["mcc"])))
}

## Reference row: the pipeline exactly as deployed -- the FROZEN shipped weights, not a
## refit. Every other row refits its combiner leave-one-out, which changes the geometry of
## the score and therefore where the cut lands; this row shows what the released tool does
## on the same mixtures, so the two are not confused.
frozen_eta <- sc$eta
out <- data.frame()
cat(sprintf("=== Study 1 (N=%d), leave-one-out scores, cut by the shipped label-free calibration ===\n", n))
cat("                                                                   full sample (55%% careless) |  realistic 10-25%% band\n")
cat("model                                          AUC [95% CI]         flagged  MCC   |  MCC [min-max]   %flagged\n")
for (k in names(MODELS)) {
  eta <- loo(MODELS[[k]])
  a <- ciauc(D$y, eta)
  ac <- autocut(eta); m <- mccf(ac$flag, D$y)
  bd <- band(eta)
  cat(sprintf("%-46s %.3f [%.3f, %.3f]  %3d/%d  %.3f |  %.3f [%.2f-%.2f]   %.0f%%\n",
              k, a[1], a[2], a[3], ac$nflag, n, m, bd["mcc"], bd["lo"], bd["hi"], 100 * bd["flag"]))
  out <- rbind(out, data.frame(model = k, auc = round(a[1], 4), auc_lo = round(a[2], 4), auc_hi = round(a[3], 4),
    full_gate_open = ac$gate, full_n_flagged = ac$nflag, full_mcc_auto = round(m, 4),
    band_mcc_auto = round(bd["mcc"], 4), band_mcc_min = round(bd["lo"], 4),
    band_mcc_max = round(bd["hi"], 4), band_flag_share = round(bd["flag"], 4)))
}
## frozen-weight reference
{
  a <- ciauc(D$y, frozen_eta); ac <- autocut(frozen_eta); m <- mccf(ac$flag, D$y); bd <- band(frozen_eta)
  cat(sprintf("%-46s %.3f [%.3f, %.3f]  %3d/%d  %.3f |  %.3f [%.2f-%.2f]   %.0f%%\n",
              "[reference] shipped weights, frozen", a[1], a[2], a[3], ac$nflag, n, m,
              bd["mcc"], bd["lo"], bd["hi"], 100 * bd["flag"]))
  out <- rbind(out, data.frame(model = "[reference] shipped weights, frozen",
    auc = round(a[1], 4), auc_lo = round(a[2], 4), auc_hi = round(a[3], 4),
    full_gate_open = ac$gate, full_n_flagged = ac$nflag, full_mcc_auto = round(m, 4),
    band_mcc_auto = round(bd["mcc"], 4), band_mcc_min = round(bd["lo"], 4),
    band_mcc_max = round(bd["hi"], 4), band_flag_share = round(bd["flag"], 4)))
}
write.csv(out, "study1_deployable.csv", row.names = FALSE)
cat("\nwrote study1_deployable.csv\n")
