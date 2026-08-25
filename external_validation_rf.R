# External validation of the FULL RF ensemble on real datasets.
#
# For each dataset with ground truth, we compute the 6 standard features
# (z_rr, z_rr_iter, irv, longstring, d2, person_tot), then train a Random
# Forest with 5-fold CV. We compare:
#   - Full ensemble (6 features)
#   - Triad (z_rr_iter + irv + d2)
#   - No-RR (4 auxiliaries only)
#   - Single z_rr_iter
#   - Mahalanobis chi-square (.001) practical baseline
#
# Operating point: 95th percentile of negative-class predictions (oracle on
# REAL data this is fine because we have labels).

suppressPackageStartupMessages({
  library(dplyr); library(MASS); library(randomForest)
  library(pROC); library(haven)
})
source("Synthetic_Good_Responses_2.R")  # not strictly needed, but harmless
source("ReReReRe.R")
source("ReReReRe_IterCoupled.R")

# ---- helpers ----
compute_irv <- function(d) apply(d, 1, sd, na.rm = TRUE)
compute_longstring <- function(d) apply(d, 1, function(r) {
  r <- r[!is.na(r)]; if (!length(r)) 0 else max(rle(r)$lengths)
})
compute_d2 <- function(d) {
  m <- as.matrix(d)
  m <- m[, apply(m, 2, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
  if (ncol(m) < 2) return(rep(0, nrow(d)))
  cov_mat <- cov(m, use = "pairwise.complete.obs")
  center <- colMeans(m, na.rm = TRUE)
  cov_inv <- tryCatch(ginv(cov_mat),
                      error = function(e) solve(cov_mat + diag(1e-6, ncol(cov_mat))))
  apply(m, 1, function(x) {
    dlt <- x - center
    as.numeric(t(dlt) %*% cov_inv %*% dlt)
  })
}
compute_person_total <- function(d) {
  m <- as.matrix(d); im <- colMeans(m, na.rm = TRUE)
  apply(m, 1, function(r) {
    if (sd(r, na.rm = TRUE) == 0 || sd(im) == 0) 0
    else suppressWarnings(cor(r, im, use = "pairwise.complete.obs"))
  })
}

mcc_calc <- function(tp, tn, fp, fn) {
  tp <- as.double(tp); tn <- as.double(tn)
  fp <- as.double(fp); fn <- as.double(fn)
  num <- (tp * tn) - (fp * fn)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (is.na(den) || den == 0) 0 else num / den
}
metrics_at <- function(prob, y, tau) {
  pred <- prob >= tau
  TP <- sum(pred & y == 1); FP <- sum(pred & y == 0)
  TN <- sum(!pred & y == 0); FN <- sum(!pred & y == 1)
  sens <- if (TP + FN > 0) TP / (TP + FN) else NA
  spec <- if (TN + FP > 0) TN / (TN + FP) else NA
  ppv  <- if (TP + FP > 0) TP / (TP + FP) else NA
  f1   <- if (!is.na(ppv) && !is.na(sens) && (ppv + sens) > 0) 2 * ppv * sens / (ppv + sens) else NA
  list(mcc = mcc_calc(TP, TN, FP, FN), f1 = f1, sens = sens, spec = spec)
}
cv_rf_predict <- function(X, y, K = 5, seed = 42) {
  set.seed(seed)
  n <- length(y); folds <- sample(rep(1:K, length.out = n))
  oof <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k; te <- folds == k
    rf <- randomForest(x = X[tr, , drop = FALSE],
                       y = factor(y[tr], levels = c(0, 1)),
                       ntree = 200, mtry = max(1, floor(sqrt(ncol(X)))))
    pp <- predict(rf, X[te, , drop = FALSE], type = "prob")
    oof[te] <- pp[, "1"]
  }
  oof
}

# ---- compute features for a single dataset ----
compute_features <- function(items_df) {
  items_df <- as.data.frame(lapply(as.data.frame(items_df), as.numeric))
  items_df <- items_df[, apply(items_df, 2, function(x) sd(x, na.rm = TRUE) > 0)]
  data_mat <- as.matrix(items_df)
  rr <- tryCatch(ReReReRe(data_mat, corProp = 0.03, iterations = 100,
                          align_signs = TRUE, mode = "auto",
                          auto_z = FALSE, variance_penalty = FALSE),
                 error = function(e) NULL)
  z_rr <- if (is.null(rr)) rep(0, nrow(data_mat)) else rr$z_score
  z_rr_iter <- tryCatch({
    rr_it <- score_iter_coupled(data_mat, corProp = 0.03, iterations = 100,
                                align_signs = TRUE, trim_pct = 0.20)
    rr_it$z_score
  }, error = function(e) rep(0, nrow(data_mat)))
  data.frame(
    z_rr       = -z_rr,
    z_rr_iter  = -z_rr_iter,
    irv        = -compute_irv(data_mat),
    longstring =  compute_longstring(data_mat),
    d2         =  compute_d2(data_mat),
    person_tot = -compute_person_total(data_mat)
  )
}

# ---- run RF on one set of features and compute metrics ----
run_one <- function(feats, y, name, seed = 42) {
  X <- as.matrix(feats)
  prob <- cv_rf_predict(X, y, K = 5, seed = seed)
  # operating point: 95th percentile of negative-class predictions
  tau_op <- quantile(prob[y == 0], 0.95, na.rm = TRUE)
  m <- metrics_at(prob, y, tau_op)
  auc_v <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE,
                                                     direction = "<"))),
                    error = function(e) NA)
  data.frame(method = name, n_features = ncol(X),
             mcc = m$mcc, f1 = m$f1, sens = m$sens, spec = m$spec,
             auc = auc_v, tau_op = as.numeric(tau_op))
}

# ---- bootstrap CI for delta-MCC (full vs no-rr) ----
boot_delta <- function(prob_full, prob_norr, y, n_boot = 500, seed = 42) {
  set.seed(seed)
  thr_f <- quantile(prob_full[y == 0], 0.95, na.rm = TRUE)
  thr_n <- quantile(prob_norr[y == 0], 0.95, na.rm = TRUE)
  pred_f <- as.integer(prob_full >= thr_f)
  pred_n <- as.integer(prob_norr >= thr_n)
  n <- length(y); deltas <- numeric(0)
  for (b in 1:n_boot) {
    idx <- sample.int(n, n, replace = TRUE)
    yb <- y[idx]; pf <- pred_f[idx]; pn <- pred_n[idx]
    if (length(unique(yb)) < 2) next
    tp_f <- sum(pf == 1 & yb == 1); fp_f <- sum(pf == 1 & yb == 0)
    tn_f <- sum(pf == 0 & yb == 0); fn_f <- sum(pf == 0 & yb == 1)
    tp_n <- sum(pn == 1 & yb == 1); fp_n <- sum(pn == 1 & yb == 0)
    tn_n <- sum(pn == 0 & yb == 0); fn_n <- sum(pn == 0 & yb == 1)
    deltas <- c(deltas, mcc_calc(tp_f, tn_f, fp_f, fn_f) - mcc_calc(tp_n, tn_n, fp_n, fn_n))
  }
  c(lo = quantile(deltas, 0.025, names = FALSE),
    hi = quantile(deltas, 0.975, names = FALSE),
    mean = mean(deltas),
    p_gt0 = mean(deltas > 0))
}

# ---- pipeline for one dataset ----
process_dataset <- function(name, items, labels) {
  cat(sprintf("\n=== %s (n=%d items=%d careless=%d/%.1f%%) ===\n",
              name, nrow(items), ncol(items), sum(labels), mean(labels) * 100))
  feats <- compute_features(items)
  if (any(is.na(feats))) {
    cat("  WARN: NAs in features; replacing with 0\n")
    feats[is.na(feats)] <- 0
  }
  FEATS_FULL  <- c("z_rr", "z_rr_iter", "irv", "longstring", "d2", "person_tot")
  FEATS_TRIAD <- c("z_rr_iter", "irv", "d2")
  FEATS_NORR  <- c("irv", "longstring", "d2", "person_tot")
  FEATS_ZRR   <- c("z_rr_iter")

  rows <- rbind(
    cbind(dataset = name, run_one(feats[, FEATS_FULL,  drop = FALSE], labels, "Full_6")),
    cbind(dataset = name, run_one(feats[, FEATS_TRIAD, drop = FALSE], labels, "Triad_3")),
    cbind(dataset = name, run_one(feats[, FEATS_NORR,  drop = FALSE], labels, "NoRR_4")),
    cbind(dataset = name, run_one(feats[, FEATS_ZRR,   drop = FALSE], labels, "z_rr_iter_only"))
  )
  # Mahalanobis chi-sq baseline
  d2 <- feats$d2; p_items <- ncol(items)
  thr_chi <- qchisq(1 - 0.001, df = p_items)
  pred_mah <- as.integer(d2 > thr_chi)
  m_mah <- metrics_at(d2, labels, thr_chi)
  auc_mah <- tryCatch(as.numeric(pROC::auc(pROC::roc(labels, d2, quiet = TRUE,
                                                       direction = "<"))),
                      error = function(e) NA)
  rows <- rbind(rows, data.frame(
    dataset = name, method = "Mah_chi2_001", n_features = 1,
    mcc = m_mah$mcc, f1 = m_mah$f1, sens = m_mah$sens, spec = m_mah$spec,
    auc = auc_mah, tau_op = as.numeric(thr_chi)
  ))
  print(rows[, c("method", "n_features", "mcc", "f1", "sens", "spec", "auc")])

  # bootstrap CI on Full vs NoRR delta
  X_full <- as.matrix(feats[, FEATS_FULL])
  X_norr <- as.matrix(feats[, FEATS_NORR])
  prob_full <- cv_rf_predict(X_full, labels, K = 5, seed = 42)
  prob_norr <- cv_rf_predict(X_norr, labels, K = 5, seed = 42)
  bci <- boot_delta(prob_full, prob_norr, labels)
  cat(sprintf("  Δ MCC bootstrap: mean=%.3f [95%% CI %.3f to %.3f] P(Δ>0)=%.3f\n",
              bci["mean"], bci["lo"], bci["hi"], bci["p_gt0"]))

  list(rows = rows, ci = data.frame(dataset = name, delta_mean = bci["mean"],
                                     delta_lo = bci["lo"], delta_hi = bci["hi"],
                                     p_gt0 = bci["p_gt0"]))
}

cat("\n====================================================\n")
cat(" External validation of FULL RF ensemble (revision)\n")
cat("====================================================\n")

all_rows <- list(); all_ci <- list()

# Schroeders 2022
d1 <- read.csv("external_datasets/01_Schroeders_2022/data_mod_resp.csv", sep = ";")
d1_label <- d1$Careless
d1_items <- dplyr::select(d1, dplyr::starts_with("HE"))
ok <- complete.cases(d1_items) & !is.na(d1_label)
res <- process_dataset("Schroeders2022", d1_items[ok, ], d1_label[ok])
all_rows[[length(all_rows) + 1]] <- res$rows; all_ci[[length(all_ci) + 1]] <- res$ci

# Schneider QoL
d2 <- read.csv("external_datasets/06_Schneider_QoL/carersp.csv")
d2_label <- d2$c01
item_cols <- grep("^(dep|pain|cog|fat)\\d", names(d2), value = TRUE)
d2_items <- d2[, item_cols]
ok <- complete.cases(d2_items) & !is.na(d2_label)
res <- process_dataset("Schneider_QoL", d2_items[ok, ], d2_label[ok])
all_rows[[length(all_rows) + 1]] <- res$rows; all_ci[[length(all_ci) + 1]] <- res$ci

# Niessen 2016
d3 <- read_sav("external_datasets/08_Niessen_2016/Raw_data.sav")
d3 <- dplyr::filter(d3, Use_me == 1)
d3_label <- as.integer(d3$Conditon)
item_cols <- grep("^[ECANO]\\d+$", names(d3), value = TRUE)
d3_items <- as.data.frame(d3[, item_cols])
ok <- complete.cases(d3_items) & !is.na(d3_label)
res <- process_dataset("Niessen2016", d3_items[ok, ], d3_label[ok])
all_rows[[length(all_rows) + 1]] <- res$rows; all_ci[[length(all_ci) + 1]] <- res$ci

# Goldammer S1 (BFI-2 bidir, 60 items)
g1 <- read.csv("external_datasets/13_Goldammer_2024/Study_1.csv")
labels_s1 <- as.integer(g1$condition > 0)
# BFI-2 items: typically 60 items named bfi_1 .. bfi_60 (verify)
# Goldammer items: e_gs1, a_mt1, c_ol1, n_an1, o_ae1 etc. -- 5 traits x scales x items
# Use only non-reverse versions (no R suffix) to test align_signs
g1_item_cols <- grep("^[eacno]_[a-z]+[0-9]+$", names(g1), value = TRUE)
g1_items <- g1[, g1_item_cols]
ok <- complete.cases(g1_items) & !is.na(labels_s1)
g1_items <- g1_items[ok, ]; labels_s1 <- labels_s1[ok]
cat(sprintf("  Goldammer S1 columns kept: %d\n", ncol(g1_items)))
if (ncol(g1_items) >= 30 && length(unique(labels_s1)) == 2) {
  res <- process_dataset("Goldammer_S1", g1_items, labels_s1)
  all_rows[[length(all_rows) + 1]] <- res$rows; all_ci[[length(all_ci) + 1]] <- res$ci
} else {
  cat("  Skipping Goldammer S1 (couldn't identify items)\n")
}

# Goldammer S2 (IPIP unidir, 60 items)
g2 <- read.csv("external_datasets/13_Goldammer_2024/Study_2.csv")
labels_s2 <- as.integer(g2$condition > 0)
g2_item_cols <- grep("^[eacno]_[a-z]+[0-9]+$", names(g2), value = TRUE)
g2_items <- g2[, g2_item_cols]
ok <- complete.cases(g2_items) & !is.na(labels_s2)
g2_items <- g2_items[ok, ]; labels_s2 <- labels_s2[ok]
cat(sprintf("  Goldammer S2 columns kept: %d\n", ncol(g2_items)))
if (ncol(g2_items) >= 30 && length(unique(labels_s2)) == 2) {
  res <- process_dataset("Goldammer_S2", g2_items, labels_s2)
  all_rows[[length(all_rows) + 1]] <- res$rows; all_ci[[length(all_ci) + 1]] <- res$ci
} else {
  cat("  Skipping Goldammer S2 (couldn't identify items)\n")
}

# ---- write outputs ----
out_rows <- do.call(rbind, all_rows)
write.csv(out_rows, "external_rf_ensemble.csv", row.names = FALSE)
cat("\nSaved external_rf_ensemble.csv\n")

out_ci <- do.call(rbind, all_ci)
write.csv(out_ci, "external_rf_bootstrap.csv", row.names = FALSE)
cat("Saved external_rf_bootstrap.csv\n")

cat("\nFinal table:\n")
print(out_rows)
cat("\nBootstrap CIs (Full - NoRR):\n")
print(out_ci)
