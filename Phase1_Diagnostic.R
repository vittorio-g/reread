# Phase1_Diagnostic.R
# Investigate whether MCC > 0.7 at GT>80% is achievable.
# Two GT scenarios + multiple classifiers + multiple operating points.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(MASS)
})

# Try to load extra classifiers
have_rf  <- requireNamespace("randomForest", quietly = TRUE)
have_glmnet <- requireNamespace("glmnet", quietly = TRUE)
if (have_rf)  library(randomForest)
if (have_glmnet) library(glmnet)
cat(sprintf("randomForest: %s | glmnet: %s\n", have_rf, have_glmnet))

df <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d respondents\n", nrow(df)))

# Build the two GT scenarios
# Scenario A — clean vs full-careless (>80%) only
# Scenario B — clean + middle (10-80%) as negative, full-careless (>80%) as positive
df$gt_full <- (df$pattern != "clean") & (df$corruption > 0.80)

# Best iter_EFA from previous tuning
features_full <- c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot")
features_triad <- c("z_rr", "z_rr_iter_t2", "irv", "d2")

# Helper: compute metrics at fixed FPR
compute_metrics <- function(probs, gt, fpr_target, neg_subset = NULL) {
  # Calibrate threshold using only specified negatives
  if (is.null(neg_subset)) neg_subset <- !gt
  thr <- quantile(probs[neg_subset], 1 - fpr_target)
  flag <- probs > thr
  TP <- sum(flag & gt); TN <- sum(!flag & !gt)
  FP <- sum(flag & !gt); FN <- sum(!flag & gt)
  den <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  list(
    sens = TP / max(1, TP + FN),
    spec = TN / max(1, TN + FP),
    fpr  = FP / max(1, FP + TN),
    mcc  = if (den == 0) 0 else (TP*TN - FP*FN) / den,
    threshold = thr
  )
}

# 5-fold CV with multiple classifiers
set.seed(42)
folds <- sample(seq_len(5), nrow(df), replace = TRUE)

cv_predict <- function(df, feats, gt_col, classifier = "logit", folds) {
  preds <- numeric(nrow(df))
  for (f in unique(folds)) {
    tr <- folds != f; te <- folds == f
    X_tr <- as.matrix(df[tr, feats, drop = FALSE])
    X_te <- as.matrix(df[te, feats, drop = FALSE])
    y_tr <- as.integer(df[[gt_col]][tr])
    if (classifier == "logit") {
      mod <- tryCatch(glm.fit(cbind(1, X_tr), y_tr, family = binomial()),
                      error = function(e) NULL)
      if (is.null(mod)) { preds[te] <- 0.5; next }
      lin <- cbind(1, X_te) %*% mod$coefficients
      preds[te] <- 1 / (1 + exp(-lin))
    } else if (classifier == "rf" && have_rf) {
      mod <- randomForest(X_tr, factor(y_tr, levels = c(0,1)),
                           ntree = 200, mtry = max(1, floor(sqrt(ncol(X_tr)))))
      preds[te] <- predict(mod, X_te, type = "prob")[, "1"]
    } else if (classifier == "glmnet" && have_glmnet) {
      mod <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5, nfolds = 5)
      preds[te] <- as.numeric(predict(mod, newx = X_te, s = "lambda.min",
                                        type = "response"))
    } else {
      preds[te] <- 0.5
    }
  }
  preds
}

# ============================================================
# Scenario A: clean vs corruption>80% only (exclude middle 10-80%)
# ============================================================
cat("\n===========================================================\n")
cat(" SCENARIO A: clean vs full-careless (corruption>80%) only\n")
cat(" Excludes the 'fuzzy middle' (10-80% corruption respondents)\n")
cat("===========================================================\n")

mask_A <- (df$pattern == "clean") | (df$corruption > 0.80)
df_A <- df[mask_A, ]
folds_A <- sample(seq_len(5), nrow(df_A), replace = TRUE)

cat(sprintf("N = %d (clean=%d, careless=%d)\n",
            nrow(df_A), sum(df_A$pattern == "clean"),
            sum(df_A$pattern != "clean")))

results_A <- data.frame()
for (cls in c("logit", if (have_glmnet) "glmnet", if (have_rf) "rf")) {
  for (fpr_t in c(0.05, 0.10, 0.20)) {
    # Triad
    probs_t <- cv_predict(df_A, features_triad, "gt_full", cls, folds_A)
    m_t <- compute_metrics(probs_t, df_A$gt_full, fpr_t)
    # Full ensemble
    probs_e <- cv_predict(df_A, features_full, "gt_full", cls, folds_A)
    m_e <- compute_metrics(probs_e, df_A$gt_full, fpr_t)
    results_A <- rbind(results_A, data.frame(
      classifier = cls, fpr_target = fpr_t, ensemble = "triad+iter (4)",
      sens = round(m_t$sens, 3), fpr = round(m_t$fpr, 3),
      mcc  = round(m_t$mcc,  3)))
    results_A <- rbind(results_A, data.frame(
      classifier = cls, fpr_target = fpr_t, ensemble = "full (6)",
      sens = round(m_e$sens, 3), fpr = round(m_e$fpr, 3),
      mcc  = round(m_e$mcc,  3)))
  }
}
print(as.data.frame(results_A %>% arrange(desc(mcc))), row.names = FALSE)
write.csv(results_A, "phase1_scenarioA.csv", row.names = FALSE)

# ============================================================
# Scenario B: full data, GT>80%, negative = clean + middle (10-80%)
# ============================================================
cat("\n===========================================================\n")
cat(" SCENARIO B: full data, GT>80%, all corruption ≤80% as negative\n")
cat("===========================================================\n")

cat(sprintf("N = %d (positives=%d, negatives=%d)\n",
            nrow(df), sum(df$gt_full), sum(!df$gt_full)))

results_B <- data.frame()
for (cls in c("logit", if (have_glmnet) "glmnet", if (have_rf) "rf")) {
  for (fpr_t in c(0.05, 0.10, 0.20)) {
    # FPR computed only on clean (most realistic — we don't want to flag clean)
    probs_t <- cv_predict(df, features_triad, "gt_full", cls, folds)
    m_t <- compute_metrics(probs_t, df$gt_full, fpr_t,
                            neg_subset = (df$pattern == "clean"))
    probs_e <- cv_predict(df, features_full, "gt_full", cls, folds)
    m_e <- compute_metrics(probs_e, df$gt_full, fpr_t,
                            neg_subset = (df$pattern == "clean"))
    results_B <- rbind(results_B, data.frame(
      classifier = cls, fpr_target = fpr_t, ensemble = "triad+iter (4)",
      sens = round(m_t$sens, 3), fpr_overall = round(m_t$fpr, 3),
      mcc  = round(m_t$mcc,  3)))
    results_B <- rbind(results_B, data.frame(
      classifier = cls, fpr_target = fpr_t, ensemble = "full (6)",
      sens = round(m_e$sens, 3), fpr_overall = round(m_e$fpr, 3),
      mcc  = round(m_e$mcc,  3)))
  }
}
print(as.data.frame(results_B %>% arrange(desc(mcc))), row.names = FALSE)
write.csv(results_B, "phase1_scenarioB.csv", row.names = FALSE)

# Diagnostic: best model in each scenario
cat("\n===========================================================\n")
cat(" SUMMARY: best models\n")
cat("===========================================================\n")
best_A <- results_A %>% arrange(desc(mcc)) %>% head(1)
best_B <- results_B %>% arrange(desc(mcc)) %>% head(1)
cat(sprintf("Best Scenario A: %s + %s @ FPR=%.2f → MCC=%.3f (sens=%.3f)\n",
            best_A$classifier, best_A$ensemble, best_A$fpr_target,
            best_A$mcc, best_A$sens))
cat(sprintf("Best Scenario B: %s + %s @ FPR=%.2f → MCC=%.3f (sens=%.3f)\n",
            best_B$classifier, best_B$ensemble, best_B$fpr_target,
            best_B$mcc, best_B$sens))

cat("\n=== DONE Phase 1 ===\n")
