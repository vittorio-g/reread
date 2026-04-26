# Phase5b_MultiMetric_tau60.R
# Re-run Phase 5 (multi-metric + ablation) under the NEW operational definition:
#   careless = corruption > 0.60   (Scenario A, the argmax from Phase 6)
#
# This replaces the old tau=0.80 framing. Outputs are saved with "_tau60" suffix
# so the original tau=0.80 files remain available for comparison.

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
  library(randomForest); library(pROC); library(PRROC)
})
select <- dplyr::select

TAU_GT <- 0.60   # <-- the only change versus Phase5_MultiMetric.R

cat("=== Phase 5b: re-evaluation at tau_GT =", TAU_GT, "(Scenario A) ===\n\n")

# Comprehensive metric calculator
compute_all_metrics <- function(probs, gt, threshold = NULL,
                                 fpr_target = NULL, neg_subset = NULL) {
  if (is.null(threshold)) {
    if (is.null(fpr_target)) stop("Need either threshold or fpr_target")
    if (is.null(neg_subset)) neg_subset <- !gt
    threshold <- quantile(probs[neg_subset], 1 - fpr_target)
  }
  flag <- probs > threshold
  TP <- sum( flag &  gt); TN <- sum(!flag & !gt)
  FP <- sum( flag & !gt); FN <- sum(!flag &  gt)
  P <- TP + FN; N_neg <- TN + FP
  sens <- TP / max(1, TP + FN)
  spec <- TN / max(1, TN + FP)
  fpr  <- FP / max(1, FP + TN)
  ppv  <- TP / max(1, TP + FP)
  npv  <- TN / max(1, TN + FN)
  bal_acc <- (sens + spec) / 2
  youden_j <- sens + spec - 1
  g_mean <- sqrt(sens * spec)
  f1 <- 2 * (ppv * sens) / max(1e-12, ppv + sens)
  f2 <- 5 * (ppv * sens) / max(1e-12, 4 * ppv + sens)
  den_mcc <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den_mcc == 0) 0 else (TP*TN - FP*FN) / den_mcc
  total <- length(gt)
  po <- (TP + TN) / total
  pe <- ((TP+FP)*(TP+FN) + (TN+FN)*(TN+FP)) / (total^2)
  kappa <- if (1 - pe < 1e-12) 0 else (po - pe) / (1 - pe)
  auc_val <- tryCatch(as.numeric(auc(roc(gt, probs, quiet = TRUE))),
                      error = function(e) NA)
  pr_obj <- tryCatch(pr.curve(scores.class0 = probs[gt], scores.class1 = probs[!gt],
                              curve = FALSE),
                     error = function(e) NULL)
  auprc <- if (!is.null(pr_obj)) pr_obj$auc.integral else NA
  list(
    threshold = as.numeric(threshold),
    sens = sens, spec = spec, fpr = fpr, ppv = ppv, npv = npv,
    bal_acc = bal_acc, youden_j = youden_j, g_mean = g_mean,
    f1 = f1, f2 = f2, mcc = mcc, kappa = kappa,
    auc = auc_val, auprc = auprc,
    TP = TP, TN = TN, FP = FP, FN = FN
  )
}

# ============================================================
# Training data (Scenario A at tau=0.60): re-train RF + ensemble
# ============================================================
df <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
mask_A <- (df$pattern == "clean") | (df$corruption > TAU_GT)
df_A <- df[mask_A, ]
df_A$gt <- df_A$pattern != "clean" & df_A$corruption > TAU_GT
features <- c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot")

cat(sprintf("Scenario A subset: n=%d (n_pos=%d, n_neg=%d)\n\n",
            nrow(df_A), sum(df_A$gt), sum(!df_A$gt)))

set.seed(42)
df_A$fold <- sample(seq_len(5), nrow(df_A), replace = TRUE)
df_A$prob <- NA_real_
for (f in 1:5) {
  tr <- df_A$fold != f; te <- df_A$fold == f
  X_tr <- as.matrix(df_A[tr, features, drop = FALSE])
  X_te <- as.matrix(df_A[te, features, drop = FALSE])
  y_tr <- factor(as.integer(df_A$gt[tr]), levels = c(0, 1))
  mod <- randomForest(X_tr, y_tr, ntree = 200,
                       mtry = max(1, floor(sqrt(ncol(X_tr)))))
  df_A$prob[te] <- predict(mod, X_te, type = "prob")[, "1"]
}

cat("===========================================================\n")
cat(sprintf(" PHASE 5b — Multi-metric @ tau=%.2f (Scenario A)\n", TAU_GT))
cat("===========================================================\n\n")

# Find optimal thresholds
roc_train <- roc(df_A$gt, df_A$prob, quiet = TRUE)
j_vals <- roc_train$sensitivities + roc_train$specificities - 1
youden_thr <- roc_train$thresholds[which.max(j_vals)]
candidates <- seq(0.01, 0.99, 0.01)
f1_search <- sapply(candidates, function(t) compute_all_metrics(df_A$prob, df_A$gt, threshold = t)$f1)
mcc_search <- sapply(candidates, function(t) compute_all_metrics(df_A$prob, df_A$gt, threshold = t)$mcc)
f1_thr <- candidates[which.max(f1_search)]
mcc_thr <- candidates[which.max(mcc_search)]

results_train <- list(
  "FPR=5%"             = compute_all_metrics(df_A$prob, df_A$gt, fpr_target = 0.05),
  "FPR=10%"            = compute_all_metrics(df_A$prob, df_A$gt, fpr_target = 0.10),
  "Youden's J optimal" = compute_all_metrics(df_A$prob, df_A$gt, threshold = youden_thr),
  "F1 optimal"         = compute_all_metrics(df_A$prob, df_A$gt, threshold = f1_thr),
  "MCC optimal"        = compute_all_metrics(df_A$prob, df_A$gt, threshold = mcc_thr)
)

metric_names <- c("threshold", "sens", "spec", "ppv", "npv",
                  "bal_acc", "youden_j", "g_mean",
                  "f1", "f2", "mcc", "kappa", "auc", "auprc")
summary_train <- do.call(rbind, lapply(names(results_train), function(nm) {
  m <- results_train[[nm]]
  vals <- sapply(metric_names, function(mn) round(m[[mn]], 3))
  data.frame(operating_point = nm, t(vals), stringsAsFactors = FALSE)
}))
print(summary_train, row.names = FALSE)
write.csv(summary_train, "phase5b_train_metrics_tau60.csv", row.names = FALSE)

# ============================================================
# Compare classifiers @ FPR=5%
# ============================================================
cat("\n===========================================================\n")
cat(" Classifier comparison: multi-metric @ FPR=5%\n")
cat("===========================================================\n\n")

library(glmnet)
classifier_results <- list()
for (cls in c("logit", "glmnet", "RF")) {
  preds <- numeric(nrow(df_A))
  for (f in 1:5) {
    tr <- df_A$fold != f; te <- df_A$fold == f
    X_tr <- as.matrix(df_A[tr, features, drop = FALSE])
    X_te <- as.matrix(df_A[te, features, drop = FALSE])
    y_tr <- as.integer(df_A$gt[tr])
    if (cls == "logit") {
      mod <- glm.fit(cbind(1, X_tr), y_tr, family = binomial())
      lin <- cbind(1, X_te) %*% mod$coefficients
      preds[te] <- 1 / (1 + exp(-lin))
    } else if (cls == "glmnet") {
      mod <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5, nfolds = 5)
      preds[te] <- as.numeric(predict(mod, newx = X_te, s = "lambda.min", type = "response"))
    } else if (cls == "RF") {
      y_tr_f <- factor(y_tr, levels = c(0, 1))
      mod <- randomForest(X_tr, y_tr_f, ntree = 200,
                           mtry = max(1, floor(sqrt(ncol(X_tr)))))
      preds[te] <- predict(mod, X_te, type = "prob")[, "1"]
    }
  }
  m <- compute_all_metrics(preds, df_A$gt, fpr_target = 0.05)
  vals <- sapply(metric_names, function(mn) round(m[[mn]], 3))
  classifier_results[[cls]] <- data.frame(classifier = cls, t(vals),
                                            stringsAsFactors = FALSE)
}
classifier_summary <- do.call(rbind, classifier_results)
print(classifier_summary, row.names = FALSE)
write.csv(classifier_summary, "phase5b_classifier_metrics_tau60.csv", row.names = FALSE)

# ============================================================
# Ablation
# ============================================================
cat("\n===========================================================\n")
cat(" Ablation @ tau=0.60\n")
cat("===========================================================\n\n")

ablation_configs <- list(
  "Full ensemble (6)"         = c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot"),
  "Without z_rr_iter_efa (5)" = c("z_rr", "irv", "longstring", "d2", "person_tot"),
  "Without ANY z_rr (4)"       = c("irv", "longstring", "d2", "person_tot"),
  "Only z_rr family (2)"       = c("z_rr", "z_rr_iter_t2"),
  "Triad: iter+IRV+D2 (3)"     = c("z_rr_iter_t2", "irv", "d2"),
  "Aux triad: IRV+D2+LS (3)"   = c("irv", "d2", "longstring"),
  "Aux quad: IRV+D2+LS+PT (4)" = c("irv", "d2", "longstring", "person_tot")
)

ablation_results <- list()
for (cfg_name in names(ablation_configs)) {
  feats_cfg <- ablation_configs[[cfg_name]]
  preds <- numeric(nrow(df_A))
  for (f in 1:5) {
    tr <- df_A$fold != f; te <- df_A$fold == f
    X_tr <- as.matrix(df_A[tr, feats_cfg, drop = FALSE])
    X_te <- as.matrix(df_A[te, feats_cfg, drop = FALSE])
    y_tr <- factor(as.integer(df_A$gt[tr]), levels = c(0, 1))
    mod <- randomForest(X_tr, y_tr, ntree = 200,
                         mtry = max(1, floor(sqrt(ncol(X_tr)))))
    preds[te] <- predict(mod, X_te, type = "prob")[, "1"]
  }
  m <- compute_all_metrics(preds, df_A$gt, fpr_target = 0.05)
  vals <- sapply(metric_names, function(mn) round(m[[mn]], 3))
  ablation_results[[cfg_name]] <- data.frame(
    config = cfg_name, n_features = length(feats_cfg),
    has_RR = any(grepl("z_rr", feats_cfg)),
    t(vals), stringsAsFactors = FALSE)
}
ablation_summary <- do.call(rbind, ablation_results)
ablation_summary <- ablation_summary[order(-ablation_summary$mcc), ]
print(ablation_summary, row.names = FALSE)
write.csv(ablation_summary, "phase5b_ablation_tau60.csv", row.names = FALSE)

full_mcc <- ablation_summary$mcc[ablation_summary$config == "Full ensemble (6)"]
aux_mcc  <- ablation_summary$mcc[ablation_summary$config == "Without ANY z_rr (4)"]
delta_mcc <- full_mcc - aux_mcc
cat(sprintf("\n=== Δ MCC (Full - NoRR) = %.3f - %.3f = %.3f ===\n",
            full_mcc, aux_mcc, delta_mcc))

# Per-pattern: where does ReReReRe help most?
df_A$pred_full <- df_A$prob
preds_no_rr <- numeric(nrow(df_A))
feats_no_rr <- c("irv", "longstring", "d2", "person_tot")
for (f in 1:5) {
  tr <- df_A$fold != f; te <- df_A$fold == f
  X_tr <- as.matrix(df_A[tr, feats_no_rr, drop = FALSE])
  X_te <- as.matrix(df_A[te, feats_no_rr, drop = FALSE])
  y_tr <- factor(as.integer(df_A$gt[tr]), levels = c(0, 1))
  mod <- randomForest(X_tr, y_tr, ntree = 200,
                       mtry = max(1, floor(sqrt(ncol(X_tr)))))
  preds_no_rr[te] <- predict(mod, X_te, type = "prob")[, "1"]
}
df_A$pred_no_rr <- preds_no_rr
thr_full  <- quantile(df_A$pred_full[df_A$pattern == "clean"], 0.95)
thr_no_rr <- quantile(df_A$pred_no_rr[df_A$pattern == "clean"], 0.95)
df_A$flag_full  <- df_A$pred_full  > thr_full
df_A$flag_no_rr <- df_A$pred_no_rr > thr_no_rr

per_pat_ablation <- df_A %>%
  filter(gt) %>%
  group_by(pattern, corruption) %>%
  summarise(
    n = n(),
    rate_full  = round(mean(flag_full), 3),
    rate_no_rr = round(mean(flag_no_rr), 3),
    delta = round(rate_full - rate_no_rr, 3),
    .groups = "drop"
  )
cat("\n=== Per-pattern detection (tau=0.60): full vs no-RR ===\n")
print(as.data.frame(per_pat_ablation), row.names = FALSE)
write.csv(per_pat_ablation, "phase5b_per_pattern_ablation_tau60.csv", row.names = FALSE)

cat("\n=== DONE Phase 5b (tau=0.60) ===\n")
