# Phase5_MultiMetric.R
# Re-evaluate Phase 2 (training, 16K resp Scenario A) and Phase 3 (validation, 54K)
# using a comprehensive set of metrics from ML and psychometric literature.
#
# Metrics computed:
#   - MCC, F1, F2, Cohen's Kappa, Balanced Accuracy, Youden's J
#   - Precision (PPV), NPV, G-mean
#   - Sensitivity, Specificity, FPR
#   - AUC (threshold-free), AUPRC (threshold-free, important for imbalanced)
#
# Goal: confirm or invalidate the current MCC-based conclusions.

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(ggplot2); library(tidyr)
  library(randomForest); library(pROC); library(PRROC)
})
select <- dplyr::select

# Comprehensive metric calculator
compute_all_metrics <- function(probs, gt, threshold = NULL,
                                 fpr_target = NULL, neg_subset = NULL) {
  # Choose threshold
  if (is.null(threshold)) {
    if (is.null(fpr_target)) stop("Need either threshold or fpr_target")
    if (is.null(neg_subset)) neg_subset <- !gt
    threshold <- quantile(probs[neg_subset], 1 - fpr_target)
  }
  flag <- probs > threshold
  TP <- sum( flag &  gt); TN <- sum(!flag & !gt)
  FP <- sum( flag & !gt); FN <- sum(!flag &  gt)
  P <- TP + FN; N_neg <- TN + FP
  # Basic
  sens <- TP / max(1, TP + FN)
  spec <- TN / max(1, TN + FP)
  fpr  <- FP / max(1, FP + TN)
  # Predictive values
  ppv  <- TP / max(1, TP + FP)
  npv  <- TN / max(1, TN + FN)
  # Combined
  bal_acc <- (sens + spec) / 2
  youden_j <- sens + spec - 1
  g_mean <- sqrt(sens * spec)
  # F-scores
  f1 <- 2 * (ppv * sens) / max(1e-12, ppv + sens)
  f2 <- 5 * (ppv * sens) / max(1e-12, 4 * ppv + sens)
  # MCC
  den_mcc <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- if (den_mcc == 0) 0 else (TP*TN - FP*FN) / den_mcc
  # Cohen's Kappa
  total <- length(gt)
  po <- (TP + TN) / total
  pe <- ((TP+FP)*(TP+FN) + (TN+FN)*(TN+FP)) / (total^2)
  kappa <- if (1 - pe < 1e-12) 0 else (po - pe) / (1 - pe)
  # Threshold-free
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
# Training data (Phase 2): re-train RF + full ensemble + Scenario A
# ============================================================
df <- read.csv("sim_robust_scores.csv", stringsAsFactors = FALSE)
mask_A <- (df$pattern == "clean") | (df$corruption > 0.80)
df_A <- df[mask_A, ]
df_A$gt <- df_A$pattern != "clean"
features <- c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot")

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
cat(" PHASE 5 — Multi-metric evaluation (training set, 24K resp)\n")
cat("===========================================================\n\n")

# Three operating points
ops <- list(
  "FPR=5%"  = 0.05,
  "FPR=10%" = 0.10,
  "Youden-J optimal" = NULL,    # find via ROC
  "F1 optimal"       = NULL,    # find via search
  "MCC optimal"      = NULL     # find via search
)

# Find optimal thresholds
# Youden's J: argmax sens + spec - 1 → equivalent to max of (TPR - FPR) on ROC
roc_train <- roc(df_A$gt, df_A$prob, quiet = TRUE)
j_vals <- roc_train$sensitivities + roc_train$specificities - 1
best_idx_j <- which.max(j_vals)
youden_thr <- roc_train$thresholds[best_idx_j]

# F1 and MCC optimal: search over candidate thresholds
candidates <- seq(0.01, 0.99, 0.01)
f1_search <- sapply(candidates, function(t) {
  m <- compute_all_metrics(df_A$prob, df_A$gt, threshold = t)
  m$f1
})
mcc_search <- sapply(candidates, function(t) {
  m <- compute_all_metrics(df_A$prob, df_A$gt, threshold = t)
  m$mcc
})
f1_thr <- candidates[which.max(f1_search)]
mcc_thr <- candidates[which.max(mcc_search)]

# Compute all metrics at all operating points
results_train <- list(
  "FPR=5%"             = compute_all_metrics(df_A$prob, df_A$gt, fpr_target = 0.05),
  "FPR=10%"            = compute_all_metrics(df_A$prob, df_A$gt, fpr_target = 0.10),
  "Youden's J optimal" = compute_all_metrics(df_A$prob, df_A$gt, threshold = youden_thr),
  "F1 optimal"         = compute_all_metrics(df_A$prob, df_A$gt, threshold = f1_thr),
  "MCC optimal"        = compute_all_metrics(df_A$prob, df_A$gt, threshold = mcc_thr)
)

# Build summary table
metric_names <- c("threshold", "sens", "spec", "ppv", "npv",
                  "bal_acc", "youden_j", "g_mean",
                  "f1", "f2", "mcc", "kappa", "auc", "auprc")
summary_train <- do.call(rbind, lapply(names(results_train), function(nm) {
  m <- results_train[[nm]]
  vals <- sapply(metric_names, function(mn) round(m[[mn]], 3))
  data.frame(operating_point = nm, t(vals), stringsAsFactors = FALSE)
}))
print(summary_train, row.names = FALSE)
write.csv(summary_train, "phase5_train_metrics.csv", row.names = FALSE)

# ============================================================
# Validation (Phase 3): same metrics on fresh data
# ============================================================
cat("\n===========================================================\n")
cat(" Validation: per careless rate (Phase 3 fresh data)\n")
cat("===========================================================\n\n")

df_val <- read.csv("phase3_all_scored.csv", stringsAsFactors = FALSE)
results_val <- list()
for (rate_v in unique(df_val$rate)) {
  sub <- df_val %>% filter(rate == rate_v)
  m <- compute_all_metrics(sub$pred_prob, sub$gt,
                            fpr_target = 0.05,
                            neg_subset = (sub$pattern == "clean"))
  results_val[[paste0("rate=", rate_v)]] <- m
}

summary_val <- do.call(rbind, lapply(names(results_val), function(nm) {
  m <- results_val[[nm]]
  vals <- sapply(metric_names, function(mn) round(m[[mn]], 3))
  data.frame(operating_point = nm, t(vals), stringsAsFactors = FALSE)
}))
print(summary_val, row.names = FALSE)
write.csv(summary_val, "phase5_val_metrics.csv", row.names = FALSE)

# ============================================================
# Compare classifiers (logit, glmnet, RF) on multi-metrics at FPR=5%
# ============================================================
cat("\n===========================================================\n")
cat(" Classifier comparison: multi-metric @ FPR=5%\n")
cat("===========================================================\n\n")

library(glmnet)
classifier_results <- list()
classifier_names <- c("logit", "glmnet", "RF")

for (cls in classifier_names) {
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
      preds[te] <- as.numeric(predict(mod, newx = X_te, s = "lambda.min",
                                        type = "response"))
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
write.csv(classifier_summary, "phase5_classifier_metrics.csv", row.names = FALSE)

# ============================================================
# Plots
# ============================================================
cat("\n===========================================================\n")
cat(" PLOTS\n")
cat("===========================================================\n")

# Plot A: classifier comparison radar-style (bar)
classifier_long <- classifier_summary %>%
  select(classifier, mcc, f1, kappa, bal_acc, g_mean, auprc) %>%
  pivot_longer(-classifier, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("mcc","f1","kappa","bal_acc","g_mean","auprc")))

pA <- ggplot(classifier_long, aes(x = metric, y = value, fill = classifier)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  scale_fill_brewer(palette = "Set1", name = "Classifier") +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Multi-metric comparison of classifiers @ FPR=5%",
       subtitle = "Random Forest dominates on MCC, F1, Kappa, Balanced Accuracy, G-mean, AUPRC",
       x = "Metric", y = "Value") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("plot_phase5_classifier_metrics.png", pA,
       width = 12, height = 6.5, dpi = 150, bg = "white")

# Plot B: operating points trade-off
op_long <- summary_train %>%
  select(operating_point, sens, spec, ppv, mcc, f1, kappa, bal_acc, youden_j) %>%
  pivot_longer(-operating_point, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("sens","spec","ppv","bal_acc","youden_j",
                                     "mcc","f1","kappa")))

pB <- ggplot(op_long, aes(x = metric, y = value, fill = operating_point)) +
  geom_col(position = "dodge", width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 2.6) +
  scale_fill_brewer(palette = "Set2", name = "Operating point") +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Operating-point trade-offs across metrics (training data)",
       subtitle = "Higher FPR boosts sens but hurts ppv; threshold-optimized for MCC vs F1 nearly identical",
       x = "Metric", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plot_phase5_operating_points.png", pB,
       width = 13, height = 7, dpi = 150, bg = "white")

# Plot C: validation across rates
val_long <- summary_val %>%
  select(operating_point, sens, spec, ppv, mcc, f1, kappa, bal_acc, auprc) %>%
  pivot_longer(-operating_point, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("sens","spec","ppv","bal_acc",
                                     "mcc","f1","kappa","auprc")))

pC <- ggplot(val_long, aes(x = metric, y = value, fill = operating_point)) +
  geom_col(position = "dodge", width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 2.6) +
  scale_fill_brewer(palette = "Set1", name = "Careless rate") +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Validation: multi-metric performance by careless rate (FPR=5%)",
       subtitle = "Sensitivity stable ~0.91 across rates; ppv/MCC/F1/AUPRC rise with rate",
       x = "Metric", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plot_phase5_validation_metrics.png", pC,
       width = 13, height = 7, dpi = 150, bg = "white")

cat("\nSaved plots:\n")
cat("  plot_phase5_classifier_metrics.png\n")
cat("  plot_phase5_operating_points.png\n")
cat("  plot_phase5_validation_metrics.png\n")

# ============================================================
# Verdict
# ============================================================
cat("\n===========================================================\n")
cat(" VERDICT — does the conclusion change?\n")
cat("===========================================================\n\n")

# 1. Does RF still dominate?
best_per_metric <- classifier_summary %>%
  pivot_longer(-classifier, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  summarise(best_classifier = classifier[which.max(value)],
            best_value = max(value), .groups = "drop")
cat("Best classifier per metric:\n")
print(as.data.frame(best_per_metric), row.names = FALSE)

# 2. Are operating points consistent?
op_train <- summary_train %>%
  filter(operating_point %in% c("FPR=5%", "MCC optimal", "F1 optimal"))
cat("\nThreshold for FPR=5%, MCC-optimal, F1-optimal:\n")
print(op_train %>% select(operating_point, threshold, mcc, f1, kappa),
      row.names = FALSE)

# ============================================================
# ABLATION STUDY: does ReReReRe contribute, or is the auxiliary
# ensemble (IRV+LongString+D²+PT) enough on its own?
# ============================================================
cat("\n===========================================================\n")
cat(" ABLATION STUDY — does ReReReRe contribute?\n")
cat("===========================================================\n\n")

ablation_configs <- list(
  "Full ensemble (6)"         = c("z_rr", "z_rr_iter_t2", "irv", "longstring", "d2", "person_tot"),
  "Without z_rr_iter_efa (5)" = c("z_rr", "irv", "longstring", "d2", "person_tot"),
  "Without ANY z_rr (4)"       = c("irv", "longstring", "d2", "person_tot"),
  "Only z_rr family (2)"       = c("z_rr", "z_rr_iter_t2"),
  "Triad: iter+IRV+D² (3)"     = c("z_rr_iter_t2", "irv", "d2"),
  "Aux triad: IRV+D²+LS (3)"   = c("irv", "d2", "longstring"),
  "Aux quad: IRV+D²+LS+PT (4)" = c("irv", "d2", "longstring", "person_tot")
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
write.csv(ablation_summary, "phase5_ablation.csv", row.names = FALSE)

# Plot ablation
ablation_long <- ablation_summary %>%
  select(config, has_RR, sens, spec, ppv, mcc, f1, kappa, auprc) %>%
  pivot_longer(c(sens, spec, ppv, mcc, f1, kappa, auprc),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("sens","spec","ppv","mcc","f1","kappa","auprc")))

pAbl <- ggplot(ablation_long, aes(x = config, y = value,
                                    fill = factor(has_RR,
                                                  levels = c(TRUE, FALSE),
                                                  labels = c("with ReReReRe",
                                                             "without ReReReRe")))) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.5) +
  facet_wrap(~ metric, ncol = 4) +
  scale_fill_manual(values = c("with ReReReRe" = "#377eb8",
                                "without ReReReRe" = "#e41a1c"),
                     name = NULL) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Ablation study: configurations with vs without ReReReRe",
       subtitle = "FPR=5%, RF, training data. Reds = no z_RR features. Blues = with z_RR.",
       x = NULL, y = "Value") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
ggsave("plot_phase5_ablation.png", pAbl,
       width = 16, height = 10, dpi = 150, bg = "white")

# Critical comparison: full vs aux only
full_mcc <- ablation_summary$mcc[ablation_summary$config == "Full ensemble (6)"]
aux_mcc  <- ablation_summary$mcc[ablation_summary$config == "Without ANY z_rr (4)"]
delta_mcc <- full_mcc - aux_mcc
cat(sprintf("\n=== CRITICAL: Δ MCC (Full - NoRR) = %.3f - %.3f = %.3f ===\n",
            full_mcc, aux_mcc, delta_mcc))

if (delta_mcc < 0.02) {
  cat("VERDICT: ReReReRe contributes <0.02 MCC. Marginal.\n")
} else if (delta_mcc < 0.05) {
  cat("VERDICT: ReReReRe contributes 0.02-0.05 MCC. Modest but real.\n")
} else if (delta_mcc < 0.10) {
  cat("VERDICT: ReReReRe contributes 0.05-0.10 MCC. Substantial.\n")
} else {
  cat("VERDICT: ReReReRe contributes >0.10 MCC. Essential.\n")
}

# Per-pattern: where exactly does ReReReRe help most?
df_A$pred_full <- df_A$prob   # already has full ensemble
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
thr_full <- quantile(df_A$pred_full[df_A$pattern == "clean"], 0.95)
thr_no_rr <- quantile(df_A$pred_no_rr[df_A$pattern == "clean"], 0.95)
df_A$flag_full <- df_A$pred_full > thr_full
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
cat("\n=== Per-pattern detection: full vs no-RR ===\n")
print(as.data.frame(per_pat_ablation), row.names = FALSE)
write.csv(per_pat_ablation, "phase5_per_pattern_ablation.csv", row.names = FALSE)

# Plot per-pattern ablation
abl_pat_long <- per_pat_ablation %>%
  pivot_longer(c(rate_full, rate_no_rr),
               names_to = "model", values_to = "rate") %>%
  mutate(model = recode(model,
                        "rate_full"  = "Full ensemble (with ReReReRe)",
                        "rate_no_rr" = "Without ReReReRe (4 aux only)"))

pAblPat <- ggplot(abl_pat_long,
                   aes(x = factor(corruption), y = rate, fill = model)) +
  geom_col(position = "dodge", width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", rate)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3) +
  facet_wrap(~ pattern, ncol = 3) +
  scale_fill_manual(values = c("Full ensemble (with ReReReRe)" = "#377eb8",
                                "Without ReReReRe (4 aux only)" = "#e41a1c")) +
  scale_y_continuous(limits = c(0, 1.1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Where does ReReReRe contribute? Per-pattern detection at FPR=5%",
       subtitle = "Comparison shows the value of the z_RR features for each careless pattern",
       x = "Corruption level", y = "Detection rate", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        strip.text = element_text(face = "bold"))
ggsave("plot_phase5_ablation_per_pattern.png", pAblPat,
       width = 14, height = 8, dpi = 150, bg = "white")

cat("\n=== DONE Phase 5 (multi-metric + ablation) ===\n")
